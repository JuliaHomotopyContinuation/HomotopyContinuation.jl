export InterpretedSystem

mutable struct TaylorInterpreters{T}
    order_1::Union{Nothing,Interpreter{Vector{TruncatedTaylorSeries{2,T}}}}
    order_2::Union{Nothing,Interpreter{Vector{TruncatedTaylorSeries{3,T}}}}
    order_3::Union{Nothing,Interpreter{Vector{TruncatedTaylorSeries{4,T}}}}
    order_4::Union{Nothing,Interpreter{Vector{TruncatedTaylorSeries{5,T}}}}
    predictor::Union{Nothing,IncrementalTaylorInterpreter{T}}
end
TaylorInterpreters{T}() where {T} =
    TaylorInterpreters{T}(nothing, nothing, nothing, nothing, nothing)

"""
    InterpretedSystem <: AbstractSystem

An [`AbstractSystem`](@ref) which automatically generates a program for the
fast evaluation of `F` and its Jacobian. The program is however, not compiled
but rather interpreted. See also [`CompiledSystem`](@ref).

    InterpretedSystem(F::System)

Construct an `InterpretedSystem` from the given [`System`](@ref) `F`.
"""
mutable struct InterpretedSystem <: AbstractSystem
    system::System
    eval_ComplexF64::Interpreter{Vector{ComplexF64}}
    eval_ComplexDF64::Interpreter{Vector{ComplexDF64}}
    # Construct acb lazily since its often not needed
    eval_acb::Union{Nothing,Interpreter{AcbRefVector}}
    taylor_ComplexF64::TaylorInterpreters{ComplexF64}
    taylor_ComplexDF64::TaylorInterpreters{ComplexDF64}
    jac_ComplexF64::Interpreter{Vector{ComplexF64}}
    jac_ComplexDF64::Union{Nothing,Interpreter{Vector{ComplexDF64}}}
    jac_acb::Union{Nothing,Interpreter{AcbRefVector}}
end

function InterpretedSystem(F::System)
    eval_ComplexF64 = interpreter(ComplexF64, F)
    eval_ComplexDF64 = interpreter(ComplexDF64, eval_ComplexF64)
    eval_acb = nothing
    taylor_ComplexF64 = TaylorInterpreters{ComplexF64}()
    jac_ComplexF64 = jacobian_interpreter(ComplexF64, F)
    jac_acb = nothing

    InterpretedSystem(
        F,
        eval_ComplexF64,
        eval_ComplexDF64,
        eval_acb,
        taylor_ComplexF64,
        TaylorInterpreters{ComplexDF64}(),
        jac_ComplexF64,
        nothing,
        jac_acb,
    )
end

Base.size(F::InterpretedSystem) = size(F.system)
variables(F::InterpretedSystem) = variables(F.system)
parameters(F::InterpretedSystem) = parameters(F.system)
variable_groups(F::InterpretedSystem) = variable_groups(F.system)
Base.:(==)(F::InterpretedSystem, G::InterpretedSystem) = F.system == G.system

function Base.show(io::IO, F::InterpretedSystem)
    print(io, "Interpreted: ")
    show(io, F.system)
end

(F::InterpretedSystem)(x, p = nothing) = F.system(x, p)

# Choose arithmetic from the inputs, not the output buffer: mixed-precision
# refinement intentionally stores an accurate residual in a Float64 buffer.
_interpreter_precision(::Type) = ComplexF64
_interpreter_precision(::Type{DoubleF64}) = ComplexDF64
_interpreter_precision(::Type{ComplexDF64}) = ComplexDF64
_interpreter_precision(::Type{TruncatedTaylorSeries{N,T}}) where {N,T} =
    _interpreter_precision(T)
_interpreter_precision(::Type{NTuple{N,T}}) where {N,T} = _interpreter_precision(T)
_interpreter_precision(x::AbstractArray{T}) where {T} = _interpreter_precision(T)
_interpreter_precision(x) = _interpreter_precision(typeof(x))
_interpreter_precision(x, p) =
    promote_type(_interpreter_precision(x), _interpreter_precision(p))

_evaluation_interpreter(obj, ::Type{ComplexF64}) = obj.eval_ComplexF64
_evaluation_interpreter(obj, ::Type{ComplexDF64}) = obj.eval_ComplexDF64
_taylor_interpreters(obj, ::Type{ComplexF64}) = obj.taylor_ComplexF64
_taylor_interpreters(obj, ::Type{ComplexDF64}) = obj.taylor_ComplexDF64
function _predictor_taylor_interpreter!(obj, ::Type{T}) where {T}
    caches = _taylor_interpreters(obj, T)
    I = caches.predictor
    if I === nothing
        I = IncrementalTaylorInterpreter(T, _evaluation_interpreter(obj, T))
        caches.predictor = I
    end
    I
end

predictor_taylor!(u, v, F::AbstractSystem, x, p, incremental::Bool) = taylor!(u, v, F, x, p)
function predictor_taylor!(
    u,
    v::Val{K},
    F::InterpretedSystem,
    x,
    p,
    incremental::Bool,
) where {K}
    K > 3 && return taylor!(u, v, F, x, p)
    I = _predictor_taylor_interpreter!(F, _interpreter_precision(x, p))
    predictor_taylor!(u, v, I, x, nothing, p, incremental)
end
_jacobian_interpreter!(obj, ::Type{ComplexF64}) = obj.jac_ComplexF64
function _jacobian_interpreter!(obj, ::Type{ComplexDF64})
    I = obj.jac_ComplexDF64
    if isnothing(I)
        I = interpreter(ComplexDF64, obj.jac_ComplexF64)
        obj.jac_ComplexDF64 = I
    end
    I
end

function evaluate!(u, F::InterpretedSystem, x, p = nothing)
    execute!(u, _evaluation_interpreter(F, _interpreter_precision(x, p)), x, p)
end


function evaluate_and_jacobian!(u, U, F::InterpretedSystem, x, p = nothing)
    execute!(u, U, _jacobian_interpreter!(F, _interpreter_precision(x, p)), x, p)
    nothing
end
function jacobian!(U, F::InterpretedSystem, x, p = nothing)
    execute!(nothing, U, _jacobian_interpreter!(F, _interpreter_precision(x, p)), x, p)
    U
end

# Helper: lazily create and cache the Taylor interpreter for order M.
# @generated so that field access (order_1..order_4) and the TruncatedTaylorSeries
# arity (M+1) are resolved at compile time, giving a type-stable return value.
@generated function _get_or_create_taylor_interpreter!(
    obj,
    ::Val{M},
    ::Type{T} = ComplexF64,
) where {M,T}
    1 <= M <= 4 || error("Taylor order M must be between 1 and 4, got $M")
    field = [:order_1, :order_2, :order_3, :order_4][M]
    N = M + 1
    cache = T === ComplexDF64 ? :taylor_ComplexDF64 : :taylor_ComplexF64
    evaluator = T === ComplexDF64 ? :eval_ComplexDF64 : :eval_ComplexF64
    quote
        I = obj.$cache.$field
        if isnothing(I)
            I′ = interpreter(TruncatedTaylorSeries{$N,$T}, obj.$evaluator)
            obj.$cache.$field = I′
            I′
        else
            I
        end
    end
end

function taylor!(
    u::AbstractVector,
    Order::Val{M},
    F::InterpretedSystem,
    x,
    p = nothing;
    assign_highest_order_only::Bool = u isa Vector,
) where {M}
    I = _get_or_create_taylor_interpreter!(F, Order, _interpreter_precision(x, p))
    execute_taylor!(
        u,
        Order,
        I,
        x,
        p;
        assign_highest_order_only = assign_highest_order_only,
    )
    u
end

# Acb
function evaluate!(
    u::AbstractArray{<:Union{Acb,AcbRef}},
    F::InterpretedSystem,
    x::AbstractArray{<:Union{Acb,AcbRef}},
    p = nothing;
    prec = max(precision(first(u)), precision(first(x))),
)
    if isnothing(F.eval_acb)
        F.eval_acb = interpreter(AcbRefVector, F.eval_ComplexF64)
    end
    setprecision!(F.eval_acb, prec)
    execute!(u, F.eval_acb, x, p)
end

function evaluate_and_jacobian!(
    u,
    U,
    F::InterpretedSystem,
    x::AbstractArray{<:Union{Acb,AcbRef}},
    p = nothing;
    prec = max(precision(first(u)), precision(first(U)), precision(first(x))),
)
    if isnothing(F.jac_acb)
        F.jac_acb = interpreter(AcbRefVector, F.jac_ComplexF64)
    end
    setprecision!(F.jac_acb, prec)

    execute!(u, U, F.jac_acb, x, p)
    nothing
end


function (F::InterpretedSystem)(x::AbstractArray{<:Union{Acb,AcbRef}}, p = nothing)
    u = Arblib.AcbVector(first(size(F)); prec = precision(first(x)))
    evaluate!(u, F, x, p)
    u
end
(F::System)(x::AbstractVector{<:Union{Acb,AcbRef}}, p::Nothing = nothing) =
    InterpretedSystem(F)(x, p)
(F::System)(x::AbstractVector{<:Union{Acb,AcbRef}}, p::AbstractArray) =
    InterpretedSystem(F)(x, p)
(F::System)(x::AbstractMatrix{<:Union{Acb,AcbRef}}, p = nothing) =
    InterpretedSystem(F)(x, p)
