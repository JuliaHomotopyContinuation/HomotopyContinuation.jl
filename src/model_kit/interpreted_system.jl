export InterpretedSystem

mutable struct TaylorInterpreters{T}
    order_1::Union{Nothing,Interpreter{Vector{TruncatedTaylorSeries{2,T}}}}
    order_2::Union{Nothing,Interpreter{Vector{TruncatedTaylorSeries{3,T}}}}
    order_3::Union{Nothing,Interpreter{Vector{TruncatedTaylorSeries{4,T}}}}
    order_4::Union{Nothing,Interpreter{Vector{TruncatedTaylorSeries{5,T}}}}
end
TaylorInterpreters{T}() where {T} =
    TaylorInterpreters{T}(nothing, nothing, nothing, nothing)

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
    jac_ComplexF64::Interpreter{Vector{ComplexF64}}
    jac_acb::Union{Nothing,Interpreter{AcbRefVector}}
end

function InterpretedSystem(F::System; kwargs...)
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
        jac_ComplexF64,
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
function evaluate!(u, F::InterpretedSystem, x, p = nothing)
    execute!(u, F.eval_ComplexF64, x, p)
end
function evaluate!(u, F::InterpretedSystem, x::AbstractVector{ComplexDF64}, p = nothing)
    execute!(u, F.eval_ComplexDF64, x, p)
end


function evaluate_and_jacobian!(u, U, F::InterpretedSystem, x, p = nothing)
    execute!(u, U, F.jac_ComplexF64, x, p)
    nothing
end
function jacobian!(U, F::InterpretedSystem, x, p = nothing)
    execute!(nothing, U, F.jac_ComplexF64, x, p;)
    U
end
jacobian!(U, F::InterpretedSystem, x, p, cache) = jacobian!(U, F, x, p)

# Helper: lazily create and cache the Taylor interpreter for order M.
# @generated so that field access (order_1..order_4) and the TruncatedTaylorSeries
# arity (M+1) are resolved at compile time, giving a type-stable return value.
@generated function _get_or_create_taylor_interpreter!(obj, ::Val{M}) where {M}
    1 <= M <= 4 || error("Taylor order M must be between 1 and 4, got $M")
    field = [:order_1, :order_2, :order_3, :order_4][M]
    N = M + 1
    quote
        I = obj.taylor_ComplexF64.$field
        if isnothing(I)
            I′ = interpreter(TruncatedTaylorSeries{$N,ComplexF64}, obj.eval_ComplexF64)
            obj.taylor_ComplexF64.$field = I′
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
    I = _get_or_create_taylor_interpreter!(F, Order)
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
    I = F.eval_acb::Interpreter{AcbRefVector}
    setprecision!(I, prec)
    execute!(u, I, x, p)
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
    I = F.jac_acb::Interpreter{AcbRefVector}
    setprecision!(I, prec)

    execute!(u, U, I, x, p)
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
