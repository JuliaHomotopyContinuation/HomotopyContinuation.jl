export InterpretedHomotopy, InterpretedHomotopy


"""
    InterpretedHomotopy <: AbstractHomotopy

An [`AbstractHomotopy`](@ref) which automatically generates a program for the
fast evaluation of `H` and its Jacobian. The program is however, not compiled
but rather interpreted. See also [`CompiledHomotopy`](@ref).

    InterpretedHomotopy(H::Homotopy)

Construct an `InterpretedHomotopy` from the given [`Homotopy`](@ref) `H`.
"""
mutable struct InterpretedHomotopy <: AbstractHomotopy
    homotopy::Homotopy
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

function InterpretedHomotopy(H::Homotopy)
    eval_ComplexF64 = interpreter(ComplexF64, H)
    eval_ComplexDF64 = interpreter(ComplexDF64, eval_ComplexF64)
    eval_acb = nothing
    taylor_ComplexF64 = TaylorInterpreters{ComplexF64}()
    jac_ComplexF64 = jacobian_interpreter(ComplexF64, H)
    jac_acb = nothing

    InterpretedHomotopy(
        H,
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

Base.size(H::InterpretedHomotopy) = size(H.homotopy)
variables(H::InterpretedHomotopy) = variables(H.homotopy)
parameters(H::InterpretedHomotopy) = parameters(H.homotopy)
variable_groups(H::InterpretedHomotopy) = variable_groups(H.homotopy)
Base.:(==)(H::InterpretedHomotopy, G::InterpretedHomotopy) = H.homotopy == G.homotopy

function Base.show(io::IO, H::InterpretedHomotopy)
    print(io, "Interpreted: ")
    show(io, H.homotopy)
end

(H::InterpretedHomotopy)(x, t, p = nothing) = H.homotopy(x, t, p)
function evaluate!(u, H::InterpretedHomotopy, x, t, p = nothing)
    T = promote_type(_interpreter_precision(x, p), _interpreter_precision(t))
    execute!(u, _evaluation_interpreter(H, T), x, t, p)
end


function evaluate_and_jacobian!(u, U, H::InterpretedHomotopy, x, t, p = nothing)
    T = promote_type(_interpreter_precision(x, p), _interpreter_precision(t))
    execute!(u, U, _jacobian_interpreter!(H, T), x, t, p)
    nothing
end
function jacobian!(U, H::InterpretedHomotopy, x, t, p = nothing)
    T = promote_type(_interpreter_precision(x, p), _interpreter_precision(t))
    execute!(nothing, U, _jacobian_interpreter!(H, T), x, t, p)
    U
end

function taylor!(
    u::AbstractVector,
    Order::Val{M},
    H::InterpretedHomotopy,
    x,
    t_,
    p::Union{Nothing,AbstractVector} = nothing;
    assign_highest_order_only = u isa Vector,
) where {M}
    t = (t_, 1)
    T = promote_type(_interpreter_precision(x, p), _interpreter_precision(t_))
    I = _get_or_create_taylor_interpreter!(H, Order, T)
    execute_taylor!(
        u,
        Order,
        I,
        x,
        t,
        p;
        assign_highest_order_only = assign_highest_order_only,
    )
    u
end

# Acb
function taylor!(u, v::Val{K}, H::InterpretedHomotopy, x, t, p, incremental::Bool) where {K}
    K > 3 && return taylor!(u, v, H, x, t, p)
    T = promote_type(_interpreter_precision(x, p), _interpreter_precision(t))
    I = _predictor_taylor_interpreter!(H, T)
    predictor_taylor!(u, v, I, x, (t, 1), p, incremental)
end
taylor!(u, v::Val, H::InterpretedHomotopy, x, t, incremental::Bool) =
    taylor!(u, v, H, x, t, nothing, incremental)

function evaluate!(
    u::AbstractArray{<:Union{Acb,AcbRef}},
    H::InterpretedHomotopy,
    x::AbstractArray{<:Union{Acb,AcbRef}},
    t,
    p = nothing;
    prec = max(precision(first(u)), precision(first(x))),
)
    if isnothing(H.eval_acb)
        H.eval_acb = interpreter(AcbRefVector, H.eval_ComplexF64)
    end
    setprecision!(H.eval_acb, prec)
    execute!(u, H.eval_acb, x, t, p)
end

function evaluate_and_jacobian!(
    u,
    U,
    H::InterpretedHomotopy,
    x::AbstractArray{<:Union{Acb,AcbRef}},
    t,
    p = nothing;
    prec = max(precision(first(u)), precision(first(U)), precision(first(x))),
)
    if isnothing(H.jac_acb)
        H.jac_acb = interpreter(AcbRefVector, H.jac_ComplexF64)
    end
    setprecision!(H.jac_acb, prec)

    execute!(u, U, H.jac_acb, x, t, p)
    nothing
end


function (H::InterpretedHomotopy)(x::AbstractArray{<:Union{Acb,AcbRef}}, t, p = nothing)
    u = Arblib.AcbVector(first(size(H)); prec = precision(first(x)))
    evaluate!(u, H, x, t, p)
    u
end
(H::Homotopy)(x::AbstractVector{<:Union{Acb,AcbRef}}, t, p::Nothing = nothing) =
    InterpretedHomotopy(H)(x, t, p)
(H::Homotopy)(x::AbstractVector{<:Union{Acb,AcbRef}}, t, p::AbstractArray) =
    InterpretedHomotopy(H)(x, t, p)
(H::Homotopy)(x::AbstractMatrix{<:Union{Acb,AcbRef}}, t, p = nothing) =
    InterpretedHomotopy(H)(x, t, p)
