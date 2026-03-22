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
    jac_ComplexF64::Interpreter{Vector{ComplexF64}}
    jac_acb::Union{Nothing,Interpreter{AcbRefVector}}
end

function InterpretedHomotopy(H::Homotopy; kwargs...)
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
        jac_ComplexF64,
        jac_acb,
    )
end

Base.size(H::InterpretedHomotopy) = size(H.homotopy)
variables(H::InterpretedHomotopy) = variables(H.homotopy)
parameters(H::InterpretedHomotopy) = parameters(H.homotopy)
variable_groups(H::InterpretedHomotopy) = variable_groups(H.homotopy)
Homotopy(H::InterpretedHomotopy) = H.homotopy
Base.:(==)(H::InterpretedHomotopy, G::InterpretedHomotopy) = H.homotopy == G.homotopy

function Base.show(io::IO, H::InterpretedHomotopy)
    print(io, "Interpreted: ")
    show(io, H.homotopy)
end

(H::InterpretedHomotopy)(x, t, p = nothing) = H.homotopy(x, t, p)
function evaluate!(u, H::InterpretedHomotopy, x, t, p = nothing)
    execute!(u, H.eval_ComplexF64, x, t, p)
end
function evaluate!(
    u,
    H::InterpretedHomotopy,
    x::AbstractVector{ComplexDF64},
    t,
    p = nothing,
)
    execute!(u, H.eval_ComplexDF64, x, t, p)
end


function evaluate_and_jacobian!(u, U, H::InterpretedHomotopy, x, t, p = nothing)
    execute!(u, U, H.jac_ComplexF64, x, t, p)
    nothing
end
function jacobian!(U, H::InterpretedHomotopy, x, t, p = nothing)
    execute!(nothing, U, H.jac_ComplexF64, x, t, p;)
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
    I = _get_or_create_taylor_interpreter!(H, Order)
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
    I = H.eval_acb::Interpreter{AcbRefVector}
    setprecision!(I, prec)
    execute!(u, I, x, t, p)
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
    I = H.jac_acb::Interpreter{AcbRefVector}
    setprecision!(I, prec)

    execute!(u, U, I, x, t, p)
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
