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
        jac_ComplexF64,
        jac_acb,
    )
end

Base.size(H::InterpretedHomotopy) = size(H.homotopy)
variables(H::InterpretedHomotopy) = variables(H.homotopy)
parameters(H::InterpretedHomotopy) = parameters(H.homotopy)
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

@generated function taylor!(
    u::AbstractVector,
    Order::Val{M},
    H::InterpretedHomotopy,
    x,
    t_,
    p::Union{Nothing,AbstractVector} = nothing;
    assign_highest_order_only = u isa Vector,
) where {M}
    if M == 1
        quote
            t = (t_, 1)
            I = H.taylor_ComplexF64.order_1
            if isnothing(I)
                I′ = interpreter(TruncatedTaylorSeries{2,ComplexF64}, H.eval_ComplexF64)
                H.taylor_ComplexF64.order_1 = I′
                execute_taylor!(
                    u,
                    Order,
                    I′,
                    x,
                    t,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            else
                execute_taylor!(
                    u,
                    Order,
                    I,
                    x,
                    t,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            end
            u
        end
    elseif M == 2
        quote
            t = (t_, 1)
            I = H.taylor_ComplexF64.order_2
            if isnothing(I)
                I′ = interpreter(TruncatedTaylorSeries{3,ComplexF64}, H.eval_ComplexF64)
                H.taylor_ComplexF64.order_2 = I′
                execute_taylor!(
                    u,
                    Order,
                    I′,
                    x,
                    t,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            else
                execute_taylor!(
                    u,
                    Order,
                    I,
                    x,
                    t,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            end
            u
        end
    elseif M == 3
        quote
            t = (t_, 1)
            I = H.taylor_ComplexF64.order_3
            if isnothing(I)
                I′ = interpreter(TruncatedTaylorSeries{4,ComplexF64}, H.eval_ComplexF64)
                H.taylor_ComplexF64.order_3 = I′
                execute_taylor!(
                    u,
                    Order,
                    I′,
                    x,
                    t,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            else
                execute_taylor!(
                    u,
                    Order,
                    I,
                    x,
                    t,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            end
            u
        end
    elseif M == 4
        quote
            t = (t_, 1)
            I = H.taylor_ComplexF64.order_4
            if isnothing(I)
                I′ = interpreter(TruncatedTaylorSeries{5,ComplexF64}, H.eval_ComplexF64)
                H.taylor_ComplexF64.order_4 = I′
                execute_taylor!(
                    u,
                    Order,
                    I′,
                    x,
                    t,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            else
                execute_taylor!(
                    u,
                    Order,
                    I,
                    x,
                    t,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            end
            u
        end
    end

end

# Acb
function evaluate!(
    u::AbstractArray{<:Union{Acb,AcbRef}},
    H::InterpretedHomotopy,
    x::AbstractArray{<:Union{Acb,AcbRef}},
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
