export InterpretedSystem

mutable struct TaylorInterpreters{T}
    order_1::Union{Nothing,Interpreter{Vector{TruncatedTaylorSeries{2,T}}}}
    order_2::Union{Nothing,Interpreter{Vector{TruncatedTaylorSeries{3,T}}}}
    order_3::Union{Nothing,Interpreter{Vector{TruncatedTaylorSeries{4,T}}}}
    order_4::Union{Nothing,Interpreter{Vector{TruncatedTaylorSeries{5,T}}}}
end
TaylorInterpreters{T}() where {T} =
    TaylorInterpreters{ComplexF64}(nothing, nothing, nothing, nothing)

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
        jac_ComplexF64,
        jac_acb,
    )
end

Base.size(F::InterpretedSystem) = size(F.system)
variables(F::InterpretedSystem) = variables(F.system)
parameters(F::InterpretedSystem) = parameters(F.system)
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

@generated function taylor!(
    u::AbstractVector,
    Order::Val{M},
    F::InterpretedSystem,
    x,
    p = nothing;
    assign_highest_order_only::Bool = u isa Vector,
) where {M}
    if M == 1
        quote
            I = F.taylor_ComplexF64.order_1
            if isnothing(I)
                I′ = interpreter(TruncatedTaylorSeries{2,ComplexF64}, F.eval_ComplexF64)
                F.taylor_ComplexF64.order_1 = I′
                execute_taylor!(
                    u,
                    Order,
                    I′,
                    x,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            else
                execute_taylor!(
                    u,
                    Order,
                    I,
                    x,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            end
            u
        end
    elseif M == 2
        quote
            I = F.taylor_ComplexF64.order_2
            if isnothing(I)
                I′ = interpreter(TruncatedTaylorSeries{3,ComplexF64}, F.eval_ComplexF64)
                F.taylor_ComplexF64.order_2 = I′
                execute_taylor!(
                    u,
                    Order,
                    I′,
                    x,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            else
                execute_taylor!(
                    u,
                    Order,
                    I,
                    x,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            end
            u
        end
    elseif M == 3
        quote
            I = F.taylor_ComplexF64.order_3
            if isnothing(I)
                I′ = interpreter(TruncatedTaylorSeries{4,ComplexF64}, F.eval_ComplexF64)
                F.taylor_ComplexF64.order_3 = I′
                execute_taylor!(
                    u,
                    Order,
                    I′,
                    x,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            else
                execute_taylor!(
                    u,
                    Order,
                    I,
                    x,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            end
            u
        end
    elseif M == 4
        quote
            I = F.taylor_ComplexF64.order_4
            if isnothing(I)
                I′ = interpreter(TruncatedTaylorSeries{5,ComplexF64}, F.eval_ComplexF64)
                F.taylor_ComplexF64.order_4 = I′
                execute_taylor!(
                    u,
                    Order,
                    I′,
                    x,
                    p;
                    assign_highest_order_only = assign_highest_order_only,
                )
            else
                execute_taylor!(
                    u,
                    Order,
                    I,
                    x,
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
