export ParameterHomotopy

"""
    ParameterHomotopy(F::Union{AbstractSystem,System}; start_parameters, target_parameters)
    ParameterHomotopy(F::Union{AbstractSystem,System}, start_parameters, target_parameters)

Construct the parameter homotopy ``H(x,t) = F(x; t p + (1 - t) q)`` where ``p`` is
`start_parameters` and ``q`` is `target_parameters`.
"""
struct ParameterHomotopy{T<:AbstractSystem} <: AbstractHomotopy
    F::T
    p::Vector{ComplexF64}
    q::Vector{ComplexF64}
    general_homotopy::GeneralParameterHomotopy{T,StraightLineTEmbedding}
end

function ParameterHomotopy(
    F;
    start_parameters::AbstractVector,
    target_parameters::AbstractVector,
)
    ParameterHomotopy(F, start_parameters, target_parameters)
end
function ParameterHomotopy(
    F::ModelKit.System,
    p::AbstractVector,
    q::AbstractVector;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
)
    ParameterHomotopy(fixed(F; compile = compile), p, q)
end
function ParameterHomotopy(F::AbstractSystem, p::AbstractVector, q::AbstractVector)
    @assert length(p) == length(q) == nparameters(F)

    p̂ = Vector{ComplexF64}(p)
    q̂ = Vector{ComplexF64}(q)
    H = GeneralParameterHomotopy(F, StraightLineTEmbedding(p̂, q̂))

    ParameterHomotopy(F, p̂, q̂, H)
end

Base.size(H::ParameterHomotopy) = size(H.F)

function start_parameters!(H::ParameterHomotopy, p)
    update_γ!(H, p, γ(H).q)
    H
end
function target_parameters!(H::ParameterHomotopy, q)
    update_γ!(H, γ(H).p, q)
    H
end
function parameters!(H::ParameterHomotopy, p, q)
    update_γ!(H, p, q)
    H
end

function ModelKit.evaluate!(u, H::ParameterHomotopy, x, t, p::Nothing = nothing)
    evaluate!(u, H.general_homotopy, x, t)
end

function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::ParameterHomotopy,
    x,
    t,
    p::Nothing = nothing,
)
    evaluate_and_jacobian!(u, U, H.general_homotopy, x, t)
end

function ModelKit.taylor!(u, v::Val, H::ParameterHomotopy, tx, t, p::Nothing = nothing)
    taylor!(u, v, H.general_homotopy, tx, t)
    u
end
