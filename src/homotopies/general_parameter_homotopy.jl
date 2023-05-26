export GeneralParameterHomotopy

include("general_parameter_homotopy/t_embedding.jl")

"""
    GeneralParameterHomotopy(F::Union{AbstractSystem,System}; start_parameters, target_parameters)
    GeneralParameterHomotopy(F::Union{AbstractSystem,System}, start_parameters, target_parameters)

Construct the parameter homotopy ``H(x,t) = F(x; t p + (1 - t) q)`` where ``p`` is
`start_parameters` and ``q`` is `target_parameters`.
"""
struct GeneralParameterHomotopy{S<:AbstractSystem,T<:TEmbedding} <: AbstractHomotopy
    F::S
    γ::CacheTEmbedding{T}
    γt::Vector{ComplexF64}
    taylor_γt::TaylorVector{5,ComplexF64}
end

function GeneralParameterHomotopy(
    F::ModelKit.System,
    γ::TEmbedding;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
)
    GeneralParameterHomotopy(fixed(F; compile = compile), γ)
end
function GeneralParameterHomotopy(F::AbstractSystem, γ::TEmbedding)
    γt = zeros(ComplexF64, nparameters(F))
    taylor_γt = TaylorVector{5}(ComplexF64, nparameters(F))

    GeneralParameterHomotopy(
        F,
        γ isa CacheTEmbedding ? γ : CacheTEmbedding(γ),
        γt,
        taylor_γt,
    )
end

γ(H::GeneralParameterHomotopy) = H.γ.embedding
update_γ!(H::GeneralParameterHomotopy, args...) = update!(H.γ, args...)
Base.size(H::GeneralParameterHomotopy) = size(H.F)

function ModelKit.evaluate!(u, H::GeneralParameterHomotopy, x, t, p::Nothing = nothing)
    evaluate!(H.γt, H.γ, t)
    evaluate!(u, H.F, x, H.γt)
end

function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::GeneralParameterHomotopy,
    x,
    t,
    p::Nothing = nothing,
)
    evaluate!(H.γt, H.γ, t)
    evaluate_and_jacobian!(u, U, H.F, x, H.γt)
end

function ModelKit.taylor!(
    u,
    v::Val,
    H::GeneralParameterHomotopy,
    tx,
    tṫ,
    p::Nothing = nothing,
)
    taylor!(H.taylor_γt, H.γ, tṫ)
    taylor!(u, v, H.F, tx, H.taylor_γt)
    u
end
