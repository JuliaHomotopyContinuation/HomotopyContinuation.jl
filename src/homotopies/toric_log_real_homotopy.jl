export ToricLogRealHomotopty

Base.@kwdef struct ToricLogRealHomotopty{S<:AbstractSystem} <: AbstractHomotopy
    H::GeneralParameterHomotopy{S,ToricLogRealEmbedding}
end

function ToricLogRealHomotopty(
    system::AbstractSystem,
    system_coeffs::Vector{Vector{ComplexF64}},
    system_support::AbstractVector{<:AbstractMatrix},
)
    m = ModelKit.nparameters(system)
    m1 = sum(length, system_coeffs)
    m == m1 || throw(
        ArgumentError(
            "System parameters and coefficients do not have the same size, got $m and $m1",
        ),
    )
    γ = ToricLogRealEmbedding(reduce(vcat, system_coeffs), system_support)

    H = GeneralParameterHomotopy(system, γ)
    ToricLogRealHomotopty(H)
end

Base.size(H::ToricLogRealHomotopty) = size(H.H)

function update!(
    H::ToricLogRealHomotopty,
    lifting::AbstractVector{<:AbstractVector},
    cell::MixedSubdivisions.MixedCell;
    min_weight::Union{Nothing,Float64} = nothing,
    max_weight::Union{Nothing,Float64} = nothing,
)
    update!(γ(H.H), lifting, cell, min_weight = min_weight, max_weight = max_weight)
end

function ModelKit.evaluate!(
    u,
    H::ToricLogRealHomotopty,
    x::AbstractVector,
    t,
    p::Nothing = nothing,
)
    evaluate!(u, H.H, x, t, p)
end

function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::ToricLogRealHomotopty,
    x::AbstractVector,
    t,
    p::Nothing = nothing,
)
    evaluate_and_jacobian!(u, U, H.H, x, t, p)
    nothing
end

function ModelKit.taylor!(u, v, H::ToricLogRealHomotopty, x, t, p::Nothing = nothing)
    taylor!(u, v, H.H, x, t)
    nothing
end
