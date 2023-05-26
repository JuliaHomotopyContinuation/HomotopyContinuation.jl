export ToricHomotopy

Base.@kwdef struct ToricHomotopy{S<:AbstractSystem} <: AbstractHomotopy
    H::GeneralParameterHomotopy{S,ToricTEmbedding}
end

function ToricHomotopy(
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
    γ = ToricTEmbedding(reduce(vcat, system_coeffs), system_support)

    H = GeneralParameterHomotopy(system, γ)
    ToricHomotopy(H)
end

Base.size(H::ToricHomotopy) = size(H.H)

function update!(
    H::ToricHomotopy,
    lifting::AbstractVector{<:AbstractVector},
    cell::MixedSubdivisions.MixedCell;
    min_weight::Union{Nothing,Float64} = nothing,
    max_weight::Union{Nothing,Float64} = nothing,
)
    update!(γ(H.H), lifting, cell, min_weight = min_weight, max_weight = max_weight)
end

function ModelKit.evaluate!(u, H::ToricHomotopy, x::AbstractVector, t, p::Nothing = nothing)
    evaluate!(u, H.H, x, t, p)
end

function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::ToricHomotopy,
    x::AbstractVector,
    t,
    p::Nothing = nothing,
)
    evaluate_and_jacobian!(u, U, H.H, x, t, p)
    nothing
end

function ModelKit.taylor!(u, v, H::ToricHomotopy, x, t, p::Nothing = nothing)
    taylor!(u, v, H.H, x, t)
    nothing
end
