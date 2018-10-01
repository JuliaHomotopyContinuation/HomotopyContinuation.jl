export StaticParameterHomotopy, ParameterHomotopyCache, gamma

import MultivariatePolynomials
const MP = MultivariatePolynomials
import StaticPolynomials
const SP = StaticPolynomials
import StaticArrays: MVector, SVector

import ..Utilities

"""
    ParameterHomotopy(F, parameters;
        p₁=randn(ComplexF64, length(parameters)),
        p₀=randn(ComplexF64, length(parameters)),
        γ₁=nothing, γ₀=nothing)

Construct the homotopy
```math
H(x, t) = F(x, (tγ₁p₁+(1-t)γ₀p₀) / (tγ₁+(1-t)γ₀))
```,
where `p₁` and `p₀` are a vector of parameter values for ``F`` and
`γ₁` and `γ₀` are complex numbers. If `γ₁` or `γ₀` is `nothing`, it is assumed
that `γ₁` and `γ₀` are ``1``.
The input ``parameters`` specifies the parameter variables of ``F``.
Neccessarily, ``length(parameters) == length(p₁) == length(p₀)``.

Note that `p₁` and `p₀` are stored as a tuple `p` of `SVectors` and `γ₁` and `γ₀`
are stored as a tuple `γ` or as `γ=nothing`

    ParameterHomotopy(F, parameters;
        start_parameters=randn(ComplexF64, length(parameters)),
        target_parameters=randn(ComplexF64, length(parameters)),
        start_gamma=nothing, target_gamma=nothing)

This is a non-unicode variant. `γ₁=start_parameters`, `γ₀=target_parameters`,
`γ₁=start_gamma`, γ₀=`target_gamma`.
"""
mutable struct ParameterHomotopy{N, NVars, NParams, T<:Number, PolyTuple} <: AbstractHomotopy
    F::SP.PolynomialSystem{N, NVars, NParams, PolyTuple}
    p::NTuple{2, SVector{NParams, T}}
    γ::Union{Nothing, NTuple{2, ComplexF64}}
end

function ParameterHomotopy(F::SP.PolynomialSystem{N, NVars, NParams, T}, p₁, p₀;
    start_gamma=nothing, γ₀ = start_gamma,
    target_gamma=nothing, γ₁ = target_gamma
    ) where {N, NVars, NParams, T}
    @assert length(p₁) == length(p₀) == NParams

    if start_gamma === nothing || target_gamma === nothing
        γ = nothing
    else
        γ = (γ₁, γ₀)
    end

    ParameterHomotopy(F, promote(MVector{NParams}(p₁), MVector{NParams}(p₀)), γ)
end

function ParameterHomotopy(F::Vector{T}, parameters::AbstractVector{V};
    start_parameters=randn(ComplexF64, length(parameters)), p₁ = start_parameters,
    target_parameters=randn(ComplexF64, length(parameters)), p₀ = target_parameters,
    kwargs...) where {T<:MP.AbstractPolynomialLike, V<:MP.AbstractVariable}
    G = SP.PolynomialSystem(F; parameters=parameters)
    ParameterHomotopy(G, p₁, p₀; kwargs...)
end

(H::ParameterHomotopy)(x, t, c=NullCache()) = evaluate(H, x, t, c)

cache(::ParameterHomotopy, x, t) = NullCache()

Base.size(H::ParameterHomotopy{N, NVars}) where {N, NVars} = (N, NVars)

p₁(H::ParameterHomotopy) = SVector(H.p₁)
p₀(H::ParameterHomotopy) = SVector(H.p₀)

@inline function p(H::ParameterHomotopy, t)
    if H.γ === nothing
        return t * p₁(H) + (1 - t) * p₀(H))
    end
    ⁠γ₁, ⁠γ₀ = H.γ
    ₜγ₁, ₁₋ₜγ₀ = t * γ₁, (1 - t) * γ₀
    γ = (ₜγ₁ + ₁₋ₜγ₀)
    (@fastmath ₜγ₁ / γ) * p₁(H) + (@fastmath ₁₋ₜγ₀ / γ) * p₀(H))
end

function evaluate!(u, H::ParameterHomotopy, x, t, c::NullCache)
    SP.evaluate!(u, H.F, x, p(H, t))
end

function jacobian!(u, H::ParameterHomotopy, x, t, c::NullCache)
    SP.jacobian!(u, H.F, x, p(H, t))
end

function evaluate_and_jacobian!(u, H::ParameterHomotopy, x, t, c::NullCache)
    SP.evaluate_and_jacobian!(u, H.F, x, p(H, t))
end

function dt!(u, H::ParameterHomotopy, x, t, c::NullCache)
    # apply chain rule to H(x, p(t))
    if h === nothing
        pₜ = p(t)
        ∂pₜ∂t = p₁(H) - p₀(H)
        ∂H∂p = SP.differentiate_parameters(H.F, x, pₜ)
        u .= ∂H∂p * ∂pₜ∂t
    else
        # copy from p(t) since we need some of the substructures
        ⁠γ₁, ⁠γ₀ = H.γ
        ₜγ₁, ₁₋ₜγ₀ = t * γ₁, (1 - t) * γ₀
        γ = (ₜγ₁ + ₁₋ₜγ₀)
        pₜ = (@fastmath ₜγ₁ / γ) * p₁(H) + (@fastmath ₁₋ₜγ₀ / γ) * p₀(H))
        # quotient rule to get the derivative of p(t)
        ∂pₜ∂t = (@fastmath H.γ₁ * H.γ₀ / (y * γ)) * (p₁(H) - p₀(H))
        ∂H∂p = SP.differentiate_parameters(H.F, x, pₜ)
        u .= ∂H∂p * ∂pₜ∂t
    end
end
