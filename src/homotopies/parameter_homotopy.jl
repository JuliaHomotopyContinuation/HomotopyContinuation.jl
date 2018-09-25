export StaticParameterHomotopy, ParameterHomotopyCache, gamma

import MultivariatePolynomials
const MP = MultivariatePolynomials
import StaticPolynomials
const SP = StaticPolynomials
import StaticArrays: MVector, SVector

import ..Utilities

"""
    ParameterHomotopy(F, variables, parameters, p₁, p₀)

Construct the homotopy ``H(x, t) = F(x, t * p₁ + (1-t) * p₀)``,
where `start` and `target` are a vector of parameters of ``F``.
The input ``parameters`` specifies the parameter variables of ``F``.
Neccessarily, ``length(parameters) == length(p₁) == length(p₀)``.
"""
struct ParameterHomotopy{N, NVars, NParams, T<:Number, PolyTuple} <: AbstractHomotopy
    F::SP.PolynomialSystem{N, NVars, NParams, PolyTuple}
    p₁::MVector{NParams, T}
    p₀::MVector{NParams, T}
    γ₁::ComplexF64
    γ₀::ComplexF64
end

function ParameterHomotopy(F::SP.PolynomialSystem{N, NVars, NParams, T}, p₁, p₀;
    start_gamma=randn(ComplexF64), γ₀ = start_gamma,
    target_gamma=randn(ComplexF64), γ₁ = target_gamma
    ) where {N, NVars, NParams, T}
    @assert length(p₁) == length(p₀) == NParams

    ParameterHomotopy(F, promote(MVector{NParams}(p₁), MVector{NParams}(p₀))..., γ₁, γ₀)
end

function ParameterHomotopy(F::Vector{T}, parameters::AbstractVector{V};
    start_parameters=randn(ComplexF64, length(parameters)), p₁ = start_parameters,
    target_parameters=randn(ComplexF64, length(parameters)), p₀ = target_parameters,
    kwargs...) where {T<:MP.AbstractPolynomialLike, V<:MP.AbstractVariable}
    G = SP.PolynomialSystem(F; parameters=parameters)
    ParameterHomotopy(G, p₁, p₀; kwargs...)
end

function ParameterHomotopy(F::Vector{T}, variables, parameters, p₁, p₀) where {T<:MP.AbstractPolynomialLike}
    G = SP.PolynomialSystem(F; variables=variables, parameters=parameters)
    ParameterHomotopy(G, p₁, p₀)
end

(H::ParameterHomotopy)(x, t, c=NullCache()) = evaluate(H, x, t, c)

cache(::ParameterHomotopy, x, t) = NullCache()

Base.size(H::ParameterHomotopy{N, NVars}) where {N, NVars} = (N, NVars)

p₁(H::ParameterHomotopy) = SVector(H.p₁)
p₀(H::ParameterHomotopy) = SVector(H.p₀)

function evaluate!(u, H::ParameterHomotopy, x, t, c::NullCache)
    SP.evaluate!(u, H.F, x, t * p₁(H) + (1 - t) * p₀(H))
end

function jacobian!(u, H::ParameterHomotopy, x, t, c::NullCache)
    SP.jacobian!(u, H.F, x, t * p₁(H) + (1 - t) * p₀(H))
end

function evaluate_and_jacobian!(u, H::ParameterHomotopy, x, t, c::NullCache)
    SP.evaluate_and_jacobian!(u, H.F, x, t * p₁(H) + (1 - t) * p₀(H))
end

function dt!(u, H::ParameterHomotopy, x, t, c::NullCache)
    # apply chain rule to H(x, p(t))
    p = t * p₁(H) + (1 - t) * p₀(H)
    ∂p∂t = p₁(H) - p₀(H)
    ∂H∂p = SP.differentiate_parameters(H.F, x, p)
    u .= ∂H∂p * ∂p∂t
end
