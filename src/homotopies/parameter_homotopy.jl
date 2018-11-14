export ParameterHomotopy, nparameters

import MultivariatePolynomials
const MP = MultivariatePolynomials
import StaticPolynomials
const SP = StaticPolynomials
import StaticArrays: SVector

import ..Utilities

"""
    ParameterHomotopy(F, parameters;
        variables=setdiff(MP.variables(F), parameters),
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
        variables=setdiff(MP.variables(F), parameters),
        startparameters=randn(ComplexF64, length(parameters)),
        targetparameters=randn(ComplexF64, length(parameters)),
        startgamma=nothing, targetgamma=nothing)

This is a non-unicode variant where `γ₁=startparameters`, `γ₀=targetparameters`,
`γ₁=startgamma`, `γ₀=targetgamma`.
"""
mutable struct ParameterHomotopy{N, NVars, NParams, T<:Number, PolyTuple} <: AbstractHomotopy
    F::SP.PolynomialSystem{N, NVars, NParams, PolyTuple}
    p::NTuple{2, SVector{NParams, T}}
    γ::Union{Nothing, NTuple{2, ComplexF64}}
end

function ParameterHomotopy(F::SP.PolynomialSystem{N, NVars, NParams, T}, p₁::AbstractVector, p₀::AbstractVector;
    startgamma=nothing, γ₀ = startgamma,
    targetgamma=nothing, γ₁ = targetgamma
    ) where {N, NVars, NParams, T}

    if !(length(p₁) == length(p₀) == NParams)
        error("Length of parameters provided doesn't match the number of parameters.")
    end

    if γ₁ === nothing || γ₀ === nothing
        γ = nothing
    else
        γ = (γ₁, γ₀)
    end
    p = promote(SVector{NParams}(p₁), SVector{NParams}(p₀))
    ParameterHomotopy(F, p, γ)
end

function ParameterHomotopy(F::Vector{T},
    parameters::AbstractVector{V};
    variables=setdiff(MP.variables(F), parameters),
    startparameters=randn(ComplexF64, length(parameters)), p₁ = startparameters,
    targetparameters=randn(ComplexF64, length(parameters)), p₀ = targetparameters,
    kwargs...) where {T<:MP.AbstractPolynomialLike, V<:MP.AbstractVariable}
    G = SP.PolynomialSystem(F; variables=variables, parameters=parameters)
    ParameterHomotopy(G, p₁, p₀; kwargs...)
end

(H::ParameterHomotopy)(x, t, c=NullCache()) = evaluate(H, x, t, c)

cache(::ParameterHomotopy, x, t) = NullCache()

Base.size(H::ParameterHomotopy{N, NVars}) where {N, NVars} = (N, NVars)

"""
    nparameters(H::ParameterHomotopy)

Returns the number of parameters of `H`.
"""
function nparameters(::ParameterHomotopy{N, NVars, NParams}) where {N, NVars, NParams}
    NParams
end


"""
    set_parameters!(H::ParameterHomotopy, p::NTuple{2, <:SVector}, γ)

Update the parameters `p` and `γ` of `H`.
"""
function set_parameters!(H::Homotopies.ParameterHomotopy{N, NVars, NParams},
    p::NTuple{2, SVector{NParams, T}},
    γ::Union{Nothing, NTuple{2, ComplexF64}}) where {N, NVars, NParams, T}
    H.p = p
    H.γ = γ
    H
end

"""
    set_parameters!(H::ParameterHomotopy, p::NTuple{2, SVector})

Update the parameter `p` of `H`.
"""
function set_parameters!(H::Homotopies.ParameterHomotopy{N, NVars, NParams},
    p::NTuple{2, SVector{NParams, T}}) where {N, NVars, NParams, T}
    H.p = p
    H
end

@inline function p(H::ParameterHomotopy, t)
    p₁, p₀ = H.p
    if H.γ === nothing
        return t * p₁ + (1 - t) * p₀
    else
        # compute (tγ₁p₁+(1-t)γ₀p₀) / (tγ₁+(1-t)γ₀)
        γ₁, γ₀ = H.γ
        ₜγ₁, γ₀_₁₋ₜ = t * γ₁, (1 - t) * γ₀
        γ = (ₜγ₁ + γ₀_₁₋ₜ)
        return (@fastmath ₜγ₁ / γ) * p₁ + (@fastmath γ₀_₁₋ₜ / γ) * p₀
    end
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
    p₁, p₀ = H.p
    if H.γ === nothing
        pₜ = p(H, t)
        ∂pₜ∂t = p₁ - p₀
        ∂H∂p = SP.differentiate_parameters(H.F, x, pₜ)
        u .= ∂H∂p * ∂pₜ∂t
    else
        # copy from p(t) since we need some of the substructures
        γ₁, γ₀ = H.γ
        ₜγ₁, γ₀_₁₋ₜ = t * γ₁, (1 - t) * γ₀
        γ = (ₜγ₁ + γ₀_₁₋ₜ)
        pₜ = (@fastmath ₜγ₁ / γ) * p₁ + (@fastmath γ₀_₁₋ₜ / γ) * p₀
        # quotient rule to get the derivative of p(t)
        ∂pₜ∂t = (@fastmath γ₁ * γ₀ / (γ * γ)) * (p₁ - p₀)
        ∂H∂p = SP.differentiate_parameters(H.F, x, pₜ)
        u .= ∂H∂p * ∂pₜ∂t
    end
end
