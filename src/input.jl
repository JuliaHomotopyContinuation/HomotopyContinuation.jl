module Input

import MultivariatePolynomials
const MP = MultivariatePolynomials

import ..Homotopies
import ..Systems
using ..Utilities

# STAGE 1 exports
export AbstractInput,
    StartTarget,
    TotalDegree,
    ParameterSystem


abstract type AbstractInput end

const MPPolys = Vector{<:MP.AbstractPolynomialLike}
const Inputs = Union{Systems.AbstractSystem, MPPolys}

"""
    StartTargetProblem(start::Vector{<:MP.AbstractPolynomial}, target::Vector{<:MP.AbstractPolynomial})

Construct a `StartTargetProblem` out of two systems of polynomials.
"""
struct StartTarget{P1<:Inputs, P2<:Inputs, V<:AbstractVector} <: AbstractInput
    start::P1
    target::P2
    startsolutions::Vector{V}

    function StartTarget{P1, P2, V}(start::P1, target::P2, startsolutions::Vector{V}) where {P1<:Inputs, P2<:Inputs, V<:AbstractVector}
        @assert length(start) == length(target) "Cannot construct `StartTargetProblem` since the lengths of `start` and `target` don't match."
        new(start, target, startsolutions)
    end
end
function StartTarget(start::P1, target::P2, startsolutions::Vector{V}) where {P1<:Inputs, P2<:Inputs, V<:AbstractVector}
    StartTarget{P1, P2, V}(start, target, startsolutions)
end


struct Homotopy{Hom<:Homotopies.AbstractHomotopy, V} <: AbstractInput
    H::Hom
    startsolutions::Vector{V}
end


"""
    TotalDegree(system::Vector{<:MP.AbstractPolynomial})

Construct a `TotalDegreeProblem`. This indicates that the system `system`
is the target system and a total degree system should be assembled.
"""
struct TotalDegree{S<:Inputs} <: AbstractInput
    system::S
end


"""
    ParameterSystem(system::Vector{<:MP.AbstractPolynomial}, x::Vector{<:MP.AbstractVariable}, p::Vector{<:MP.AbstractVariable})

Construct a `ParameterSystem`. This indicates that the system `system` has variables `x` and parameters `p`.
"""
struct ParameterSystem{P<:MP.AbstractPolynomialLike, V<:MP.AbstractVariable, T1<:Number, T2<:Number, V1<:AbstractVector} <: AbstractInput
    system::Vector{P}
    parameters::Vector{V}
    p₁::Vector{T1}
    p₀::Vector{T2}
    startsolutions::Vector{V1}
end


end
