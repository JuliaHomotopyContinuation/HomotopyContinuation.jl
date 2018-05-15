module Problems

import MultivariatePolynomials
const MP = MultivariatePolynomials


import ..Homotopies: AbstractHomotopy, StraightLineHomotopy
import ..ProjectiveVectors
import ..Systems: SPSystem

using ..Utilities

# STAGE 1 exports
export AbstractDynamicProblem,
    AbstractStartTargetProblem,
    AbstractParameterProblem,
    StartTargetProblem,
    TotalDegreeProblem

# STAGE 2 exports
export AbstractProblem,
    AbstractProjectiveProblem,
    AbstractHomogenizationStrategy,
    NullHomogenization,
    DefaultHomogenization,
    ProjectiveStartTargetProblem,
    homotopy,
    homogenization_strategy,
    embed

# STAGE 1

abstract type AbstractDynamicProblem end

abstract type AbstractStartTargetProblem <: AbstractDynamicProblem end
abstract type AbstractParameterProblem <: AbstractDynamicProblem end

"""
    StartTargetProblem(start::Vector{<:MP.AbstractPolynomial}, target::Vector{<:MP.AbstractPolynomial})

Construct a `StartTargetProblem` out of two systems of polynomials.
"""
struct StartTargetProblem{P1<:MP.AbstractPolynomialLike, P2<:MP.AbstractPolynomialLike} <: AbstractStartTargetProblem
    start::Vector{P1}
    target::Vector{P2}

    function StartTargetProblem{P1, P2}(start::Vector{P1}, target::Vector{P2}) where {P1<:MP.AbstractPolynomialLike, P2<:MP.AbstractPolynomialLike}
        @assert length(start) == length(target) "Cannot construct `StartTargetProblem` since the lengths of `start` and `target` don't match."
        @assert nvariables(start) == nvariables(target) "Cannot construct `StartTargetProblem` since the number of variables of `start` and `target` don't match."
        @assert ishomogenous(start) == ishomogenous(target) "Cannot construct `StartTargetProblem` since `start` is homogenous and  `target` not or the other way around."
        new(start, target)
    end
end
function StartTargetProblem(start::Vector{P1}, target::Vector{P2}) where {P1<:MP.AbstractPolynomialLike, P2<:MP.AbstractPolynomialLike}
    StartTargetProblem{P1, P2}(start, target)
end


"""
    TotalDegreeProblem(system::Vector{<:MP.AbstractPolynomial})

Construct a `TotalDegreeProblem`. This indicates that the system `system`
is the target system and a total degree system should be assembled.
"""
struct TotalDegreeProblem{P1<:MP.AbstractPolynomialLike} <: AbstractStartTargetProblem
    system::Vector{P1}
end


# STAGE 2

abstract type AbstractProblem end
abstract type AbstractProjectiveProblem <: AbstractProblem end


# The HomogenizationStrategy will be useful if we have multi-homogenous stuff
"""
    AbstractHomogenizationStrategy

Abstract type for describing the homogenization strategy.
"""
abstract type AbstractHomogenizationStrategy end

"""
    NullHomogenization()

Homogenization was not necessary since it is already homogenous.
"""
struct NullHomogenization <: AbstractHomogenizationStrategy end

"""
    DefaultHomogenization()

Homogenization by adding a new variable with highest precedence.
"""
struct DefaultHomogenization <: AbstractHomogenizationStrategy end

"""
    ProjectiveStartTargetProblem(prob::TotalDegreeProblem; options...)
    ProjectiveStartTargetProblem(prob::StartTargetProblem; options...)
    ProjectiveStartTargetProblem(start<:Vector{<:Polynomial}, target::Vector{<:MP.Polynomial}; options...)

Construct a `ProjectiveStartTargetProblem`. This steps constructs a homotopy and homogenizes
the systems if necessary. The options are
* `system=SPSystem`: A constructor to assemble a [`Systems.AbstractSystem`](@ref). The constructor
is called with `system(polys, variables)` where `variables` determines the variable ordering.
* `homotopy=StraightLineHomotopy`: A constructor to construct a [`Homotopies.AbstractHomotopy`](@ref) an `Systems.AbstractSystem`. The constructor
is called with `homotopy(start, target)` where `start` and `target` are systems constructed
with `system`.

    ProjectiveStartTargetProblem(H::AbstractHomotopy)

Construct a `ProjectiveProblem`. The homotopy `H` needs to be homogenous.
"""
struct ProjectiveStartTargetProblem{H<:AbstractHomotopy, HS<:AbstractHomogenizationStrategy} <: AbstractProjectiveProblem
    homotopy::H
    homogenization_strategy::HS
end

function ProjectiveStartTargetProblem(prob::TotalDegreeProblem; kwargs...)
    ProjectiveStartTargetProblem(totaldegree(prob.system), prob.system; kwargs...)
end

function ProjectiveStartTargetProblem(prob::StartTargetProblem; kwargs...)
    ProjectiveStartTargetProblem(prob.start, prob.target; kwargs...)
end
function ProjectiveStartTargetProblem(
    s::Vector{<:MP.AbstractPolynomialLike},
    t::Vector{<:MP.AbstractPolynomialLike};
    system=SPSystem,
    homotopy=StraightLineHomotopy)

    if ishomogenous(s)
        homogenization_strategy = NullHomogenization()
        start = system(s)
        target = system(t)
    # TODO: Multihomogenous...
    else
        homogenization_strategy = DefaultHomogenization()
        vars = allvariables(s)
        homvar = uniquevar(s)
        var_ordering = [homvar; vars]

        start = system(homogenize(s, homvar), var_ordering)
        target = system(homogenize(t, homvar), var_ordering)
    end
    H = homotopy(start, target)

    ProjectiveStartTargetProblem(H, homogenization_strategy)
end
function ProjectiveStartTargetProblem(H::AbstractHomotopy)
    ProjectiveStartTargetProblem(H, NullHomogenization())
end

"""
    homotopy(prob::ProjectiveStartTargetProblem)

Get the homotopy stored in the problem `prob`.
"""
homotopy(prob::ProjectiveStartTargetProblem) = prob.homotopy

"""
    homogenization_strategy(prob::ProjectiveStartTargetProblem)

Get the [`HomogenizationStrategy`](@ref) stored in the problem `prob`.
"""
homogenization_strategy(prob::ProjectiveStartTargetProblem) = prob.homogenization_strategy

"""
    embed(prob::ProjectiveStartTargetProblem, x)

Embed the solution `x` into projective space if necessary.
"""
function embed(prob::ProjectiveStartTargetProblem{<:AbstractHomotopy, NullHomogenization}, x)
    M, N = size(homotopy(prob))
    length(x) != N && throw(error("The length of the intial solution is $(length(x)) but expected length $N."))
    return x
end
function embed(prob::ProjectiveStartTargetProblem{<:AbstractHomotopy, DefaultHomogenization}, x)
    M, N = size(homotopy(prob))
    if length(x) == N
        return ProjectiveVectors.PVector(x, 1)
    elseif length(x) == N - 1
        return ProjectiveVectors.embed(x, 1)
    end
    throw(error("The length of the intial solution is $(length(x)) but expected length $N or $(N-1)."))
end



end
