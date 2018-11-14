module Input

import MultivariatePolynomials
const MP = MultivariatePolynomials
using LinearAlgebra

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
struct StartTarget{P1<:Inputs, P2<:Inputs} <: AbstractInput
    start::P1
    target::P2
    startsolutions

    function StartTarget{P1, P2}(start::P1, target::P2, startsolutions) where {P1<:Inputs, P2<:Inputs}
        if length(start) ≠ length(target)
            error("Cannot construct `StartTargetProblem` since the lengths of `start` and `target` don't match.")
        end
        new(start, target, startsolutions)
    end
end
function StartTarget(start::P1, target::P2, startsolutions) where {P1<:Inputs, P2<:Inputs}
    StartTarget{P1, P2}(start, target, startsolutions)
end


struct Homotopy{Hom<:Homotopies.AbstractHomotopy} <: AbstractInput
    H::Hom
    startsolutions
end


"""
    TotalDegree(system::Vector{<:MP.AbstractPolynomial})

Construct a `TotalDegreeProblem`. This indicates that the system `system`
is the target system and a total degree system should be assembled.
"""
struct TotalDegree{S<:Inputs} <: AbstractInput
    system::S
    degrees::Vector{Int}
end
function TotalDegree(S::Vector{<:MP.AbstractPolynomialLike})
    TotalDegree(S, MP.maxdegree.(S))
end


"""
    ParameterSystem(system::Vector{<:MP.AbstractPolynomial}, x::Vector{<:MP.AbstractVariable}, p::Vector{<:MP.AbstractVariable})

Construct a `ParameterSystem`. This indicates that the system `system` has variables `x` and parameters `p`.
"""
struct ParameterSystem{P<:MP.AbstractPolynomialLike, V<:MP.AbstractVariable} <: AbstractInput
    system::Vector{P}
    parameters::Vector{V}
    p₁::AbstractVector
    p₀::AbstractVector
    startsolutions::AbstractVector
    γ₁::Union{Nothing, ComplexF64}
    γ₀::Union{Nothing, ComplexF64}
end

const overdetermined_error_msg = """
The input system is overdetermined. Therefore it is necessary to provide an explicit start system.
See
    https://www.JuliaHomotopyContinuation.org/guides/latest/overdetermined_tracking/
for details.
"""


const supported_keywords = [:parameters, :startparameters, :targetparameters,
                            :targetgamma, :startgamma, :p₁, :p₀, :γ₁, :γ₀]

"""
    input(F::Vector{<:MP.AbstractPolynomial})::TotalDegree
    input(F::Systems.AbstractSystem)::TotalDegree
    input(G::Vector{<:MP.AbstractPolynomial}, F::Vector{<:MP.AbstractPolynomial}, startsolutions)::StartTargetProblem
    input(F::Vector{<:MP.AbstractPolynomial}, parameters, startsolutions; kwargs...)::ParameterSystem
    input(H::Homotopies.AbstractHomotopy, startsolutions)::Homotopy

Construct an `AbstractInput`.
"""
function input(F::Vector{<:MP.AbstractPolynomial})
    F = filter(f -> !iszero(f), F)
    check_zero_dimensional(F)
    # square system and each polynomial is non-zero
    if length(F) == MP.nvariables(F) && ishomogenous(F)
        error("The input system is a square homogenous system. This will result in an at least 1 dimensional solution space.")
    end
    TotalDegree(F, MP.maxdegree.(F))
end

function input(F::Systems.AbstractSystem)
    n, N = size(F)
    degrees = check_homogenous_degrees(F)
    # system needs to be homogenous
    if n + 1 > N
        error(overdetermined_error_msg)
    elseif  n + 1 ≠ N
        error("Input system is not a square homogenous system!")
    end
    TotalDegree(F, degrees)
end

function input(G::Vector{<:MP.AbstractPolynomial}, F::Vector{<:MP.AbstractPolynomial}, startsolutions)
    if length(G) ≠ length(F)
        error("Start and target system don't have the same length")
    end
    check_zero_dimensional(F)
    StartTarget(G, F, startsolutions)
end

function input(F::Vector{<:MP.AbstractPolynomial}, startsolutions;
    parameters::Vector{<:MP.AbstractVariable}=error("parameters not defined"),
    startparameters=nothing, p₁ = startparameters,
    targetparameters=nothing, p₀ = targetparameters,
    startgamma=nothing, γ₁ = startgamma,
    targetgamma=nothing, γ₀ = targetgamma)

    if p₁ === nothing
        error("!`startparameters=` or `p₁=` need to be passed as argument")
    elseif p₀ === nothing
        error("!`targetparameters=` or `p₀=` need to be passed as argument")
    end

    if !(length(parameters) == length(p₁) == length(p₀))
        error("Number of parameters doesn't match!")
    end

    ParameterSystem(F, parameters, p₁, p₀, startsolutions, γ₁, γ₀)
end

function input(H::Homotopies.AbstractHomotopy, startsolutions)
    check_homogenous_degrees(Systems.FixedHomotopy(H, rand()))
    Homotopy(H, startsolutions)
end


"""
    check_homogenous_degrees(F::AbstractSystem)

Compute (numerically) the degrees of `F` and verify that `F` is homogenous,
"""
function check_homogenous_degrees(F::Systems.AbstractSystem)
    n, N = size(F)
    if n < N - 1
        error("Input system is not homogenous! It has $n polynomials in $N variables according to `size`.")
    end
    # The number of variables match, but it still cannot be homogenous.
    # We evaluate the system with y:=rand(N) and 2y. If homogenous then the output
    # scales accordingly to the degrees which we can obtain by taking logarithms.
    x = rand(ComplexF64, N)
    cache = Systems.cache(F, x)
    y = Systems.evaluate(F, x, cache)
    rmul!(x, 2)
    y2 = Systems.evaluate(F, x, cache)

    degrees = map(y2, y) do y2ᵢ, yᵢ
        # y2ᵢ = 2^dᵢ yᵢ
        float_dᵢ = log2(abs(y2ᵢ / yᵢ))
        dᵢ = round(Int, float_dᵢ)
        if abs(dᵢ - float_dᵢ) > 1e-10
            error("Input system is not homogenous by our numerical check.")
        end
        dᵢ
    end
    degrees
end

"""
    check_zero_dimensional(F::Vector{<:MP.AbstractPolynomial})

Check that the given polynomial system can have zero dimensional components.
"""
function check_zero_dimensional(F::Vector{<:MP.AbstractPolynomial})
    N = MP.nvariables(F)
    n = length(F)

    if n ≥ N || (n == N - 1 && ishomogenous(F))
        return nothing
    end
    error("The input system will not result in a finite number of solutions.")
end

end
