import .Homotopies
import .Systems
using .Utilities

export AbstractInput,
    StartTargetInput,
    TotalDegreeInput,
    HomotopyInput,
    ParameterSystemInput

const MPPolys = Vector{<:MP.AbstractPolynomialLike}
const Inputs = Union{<:Systems.AbstractSystem, <:MPPolys, <:Composition}
const MPPolyInputs = Union{<:MPPolys, <:Composition}

const input_supported_keywords = [:parameters, :startparameters, :targetparameters,
                            :targetgamma, :startgamma, :p₁, :p₀, :γ₁, :γ₀]



"""
    AbstractInput

An abstract type to model different input types.
"""
abstract type AbstractInput end

"""
    StartTargetInput(start::MPPolyInputs, target::MPPolyInputs, startsolutions)

Construct a `StartTargetInput` out of two systems of polynomials.
"""
struct StartTargetInput{P1<:Inputs, P2<:Inputs} <: AbstractInput
    start::P1
    target::P2
    startsolutions

    function StartTargetInput{P1, P2}(start::P1, target::P2, startsolutions) where {P1<:Inputs, P2<:Inputs}
        if length(start) ≠ length(target)
            error("Cannot construct `StartTargetInput` since the lengths of `start` and `target` don't match.")
        end
        new(start, target, startsolutions)
    end
end
function StartTargetInput(start::P1, target::P2, startsolutions) where {P1<:Inputs, P2<:Inputs}
    StartTargetInput{P1, P2}(start, target, startsolutions)
end


"""
    HomotopyInput(H::Homotopy, startsolutions)

Construct a `HomotopyInput` for a homotopy `H` with given `startsolutions`.
"""
struct HomotopyInput{Hom<:Homotopies.AbstractHomotopy} <: AbstractInput
    H::Hom
    startsolutions
end


"""
    TotalDegreeInput(system::MPPolyInputs)

Construct a `TotalDegreeInputProblem`. This indicates that the system `system`
is the target system and a total degree system should be assembled.
"""
struct TotalDegreeInput{S<:Inputs} <: AbstractInput
    system::S
    degrees::Vector{Int}
end
function TotalDegreeInput(S::Vector{<:MP.AbstractPolynomialLike})
    TotalDegreeInput(S, maxdegrees(S))
end


"""
    ParameterSystemInput(F, parameters, p₁, p₀, startsoluions, γ₁=nothing, γ₂=nothing)

Construct a `ParameterSystemInput`.
"""
struct ParameterSystemInput{P<:MPPolyInputs, V<:MP.AbstractVariable} <: AbstractInput
    system::P
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

"""
    input(F::MPPolyInputs)::TotalDegreeInput
    input(F::Systems.AbstractSystem)::TotalDegreeInput
    input(G::MPPolyInputs, F::MPPolyInputs, startsolutions)::StartTargetInput
    input(F::MPPolyInputs, parameters, startsolutions; kwargs...)::ParameterSystemInput
    input(H::Homotopies.AbstractHomotopy, startsolutions)::HomotopyInput

Construct an `AbstractInput`.
"""
function input(F::MPPolyInputs)
    remove_zeros!(F)
    check_zero_dimensional(F)
    # square system and each polynomial is non-zero
    if length(F) == nvariables(F) && ishomogenous(F)
        error("Cannot construct a total degree homotopy for a square homogenous system.")
    end
    TotalDegreeInput(F, maxdegrees(F))
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
    TotalDegreeInput(F, degrees)
end

function input(G::MPPolyInputs, F::MPPolyInputs, startsolutions)
    if length(G) ≠ length(F)
        error("Start and target system don't have the same length")
    end
    check_zero_dimensional(F)
    StartTargetInput(G, F, startsolutions)
end

function input(F::MPPolyInputs, startsolutions;
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

    ParameterSystemInput(F, parameters, p₁, p₀, startsolutions, γ₁, γ₀)
end

function input(H::Homotopies.AbstractHomotopy, startsolutions)
    check_homogenous_degrees(Systems.FixedHomotopy(H, rand()))
    HomotopyInput(H, startsolutions)
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
