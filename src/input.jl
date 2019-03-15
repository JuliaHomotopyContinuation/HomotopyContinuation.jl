export AbstractInput,
    StartTargetInput,
    TotalDegreeInput,
    HomotopyInput,
    ParameterSystemInput

const Inputs = Union{<:AbstractSystem, <:MPPolys, <:Composition}
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
struct HomotopyInput{Hom<:AbstractHomotopy} <: AbstractInput
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
    input(F::AbstractSystem)::TotalDegreeInput
    input(G::MPPolyInputs, F::MPPolyInputs, startsolutions)::StartTargetInput
    input(F::MPPolyInputs, parameters, startsolutions; kwargs...)::ParameterSystemInput
    input(H::AbstractHomotopy, startsolutions)::HomotopyInput

Construct an `AbstractInput`.
"""
function input(F::MPPolyInputs; parameters=nothing, kwargs...)
    # if parameters === nothing this is actually the
    # input constructor for a parameter homotopy, but no startsolutions
    # are provided
    if parameters !== nothing
        startsolutions = [randn(ComplexF64, nvariables(F, parameters=parameters))]
        return input(F, startsolutions; parameters=parameters, kwargs...)
    end

    remove_zeros!(F)
    # check_zero_dimensional(F)
    # square system and each polynomial is non-zero
    if length(F) == nvariables(F) && ishomogenous(F)
        error("Cannot construct a total degree homotopy for a square homogenous system.")
    end
    TotalDegreeInput(F)
end

function input(F::AbstractSystem)
    TotalDegreeInput(F)
end

function input(G::MPPolyInputs, F::MPPolyInputs, startsolutions=nothing)
    if length(G) ≠ length(F)
        error("Start and target system don't have the same length")
    end
    check_zero_dimensional(F)
    if startsolutions === nothing
        startsolutions = [randn(ComplexF64, nvariables(F))]
    elseif isa(startsolutions, AbstractVector{<:Number})
        startsolutions = [startsolutions]
    end
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
    if startsolutions === nothing
        startsolutions = [randn(ComplexF64, nvariables(F, parameters=parameters))]
    elseif isa(startsolutions, AbstractVector{<:Number})
        startsolutions = [startsolutions]
    end
    ParameterSystemInput(F, parameters, p₁, p₀, startsolutions, γ₁, γ₀)
end

function input(H::AbstractHomotopy, startsolutions)
    check_homogenous_degrees(FixedHomotopy(H, rand()))
    if isa(startsolutions, AbstractVector{<:Number})
        startsolutions = [startsolutions]
    end
    HomotopyInput(H, startsolutions)
end
