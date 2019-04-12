export AbstractInput,
    StartTargetInput,
    TotalDegreeInput,
    HomotopyInput,
    ParameterSystemInput

const Inputs = Union{<:AbstractSystem, <:MPPolys, <:Composition}
const MPPolyInputs = Union{<:MPPolys, <:Composition}

const input_supported_keywords = [:parameters, :generic_parameters, :startparameters, :targetparameters,
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

    function StartTargetInput{P1, P2}(start::P1, target::P2) where {P1<:Inputs, P2<:Inputs}
        if length(start) ≠ length(target)
            error("Cannot construct `StartTargetInput` since the lengths of `start` and `target` don't match.")
        end
        new(start, target)
    end
end
function StartTargetInput(start::P1, target::P2) where {P1<:Inputs, P2<:Inputs}
    StartTargetInput{P1, P2}(start, target)
end


"""
    HomotopyInput(H::Homotopy, startsolutions)

Construct a `HomotopyInput` for a homotopy `H` with given `startsolutions`.
"""
struct HomotopyInput{Hom<:AbstractHomotopy} <: AbstractInput
    H::Hom
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
struct ParameterSystemInput{S<:Inputs} <: AbstractInput
    system::S
    parameters::Union{Nothing, Vector{<:MP.AbstractVariable}}
    p₁::AbstractVector
    p₀::AbstractVector
    γ₁::Union{Nothing, ComplexF64}
    γ₀::Union{Nothing, ComplexF64}
end

"""
    input_startsolutions(F::MPPolyInputs)
    input_startsolutions(F::AbstractSystem)
    input_startsolutions(G::MPPolyInputs, F::MPPolyInputs, startsolutions)
    input_startsolutions(F::MPPolyInputs, parameters, startsolutions; kwargs...)
    input_startsolutions(H::AbstractHomotopy, startsolutions)

Construct an `AbstractInput` and pass through startsolutions if provided.
Returns a named tuple `(input=..., startsolutions=...)`.
"""
function input_startsolutions(F::MPPolyInputs; parameters=nothing, kwargs...)
    # if parameters !== nothing this is actually the
    # input constructor for a parameter homotopy, but no startsolutions
    # are provided
    if parameters !== nothing
        startsolutions = [randn(ComplexF64, nvariables(F, parameters=parameters))]
        return input_startsolutions(F, startsolutions; parameters=parameters, kwargs...)
    end

    remove_zeros!(F)
    # check_zero_dimensional(F)
    # square system and each polynomial is non-zero
    if length(F) == nvariables(F) && ishomogeneous(F)
        error("Cannot construct a total degree homotopy for a square homogeneous system.")
    end
    (input=TotalDegreeInput(F), startsoluions=nothing)
end

function input_startsolutions(F::AbstractSystem)
    (input=TotalDegreeInput(F), startsolutions=nothing)
end

function input_startsolutions(G::MPPolyInputs, F::MPPolyInputs, startsolutions=nothing)
    if length(G) ≠ length(F)
        error("Start and target system don't have the same length")
    end
    check_zero_dimensional(F)
    if startsolutions === nothing
        startsolutions = [randn(ComplexF64, nvariables(F))]
    elseif isa(startsolutions, AbstractVector{<:Number})
        startsolutions = [startsolutions]
    end
    (input=StartTargetInput(G, F), startsolutions=startsolutions)
end


# need
input_startsolutions(F::MPPolyInputs, starts; kwargs...) = parameter_homotopy(F, starts; kwargs...)
input_startsolutions(F::AbstractSystem, starts; kwargs...) = parameter_homotopy(F, starts; kwargs...)

function parameter_homotopy(F::Inputs, startsolutions;
    parameters=(isa(F, AbstractSystem) ? nothing : error(ArgumentError("You need to pass `parameters=...` as a keyword argument."))),
    generic_parameters=nothing,
    startparameters=generic_parameters, p₁ = startparameters,
    targetparameters=generic_parameters, p₀ = targetparameters,
    startgamma=nothing, γ₁ = startgamma,
    targetgamma=nothing, γ₀ = targetgamma)

    if p₁ === nothing
        error("You need to pass `generic_parameters=`, `startparameters=` or `p₁=` as a keyword argument")
    elseif p₀ === nothing
        error("`targetparameters=` or `p₀=` need to be passed as a keyword argument.")
    end

    if length(p₁) != length(p₀) || (parameters !== nothing && length(parameters) != length(p₀))
        error("Number of parameters doesn't match!")
    end
    if startsolutions === nothing && parameters !== nothing
        startsolutions = [randn(ComplexF64, nvariables(F, parameters=parameters))]
    elseif isa(startsolutions, AbstractVector{<:Number})
        startsolutions = [startsolutions]
    end
    (input=ParameterSystemInput(F, parameters, p₁, p₀, γ₁, γ₀), startsolutions=startsolutions)
end

function input_startsolutions(H::AbstractHomotopy, startsolutions)
    if isa(startsolutions, AbstractVector{<:Number})
        startsolutions = [startsolutions]
    end
    (input=HomotopyInput(H), startsolutions=startsolutions)
end
