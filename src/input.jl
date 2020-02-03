export AbstractInput,
    StartTargetInput, TargetSystemInput, HomotopyInput, ParameterSystemInput

const Inputs = Union{AbstractSystem,MPPolys,Composition,ModelKit.System}
const MPPolyInputs = Union{MPPolys,Composition}

const input_supported_keywords = [
    :parameters,
    :generic_parameters,
    :start_parameters,
    :target_parameters,
    :target_gamma,
    :start_gamma,
    :p₁,
    :p₀,
    :γ₁,
    :γ₀,
    :variable_ordering,
]



"""
    AbstractInput

An abstract type to model different input types.
"""
abstract type AbstractInput end

"""
    StartTargetInput(start::MPPolyInputs, target::MPPolyInputs, startsolutions)

Construct a `StartTargetInput` out of two systems of polynomials.
"""
struct StartTargetInput{P1<:Inputs,P2<:Inputs} <: AbstractInput
    start::P1
    target::P2

    function StartTargetInput{P1,P2}(start::P1, target::P2) where {P1<:Inputs,P2<:Inputs}
        if length(start) ≠ length(target)
            error("Cannot construct `StartTargetInput` since the lengths of `start` and `target` don't match.")
        end
        new(start, target)
    end
end
function StartTargetInput(start::P1, target::P2) where {P1<:Inputs,P2<:Inputs}
    StartTargetInput{P1,P2}(start, target)
end


"""
    HomotopyInput(H::Homotopy, startsolutions)

Construct a `HomotopyInput` for a homotopy `H` with given `startsolutions`.
"""
struct HomotopyInput{Hom<:AbstractHomotopy} <: AbstractInput
    H::Hom
end


"""
    TargetSystemInput(system::Inputs)

Construct a `TargetSystemInput`. This indicates that the system `system`
is the target system and a start system should be assembled.
"""
struct TargetSystemInput{S<:Inputs} <: AbstractInput
    system::S
end

"""
    ParameterSystemInput(F, parameters, p₁, p₀, startsolutions, γ₁=nothing, γ₂=nothing)

Construct a `ParameterSystemInput`.
"""
struct ParameterSystemInput{S<:Inputs} <: AbstractInput
    system::S
    parameters::Union{Nothing,Vector{<:MP.AbstractVariable}}
    p₁::AbstractVector
    p₀::AbstractVector
    γ₁
    γ₀
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
function input_startsolutions(
    F::Union{MPPolyInputs,ModelKit.System};
    variable_ordering = nothing,
    # parameter homotopy indication,
    generic_parameters = nothing,
    start_parameters = generic_parameters,
    p₁ = start_parameters,
    target_parameters = generic_parameters,
    p₀ = target_parameters,
    kwargs...,
)
    # if p₀ or p₁ is provided this is actually the input constructor
    # for a parameter homotopy, but no startsolutions are provided
    if !isnothing(p₁) || !isnothing(p₀)
        return parameter_homotopy(
            F,
            nothing;
            variable_ordering = variable_ordering,
            p₁ = p₁,
            p₀ = p₀,
            kwargs...,
        )
    end

    if isa(F, MPPolyInputs)
        if variable_ordering !== nothing && nvariables(F) != length(variable_ordering)
            throw(ArgumentError("Number of assigned variables is too small."))
        end
        remove_zeros!(F)
        if has_constant_polynomial(F)
            throw(ArgumentError("System contains a non-zero constant polynomial."))
        end
    end

    (input = TargetSystemInput(F), startsolutions = nothing)
end

function input_startsolutions(F::AbstractSystem; variable_ordering = nothing)
    (input = TargetSystemInput(F), startsolutions = nothing)
end

# create ModelKit.System and use this case
function input_startsolutions(
    F::Vector{<:ModelKit.Expression},
    args...;
    variable_ordering::Union{Nothing,Vector{ModelKit.Variable}} = nothing,
    parameters::Union{Nothing,Vector{ModelKit.Variable}} = nothing,
    kwargs...,
)
    if variable_ordering === nothing
        throw(ArgumentError("`variable_ordering = ...` needs to be passed as a keyword argument to indicate the order of the variables."))
    end
    if isnothing(parameters)
        system = ModelKit.System(F, variable_ordering)
    else
        system = ModelKit.System(F, variable_ordering, parameters)
    end
    input_startsolutions(system, args...; kwargs...)
end

function input_startsolutions(
    G::MPPolyInputs,
    F::MPPolyInputs,
    startsolutions = nothing;
    variable_ordering = nothing,
)
    if length(G) ≠ length(F)
        throw(ArgumentError("Start and target system don't have the same length"))
    end
    if variable_ordering !== nothing && (
        nvariables(F) != length(variable_ordering) ||
        nvariables(G) != length(variable_ordering)
    )
        throw(ArgumentError("Number of assigned variables is too small."))
    end

    check_zero_dimensional(F)
    if startsolutions === nothing
        startsolutions = [randn(ComplexF64, nvariables(F))]
    elseif isa(startsolutions, AbstractVector{<:Number})
        startsolutions = [startsolutions]
    end
    (input = StartTargetInput(G, F), startsolutions = startsolutions)
end


# need
input_startsolutions(F::Inputs, starts; kwargs...) =
    parameter_homotopy(F, starts; kwargs...)

function parameter_homotopy(
    F::Inputs,
    startsolutions;
    parameters = F isa MPPolyInputs ?
                 throw(ArgumentError("Argument `parameters = ...` needs to be provided.")) :
                 nothing,
    variable_ordering = nothing,
    generic_parameters = nothing,
    start_parameters = generic_parameters,
    p₁ = start_parameters,
    target_parameters = generic_parameters,
    p₀ = target_parameters,
    start_gamma = nothing,
    γ₁ = start_gamma,
    target_gamma = nothing,
    γ₀ = target_gamma,
)
    if isnothing(p₁)
        throw(ArgumentError(
            "You need to pass `generic_parameters = ...`, " *
            "`start_parameters = ...` or `p₁ = ...` as a keyword argument",
        ))
    end
    if isnothing(p₀)
        throw(ArgumentError(
            "`target_parameters = ...` or `p₀ = ...` need " *
            "to be passed as a keyword argument.",
        ))
    end
    if length(p₁) != length(p₀)
        throw(ArgumentError("Number of parameters doesn't match!"))
    end
    if F isa MPPolyInputs && length(parameters) != length(p₀)
        throw(ArgumentError("Number of parameters doesn't match!"))
    end
    if F isa MPPolyInputs && !isnothing(variable_ordering) &&
       nvariables(F, parameters) != length(variable_ordering)
        throw(ArgumentError("Number of assigned variables is too small."))
    end

    if F isa MPPolyInputs && variable_ordering !== nothing &&
       nvariables(F, parameters) != length(variable_ordering)
        throw(ArgumentError("Number of assigned variables is too small."))
    end

    if isnothing(startsolutions)
        startsolutions = [randn(ComplexF64, nvariables(F, parameters))]
    elseif isa(startsolutions, AbstractVector{<:Number})
        startsolutions = [startsolutions]
    end

    (
        input = ParameterSystemInput(F, parameters, p₁, p₀, γ₁, γ₀),
        startsolutions = startsolutions,
    )
end

function input_startsolutions(
    H::AbstractHomotopy,
    startsolutions;
    variable_ordering = nothing,
)
    if isa(startsolutions, AbstractVector{<:Number})
        startsolutions = [startsolutions]
    end
    (input = HomotopyInput(H), startsolutions = startsolutions)
end

function input_startsolutions(
    H::ModelKit.Homotopy,
    startsolutions;
    variable_ordering = nothing,
)
    if isa(startsolutions, AbstractVector{<:Number})
        startsolutions = [startsolutions]
    end
    (input = HomotopyInput(ModelKitHomotopy(H)), startsolutions = startsolutions)
end
