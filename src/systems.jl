export AbstractSystem, jacobian

"""
    AbstractSystem

An abstract type representing a polynomial system ``F(x)``.

The following systems are available:

* [`AffineChartSystem`](@ref)
* [`CompositionSystem`](@ref)
* [`FixedParameterSystem`](@ref)
* [`ModelKitSystem`](@ref)
* [`RandomizedSystem`](@ref)
"""
abstract type AbstractSystem end

# Base.size(F::AbstractSystem)
Base.size(F::AbstractSystem, i::Integer) = size(F)[i]

# ModelKit.variables(F::AbstractSystem)
ModelKit.parameters(F::AbstractSystem) = Variable[]
ModelKit.variable_groups(::AbstractSystem) = nothing

ModelKit.nvariables(F::AbstractSystem) = size(F, 2)
ModelKit.nparameters(F::AbstractSystem) = length(parameters(F))
function ModelKit.System(F::AbstractSystem)
    x, p = variables(F), parameters(F)
    System(F(x, p), x, p, variable_groups(F))
end

(F::AbstractSystem)(x, p = nothing) = evaluate(F, x, p)

"""
    evaluate(F::AbstractSystem, x, p = nothing)

Evaluate the given system.
"""
function ModelKit.evaluate(F::AbstractSystem, x, p = nothing)
    u = Vector{Any}(undef, size(F, 1))
    evaluate!(u, F, x, p)
    ModelKit.to_smallest_eltype(u)
end

"""
    jacobian(F::AbstractSystem, x, p = nothing)

Compute the Jacobian of the given system.
"""
function jacobian(F::AbstractSystem, x, p = nothing)
    u = Vector{Any}(undef, size(F, 1))
    U = Matrix{Any}(undef, size(F))
    evaluate_and_jacobian!(u, U, F, x, p)
    ModelKit.to_smallest_eltype(U)
end

include("systems/model_kit_system.jl")
include("systems/fixed_parameter_system.jl")
include("systems/randomized_system.jl")
include("systems/affine_chart_system.jl")
include("systems/composition_system.jl")

include("systems/start_pair_system.jl")
