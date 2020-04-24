export AbstractSystem

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

Base.size(F::AbstractSystem, i::Integer) = size(F)[i]

include("systems/model_kit_system.jl")
include("systems/fixed_parameter_system.jl")
include("systems/randomized_system.jl")
include("systems/affine_chart_system.jl")
include("systems/composition_system.jl")
