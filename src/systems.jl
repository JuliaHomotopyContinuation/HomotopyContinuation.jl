export AbstractSystem

"""
    AbstractSystem

An abstract type representing a polynomial system ``F(x)``.

The following systems are available:

* [`ModelKitSystem`](@ref)
* [`RandomizedSystem`](@ref)
"""
abstract type AbstractSystem end

include("systems/model_kit_system.jl")
include("systems/randomized_system.jl")
