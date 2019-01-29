export AbstractSystem, AbstractSystemCache,
    AbstractHomotopy, AbstractHomotopyCache

"""
    AbstractSystem

Representing a system of polynomials.
"""
abstract type AbstractSystem end

"""
    AbstractSystemCache

A cache to avoid allocations for the evaluation of an [`AbstractSystem`](@ref).
"""
abstract type AbstractSystemCache end


"""
    AbstractHomotopy

Representing a homotopy.
"""
abstract type AbstractHomotopy end


"""
    AbstractHomotopyCache

A cache to avoid allocations for the evaluation of an [`AbstractHomotopy`](@ref).
"""
abstract type AbstractHomotopyCache end

include("systems.jl")
include("homotopies.jl")
