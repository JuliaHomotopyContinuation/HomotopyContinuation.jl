module PathTrackers

export AbstractPathTracker,
    AbstractPathTrackerCache,
    AbstractPathTrackerState,
    Options

"""
     AbstractPathTrackerMethod

A `PathTrackerConfig` holds informations like the homotopy, algorithms used etc.
"""
abstract type AbstractPathTrackerMethod end

"""
     AbstractPathTrackerState

A `PathTrackerState` holds the current state. This is information which is reseted
after each run.
"""
abstract type AbstractPathTrackerState end

"""
     AbstractPathTrackerCache

A `PathTrackerCache` holds all data structures to avoid allocations.
"""
abstract type AbstractPathTrackerCache end

"""
    Options(;tolerance=1e-7, maxiters=10_000)

The options used for the pathtracker.
* `tolerance`: The precision used to track a value.
* `maxiters`: The maximal number of trackerations.
"""
mutable struct Options
    tolerance::Float64
    maxiters::Int64
end
function Options(;tolerance=1e-7, maxiters=10_000)
    Options(tolerance, maxiters)
end


# INITIAL CONSTRUCTION

"""
    state(method::AbstractPathTrackerMethod, x₀, t₁, t₀)::AbstractPathTrackerState

Initialize the necessary state to track a path from `t₁` to `t₀` starting
from `x₀`.
"""
function state end


"""
    cache(method::AbstractPathTrackerMethod, state::AbstractPathTrackerState)::AbstractPathTrackerCache

Initialize the cache to track paths.
"""
function cache end

# UPDATES

"""
    reset!(state::AbstractPathTrackerState, ::AbstractPathTrackerMethod, ::AbstractPathTrackerCache, x₀, t₁, t₀)

Reset the state to track a path from `t₁` to `t₀` starting from `x₀`.
"""
function reset! end

# ITERATION

"""
    step!(::AbstractPathTrackerMethod, ::AbstractPathTrackerState, AbstractPathTrackerCache, options::Options)

Make one step in the path tracking method.
"""
function step! end

"""
    isdone(::AbstractPathTrackerMethod, ::AbstractPathTrackerState, AbstractPathTrackerCache, options::Options)

Check whether the path tracking terminated.
"""
function isdone end

# Implementation
include("path_trackers/projective.jl")

end
