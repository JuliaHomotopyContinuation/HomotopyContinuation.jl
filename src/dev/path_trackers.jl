module PathTrackers

export AbstractPathTracker,
    AbstractPathTrackerCache,
    AbstractPathTrackerState,
    Options,
    state,
    cache,
    reset!,
    step!,
    isdone,
    current_t,
    current_x,
    current_status

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
    Options(;options..)

The possible options used for the pathtracker are:
* `tol=1e-7`: The precision used to track a value.
* `corrector_maxiters=3`: The maximal number of correction steps in a single step.
* `refinement_tol=max(1e-15, tolerance^2)`: The precision used to refine the final value.
* `refinement_maxiters=3`: The maximal number of correction steps used to refine the final value.
* `maxiters=10_000`: The maximal number of iterations.
"""
mutable struct Options
    tol::Float64
    refinement_tol::Float64
    refinement_maxiters::Float64
    corrector_maxiters::Int
    maxiters::Int
end
function Options(;tol=1e-7, refinement_tol=max(1e-15, tol^2), corrector_maxiters=3, refinement_maxiters=20, maxiters=10_000)
    Options(tol, refinement_tol, refinement_maxiters, corrector_maxiters, maxiters)
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
    step!(::AbstractPathTrackerState, ::AbstractPathTrackerMethod, AbstractPathTrackerCache, options::Options)

Make one step in the path tracking method.
"""
function step! end

"""
    check_terminated!(state::AbstractPathTrackerState, ::AbstractPathTrackerMethod, AbstractPathTrackerCache, options::Options)

Check whether the path tracking terminated and update `status` in `state` accordingly.
If the path reached its target set the state to `:success`. Every state other
than `:ok` will be interpreted as a termination of the algorithm.
"""
function check_terminated! end

# Access
"""
    current_t(state::AbstractPathTrackerState)

Get the current absolute `t` from the state `state`.
"""
function current_t end

"""
    current_x(state::AbstractPathTrackerState)

Get the current solution `x` from the state `state`.
"""
function current_x end

"""
    current_status(state::AbstractPathTrackerState)

Get the current status from the state `state`.
"""
function current_status end

# Implementation
include("path_trackers/projective.jl")

end
