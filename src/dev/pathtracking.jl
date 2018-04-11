module PathTracking

import ..PathTrackers: AbstractPathTrackerMethod, AbstractPathTrackerCache, AbstractPathTrackerState, Projective
import ..PathTrackers
import ..NewHomotopies: AbstractHomotopy

export Projective,
     PathTracker,
     PathTrackerResult,
     track,
     # advanced
     track!,
     current_t,
     current_x,
     current_status
"""
     PathTracker(H::NewHomotopies.AbstractHomotopy, x₀, t₁, t₀; options=Options(), method=Projective, method_options...)::PathTracker

Create a `PathTracker` to track `x₀` from `t₁` to `t₀`. The homotopy `H` needs to be
compatible with the chosen path tracker method `tracker`. The `method_options` depend
on the chosen path tracker for the default case see [`Projective`](@ref).

     PathTracker(method::AbstractPathTrackerMethod, x₀, t₁, t₀, options)

If a method is already assembled this constructor is beneficial.
"""
struct PathTracker{Method<:AbstractPathTrackerMethod, S<:AbstractPathTrackerState, C<:AbstractPathTrackerCache}
     method::Method
     options::PathTrackers.Options
     # TODO: The following actually depend on the current precision. So we would need to introduce
     # multiple of these and switch if necessary.
     state::S
     cache::C
end
function PathTracker(H::AbstractHomotopy, x, start, target; options=PathTrackers.Options(), method=Projective, kwargs...)
    PathTracker(method(H; kwargs...), x, start, target, options)
end

function PathTracker(method::AbstractPathTrackerMethod, x, start, target, options=PathTrackers.Options())
    tracker_state = PathTrackers.state(method, x, start, target)
    tracker_cache = PathTrackers.cache(method, tracker_state)
    PathTracker(method, options, tracker_state, tracker_cache)
end

"""
     PathTrackerResult{T}

Containing the result of a tracked path. The fields are
* `successfull::Bool` Indicating whether tracking was successfull.
* `returncode::Symbol` A return code, which gives an indication what
 happened if the tracking was not successfull
* `x::Vector{T}` The result.
* `t::Float64` The `t` when the path tracker stopped.
* `res::Float64` The residual at `(x, t)`.
"""
struct PathTrackerResult{T}
     successfull::Bool
     returncode::Symbol
     x::Vector{T}
     t::Float64
     res::Float64
end

"""
    track(tracker, x₀, t₁, t₀)::PathTrackerResult

Track a value `x₀` from `t₁` to `t₀` using the given `PathTracker` `tracker`.
This returns a `PathTrackerResult`. This modifies tracker.
"""
function track(tracker::PathTracker, x₀, t₁, t₀)
     PathTrackers.reset!(tracker.state, tracker.method, tracker.cache, x₀, t₁, t₀)
     track!(tracker)
     result(tracker)
end

"""

    track!(tracker)

Run the given `tracker` with the set values. This is useful to call directly
after construction of a `PathTracker`.
"""
function track!(tracker::PathTracker)
     method, state, cache, options = tracker.method, tracker.state, tracker.cache, tracker.options
     while PathTrackers.current_status(state) == :ok
          PathTrackers.step!(state, method, cache, options)
          PathTrackers.check_terminated!(state, method, cache, options)
     end
     nothing
end

"""
     result(tracker)

Obtain a result from a path tracker. This will also do a refinement step.
"""
function result(tracker::PathTracker)
     returncode = current_status(tracker)
     successfull = returncode == :success
     if successfull
          res = convert(Float64, PathTrackers.refine!(tracker.state, tracker.method, tracker.cache, tracker.options))
     else
          res = NaN
     end
     x = copy(current_x(tracker))
     PathTrackerResult(successfull, returncode, x, current_t(tracker), res)
end


"""
    current_t(tracker::PathTracker)

Get the current `t` from the tracker.
"""
current_t(tracker) = PathTrackers.current_t(tracker.state)

"""
    current_x(tracker::PathTracker)

Get the current solution `x` from the tracker. Note this is *not* a copy. If
you want to store it persistent you need to make a copy.
"""
current_x(tracker) = PathTrackers.current_x(tracker.state)

"""
    current_status(tracker::PathTracker)

Get the current status from the tracker.
"""
current_status(tracker) = PathTrackers.current_status(tracker.state)

# Iterator interface
Base.start(::PathTracker) = 0
function Base.next(tracker::PathTracker, k)
     method, state, cache, options = tracker.method, tracker.state, tracker.cache, tracker.options
     PathTrackers.step!(state, method, cache, options)
     PathTrackers.check_terminated!(state, method, cache, options)
     tracker, k + 1
end
function Base.done(tracker::PathTracker, k)
     PathTrackers.current_status(tracker.state) != :ok
end
Base.eltype(tracker::T) where {T<:PathTracker} = T
Base.iteratorsize(tracker::PathTracker) = Base.SizeUnknown()

end
