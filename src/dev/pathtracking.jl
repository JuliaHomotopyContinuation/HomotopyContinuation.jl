module PathTracking

import ..PathTrackers: AbstractPathTrackerMethod, AbstractPathTrackerCache, AbstractPathTrackerState, Projective
import ..PathTrackers
import ..NewHomotopies: AbstractHomotopy

export Projective,
     PathTracker,
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
     # TODO: The following actually depend on the current precision. So we would need to introduce
     # multiple of these and switch if necessary.
     state::S
     cache::C
     options::PathTrackers.Options
end
function PathTracker(H::AbstractHomotopy, x, start, target; options=PathTrackers.Options(), method=Projective, kwargs...)
    PathTracker(method(H; kwargs...), x, start, target, options)
end

function PathTracker(method::AbstractPathTrackerMethod, x, start, target, options=PathTrackers.Options())
    tracker_state = PathTrackers.state(method, x, start, target)
    tracker_cache = PathTrackers.cache(method, tracker_state)
    PathTracker(method, tracker_state, tracker_cache, options)
end

"""
    track!(tracker, x₀, t₁, t₀)

Track a value `x₀` from `t₁` to `t₀` using the given `PathTracker` `tracker`.

    track!(tracker)

Run the given `tracker` with the set values. This is useful to call directly
after construction of a `PathTracker`.
"""
function track!(tracker::PathTracker)
     while !PathTrackers.isdone(tracker.method, tracker.state, tracker.cache, tracker.options)
          PathTrackers.step!(tracker.method, tracker.state, tracker.cache, tracker.options)
     end
     tracker
end

function track!(tracker::PathTracker, x₀, t₁, t₀)
     PathTrackers.reset!(tracker.state, tracker.method, tracker.cache, x₀, t₁, t₀)
     track!(tracker)
end

"""
    current_t(tracker::PathTracker)

Get the current `t` from the tracker.
"""
current_t(tracker) = PathTrackers.current_t(tracker.state)

"""
    current_x(tracker::PathTracker)

Get the current solution `x` from the tracker.
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
     PathTrackers.step!(tracker.method, tracker.state, tracker.cache, tracker.options)
     tracker, k + 1
end
function Base.done(tracker::PathTracker, k)
     PathTrackers.isdone(tracker.method, tracker.state, tracker.cache, tracker.options)
end
Base.eltype(tracker::T) where {T<:PathTracker} = T
Base.iteratorsize(tracker::PathTracker) = Base.SizeUnknown()

end
