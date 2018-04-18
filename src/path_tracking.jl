module PathTracking

import ..PathTrackers: AbstractPathTrackerMethod, AbstractPathTrackerCache, AbstractPathTrackerState, Projective, Options
import ..PathTrackers
import ..Homotopies: AbstractHomotopy
import ..Problems
using ..Utilities

export Projective,
     Options,
     PathTracker,
     tol, refinement_tol, refinement_maxiters, corrector_maxiters,
     set_tol!, set_refinement_tol!, set_refinement_maxiters!, set_corrector_maxiters!,
     PathTrackerResult,
     track,
     # advanced
     track!,
     current_t,
     current_x,
     current_iter,
     current_status
"""
     PathTracker(H::Homotopies.AbstractHomotopy, x₁, t₁, t₀; options=Options(), method=Projective, method_options...)::PathTracker

Create a `PathTracker` to track `x₁` from `t₁` to `t₀`. The homotopy `H` needs to be
compatible with the chosen path tracker method `tracker`. The `method_options` depend
on the chosen path tracker for the default case see [`Projective`](@ref).

     PathTracker(method::AbstractPathTrackerMethod, x₁, t₁, t₀, options)

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
function PathTracker(prob::Problems.AbstractProblem, x₁::AbstractVector{<:Number}, t₁, t₀; kwargs...)
     PathTracker(prob.homotopy, Problems.embed(prob, x₁), t₁, t₀; kwargs...)
end
function PathTracker(H::AbstractHomotopy, x, start, target; options=PathTrackers.Options(), method=Projective, kwargs...)
    PathTracker(method(H; kwargs...), x, start, target, options)
end

function PathTracker(method::AbstractPathTrackerMethod, x, start, target, options=PathTrackers.Options())
    tracker_state = PathTrackers.state(method, x, start, target)
    tracker_cache = PathTrackers.cache(method, tracker_state)
    PathTracker(method, options, tracker_state, tracker_cache)
end

# Access + modifications
"""
     tol(tracker::PathTracker)

Current tolerance.
"""
tol(tracker::PathTracker) = tracker.options.tol

"""
     refinement_tol(tracker::PathTracker)

Current refinement tolerance.
"""
refinement_tol(tracker::PathTracker) = tracker.options.refinement_tol

"""
     refinement_maxiters(tracker::PathTracker)

Current refinement maxiters.
"""
refinement_maxiters(tracker::PathTracker) = tracker.options.refinement_maxiters

"""
     corrector_maxiters(tracker::PathTracker)

Current correction maxiters.
"""
corrector_maxiters(tracker::PathTracker) = tracker.options.corrector_maxiters

"""
     set_tol!(tracker::PathTracker, tol)

Set the current tolerance to `tol`.
"""
function set_tol!(tracker::PathTracker, tol)
     tracker.options.tol = tol
     tol
end

"""
     set_refinement_maxiters!(tracker::PathTracker, tol)

Set the current refinement tolerance to `tol`.
"""
function set_refinement_tol!(tracker::PathTracker, tol)
     tracker.options.refinement_tol = tol
     tol
end

"""
     set_refinement_maxiters!(tracker::PathTracker, n)

Set the current refinement maxiters to `n`.
"""
function set_refinement_maxiters!(tracker::PathTracker, n)
     tracker.options.refinement_maxiters = n
     n
end

"""
     set_corrector_maxiters!(tracker::PathTracker, n)

Set the current correction maxiters to `n`.
"""
function set_corrector_maxiters!(tracker::PathTracker, n)
     tracker.options.corrector_maxiters = n
     n
end


"""
     PathTrackerResult{V<:AbstractVector}

Containing the result of a tracked path. The fields are
* `successfull::Bool` Indicating whether tracking was successfull.
* `returncode::Symbol` If the tracking was successfull then it is `:success`.
Otherwise the return code gives an indication what happened.
* `x::V` The result.
* `t::Float64` The `t` when the path tracker stopped.
* `res::Float64` The residual at `(x, t)`.
"""
struct PathTrackerResult{T, V<:AbstractVector}
     successfull::Bool
     returncode::Symbol
     x::V
     t::T
     res::Float64
     iters::Int
end

"""
    track(tracker, x₁, t₁, t₀)::PathTrackerResult

Track a value `x₁` from `t₁` to `t₀` using the given `PathTracker` `tracker`.
This returns a `PathTrackerResult`. This modifies tracker.

    track(tracker, xs, t₁, t₀)::Vector{PathTrackerResult}

Track all values in xs from `t₁` to `t₀` using the given `PathTracker` `tracker`.
"""
function track(tracker::PathTracker, x₁::AbstractVector{<:Number}, t₁, t₀)
     PathTrackers.reset!(tracker.state, tracker.method, tracker.cache, x₁, t₁, t₀)
     track!(tracker)
     result(tracker)
end
function track(tracker::PathTracker, xs, t₁, t₀)
     map(xs) do x₁
          PathTrackers.reset!(tracker.state, tracker.method, tracker.cache, x₁, t₁, t₀)
          track!(tracker)
          result(tracker)
     end
end



"""
    track!(tracker)

Run the given `tracker` with the set values. This is useful to call directly
after construction of a `PathTracker`. Returns a `Symbol` indicating the status.
If the tracking was successfull it is `:success`.

     track!(x₀, tracker, x₁, t₁, t₀ [, reset_steplength=true])

Track a value `x₁` from `t₁` to `t₀` using the given `PathTracker` `tracker`
and store the result in `x₀`. Returns a `Symbol` indicating the status.
If the tracking was successfull it is `:success`.
"""
function track!(tracker::PathTracker)
     method, state, cache, options = tracker.method, tracker.state, tracker.cache, tracker.options
     while PathTrackers.current_status(state) == :ok
          PathTrackers.step!(state, method, cache, options)
          PathTrackers.check_terminated!(state, method, cache, options)
     end
     current_status(tracker)
end
function track!(x₀, tracker::PathTracker, x₁, t₁, t₀)
     track!(tracker::PathTracker, x₁, t₁, t₀)
     x₀ .= current_x(tracker)
     current_status(tracker)
end
function track!(tracker::PathTracker, x₁, t₁, t₀)
     PathTrackers.reset!(tracker.state, tracker.method, tracker.cache, x₁, t₁, t₀)
     track!(tracker)
     returncode = current_status(tracker)
     if returncode == :success
          PathTrackers.refine!(tracker.state, tracker.method, tracker.cache, tracker.options)
     end
     current_status(tracker)
end

function reset!(tracker::PathTracker, x₁, t₁, t₀)
     PathTrackers.reset!(tracker.state, tracker.method, tracker.cache, x₁, t₁, t₀)
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
     PathTrackerResult(successfull, returncode, copy(tracker.state.x), current_t(tracker), res, current_iter(tracker))
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
    current_iter(tracker::PathTracker)

Get the number of the current iteration from the tracker.
"""
current_iter(tracker) = PathTrackers.current_iter(tracker.state)

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
