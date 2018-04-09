"""
     pathtracker(H::NewHomotopies.AbstractHomotopy, x₀, t₁, t₀; tracker=ProjectiveTracker, options...)::TrackerIterator

Create a path tracker to track `x₀` from `t₁` to `t₀`. The homotopy `H` needs to be
compatible with the chosen path tracker `tracker`. The options depend
on the chosen path tracker for the default case see [`ProjectiveTracker`](@ref).
"""
function pathtracker end

"""
     iterator(::AbstractPathTracker, x₀, t₁, t₀)

Construct an iterator to track `x₀` from `t₁` to `t₀`.
"""
function iterator end

"""
    Options(;tolerance=1e-7, maxiters=10_000)

The options used for the pathtracker.
* `tolerance`: The precision used to track a value.
* `maxiters`: The maximal number of iterations.
"""
struct Options
    tolerance::Float64
    maxiters::Int64
end
function Options(;tolerance=1e-7, maxiters=10_000)
    Options(tolerance, maxiters)
end

struct TrackerIterator{Tracker<:AbstractPathTracker, S<:AbstractPathTrackerState, C<:AbstractPathTrackerCache}
     tracker::Tracker
     # TODO: The following actually depend on the current precision. So we would need to introduce
     # multiple of these and switch if necessary.
     state::S
     cache::C
end


function setup!(iter::TrackerIterator, x₀, t₁, t₀)
     reset!(iter.state, iter.tracker, x₀, t₁, t₀)
     iter
end

function track!(iter::TrackerIterator)
     tracker, state, cache = iter.tracker, iter.state, iter.cache
     while !isdone(tracker, state, cache)
          step!(tracker, state, cache)
     end
     iter
end

function track!(iter::TrackerIterator, x₀, t₁, t₀)
     setup!(iter, x₀, t₁, t₀)
     track!(iter)
end
# Iterator interface
Base.start(::TrackerIterator) = 0
function Base.next(iter::TrackerIterator, k)
     step!(iter.tracker, iter.state, iter.cache)
     iter, k + 1
end
Base.done(iter::TrackerIterator, k) = isdone(iter.tracker, iter.state, iter.cache)
Base.eltype(iter::T) where {T<:TrackerIterator} = T
Base.iteratorsize(iter::TrackerIterator) = Base.SizeUnknown()
