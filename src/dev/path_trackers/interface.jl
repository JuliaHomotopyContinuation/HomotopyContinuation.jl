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
