using BenchmarkTools

mutable struct Tracker{T<:Real}
    x::Vector{Complex{T}}
    iter::Int
end

function Base.convert(::Type{Tracker{S}}, b::Tracker{T}) where {S,T}
    Tracker{S}(convert(Vector{Complex{S}}, b.x), b.iter)
end

mutable struct TrackerIterator{Low,High}
    low::Tracker{Low}
    high::Tracker{High}
    lowprecision::Bool
end

function TrackerIterator(low::Tracker{T}) where T
    high = Base.convert(Tracker{widen(T)}, low)
    TrackerIterator(low, high, true)
end

function sync!(iter::TrackerIterator)
    if lowprecision
        iter.high -= iter.low
    end
end

t = Tracker(rand(Complex64, 1000), 0)
titer = TrackerIterator(t)

@inline Base.start(::TrackerIterator) = 0
@inline Base.next(x::TrackerIterator, state) = x, state + 1
@inline Base.done(x::TrackerIterator, state) = state > 100

function step!(x::AbstractVector{T}) where T
    x .= x  .+ mean(x)
end

function solve!(tracker::TrackerIterator)
    while tracker.state < 1000
        tracker.state += 1
        step2!(tracker.x)
    end
end





#
# @inline Base.start(::AdaptiveTracker) = 0
# @inline Base.next(x::AdaptiveTracker, state) = x.x, state + 1
# @inline Base.done(x::AdaptiveTracker, state) = state > 1000

mutable struct FixedTracker{AV<:AbstractVector}
    x::AV
end

@inline Base.start(::FixedTracker) = 0
@inline Base.next(x::FixedTracker, state) = x, state + 1
@inline Base.done(x::FixedTracker, state) = state > 100





function solve!(tracker::FixedTracker)
    k = 0
    for t in tracker
        k += 1
        if k == 50
            convert.(Int, round.(Int, tracker.x))
            step!(t.x)
        else
            step!(t.x)
        end
    end
end

x = rand(Complex128, 10000)
tracker = AdaptiveTracker(x, 0)
solve2!(tracker)


@benchmark solve!($tracker)


fixed = FixedTracker(x)
@btime solve!($fixed)





x = rand(4)

t = AdaptiveTracker(x)

t.x = rand(Int, 4)


a, b = rand(Complex64, 2)



@benchmark *($a, $b)
