"""
    Options(;options)

## Options
* `sampling_factor=0.3`
* `tol=1e-12`
* `minradius=1e-20`
* `maxnorm=1e6`
"""
struct Options
    sampling_factor::Float64
    tol::Float64
    minradius::Float64
    maxnorm::Float64
    maxwindingnumber::Float64
    max_extrapolation_samples::Int
end

function Options(;sampling_factor=0.6, tol=1e-10, minradius=1e-15, maxnorm=1e5, maxwindingnumber=12, max_extrapolation_samples=5)
    Options(sampling_factor, tol, minradius, maxnorm, maxwindingnumber, max_extrapolation_samples)
end

mutable struct State{V, T}
    # Samples of the solution path obtained from a geometric series
    # defined by x(λ^k R₀).
    samples::Vector{V}
    nsamples::Int
    # The data structure is that `logabs_samples` Contains a vector for each
    # coordinate.
    logabs_samples::Vector{Vector{T}}
    directions::Vector{Float64}
    directions_matching::Vector{Int}
    n_directions_matching::Int

    R::Float64

    # prediction and previous prediction
    p::V
    pprev::V
    # number of predictions so far
    npredictions::Int


    # number of iterations so far
    iters::Int

    # current status
    status::Symbol

    # current estimate of winding number
    windingnumber::Int
    # consecutive matching estimates
    cons_matching_estimates::Int
end

function State(x, R₀::Float64)
    samples = [copy(x)]
    nsamples = 1
    logabs_samples = [[logabs(x[i])] for i=1:length(x)]
    directions = fill(0.0, length(x))
    directions_matching = fill(0, length(x))
    n_directions_matching = 0

    p = copy(x)
    pprev = copy(x)

    npredictions = 0
    iters = 0

    status = :ok
    windingnumber = 1
    cons_matching_estimates = 0
    State(samples, nsamples, logabs_samples,
        directions, directions_matching, n_directions_matching,
        R₀, p, pprev, npredictions, iters, status,
        windingnumber,
        cons_matching_estimates)
end

struct Cache{V}
    # windingnumber_extrapolation_buffer::Vector{Float64}
    windingnumbers::Vector{Int}
    direction_buffer::Vector{Float64}
    fitpowerseries_buffer::Vector{V}
end

function Cache(state::State, options::Options)
    windingnumbers = zeros(Int, length(state.logabs_samples))
    direction_buffer = zeros(options.max_extrapolation_samples)
    fitpowerseries_buffer = [deepcopy(state.samples[1]) for _=1:options.max_extrapolation_samples]
    Cache(windingnumbers, direction_buffer,
        fitpowerseries_buffer)
end

"""
    Endgamer(endgame, pathtracker, R₀=0.1, options=EndgameOptions())

"""
struct Endgamer{P<:PathTracking.PathTracker, V}
    tracker::P
    state::State{V}
    cache::Cache
    options::Options
end

function Endgamer(tracker::PathTracking.PathTracker, R₀=0.1; options::Options=Options())
    x = PathTracking.currx(tracker)
    state = State(x, R₀)
    Endgamer(tracker, state, Cache(state, options), options)
end


"""
    EndgamerResult(endgamer)

## Fields
* `returncode::Symbol`
* `x::Vector{T}`: The solution or last prediction.
* `npredictions`: The number of predictions.
* `iters`: The number of iterations.
"""
struct EndgamerResult{V}
    returncode::Symbol
    x::V
    t::Float64
    windingnumber::Int
    npredictions::Int
    iters::Int
end

function EndgamerResult(endgamer::Endgamer)
    S = endgamer.state

    if S.npredictions == 0
        x = copy(S.samples[S.nsamples])
    else
        x = copy(S.p)
    end
    EndgamerResult(S.status, x, S.R, S.windingnumber, S.npredictions, S.iters)
end
