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
    truncate_infinity::Bool
end

function Options(;sampling_factor=0.5, tol=1e-10, minradius=1e-15, maxnorm=1e5,
        maxwindingnumber=12, max_extrapolation_samples=4,
        truncate_infinity=true)
    Options(sampling_factor, tol, minradius, maxnorm, maxwindingnumber,
        max_extrapolation_samples, truncate_infinity)
end

mutable struct State{V<:ProjectiveVectors.AbstractProjectiveVector, C, T}
    # Current value
    x::V
    # prediction and previous prediction
    p::V
    pprev::V
    pbest::V
    pbest_delta::Float64
    # Samples of the solution path obtained from a geometric series
    # defined by x(λ^k R₀).
    # The raw vectors are stored.
    samples::Matrix{C}
    nsamples::Int
    # The data structure is that `logabs_samples` Contains a vector for each
    # coordinate.
    logabs_samples::Matrix{T}
    directions::Matrix{T}

    R::Float64
    npredictions::Int
    iters::Int
    # current status
    status::Symbol
    # current estimate of winding number
    windingnumber::Int
    in_endgame_zone::Bool
end

function State(x, R₀, options)
    n = maxsamples(R₀, options)

    R = real(R₀)

    samples = fill(zero(eltype(x)), length(x), n)
    nsamples = 1
    directions = fill(NaN, length(x), n)
    logabs_samples = copy(directions)

    p = copy(x)
    pprev = copy(x)
    pbest = copy(x)
    pbest_delta = Inf

    npredictions = 0
    iters = 0

    status = :ok
    windingnumber = 1
    in_endgame_zone = false

    State(copy(x), p, pprev, pbest, pbest_delta, samples, nsamples, logabs_samples,
        directions,
        R, npredictions, iters, status,
        windingnumber, in_endgame_zone)
end

function maxsamples(R₀, options)
    ceil(Int, log(options.sampling_factor, options.minradius / real(R₀))) + 2
end


function reset!(state::State, x, R)
    state.nsamples = 1
    state.x .= x
    for i = 1:length(x)
        state.samples[i, 1] = x[i]
        state.logabs_samples[i, 1] = logabs(x[i])
    end

    state.pbest_delta = Inf

    state.npredictions = 0
    state.iters = 0
    state.R = R
    state.status = :ok

    state.windingnumber = 0
    state.in_endgame_zone = false
    nothing
end

struct Cache{V}
    log_sampling_factor::Float64
    windingnumbers::Vector{Int}
    unitroots::Vector{Complex128}
    direction_buffer::Matrix{Float64}
    xbuffer::V
    pbuffer::V
end

function Cache(state::State, options::Options)
    log_sampling_factor = log(options.sampling_factor)
    windingnumbers = zeros(Int, length(state.x))
    unitroots = Vector{Complex128}()
    direction_buffer = zeros(length(state.x), options.max_extrapolation_samples)
    pbuffer = copy(state.x)
    xbuffer = copy(state.x)
    Cache(log_sampling_factor, windingnumbers, unitroots, direction_buffer,
        pbuffer, xbuffer)
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
    state = State(x, R₀, options)
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
    res::Float64
    windingnumber::Int
    npredictions::Int
    iters::Int
end

function EndgamerResult(endgamer::Endgamer)
    state = endgamer.state

    if state.npredictions == 0 || state.status == :at_infinity
        x = copy(state.x)
    else
        x = copy(state.pbest)
    end
    EndgamerResult(state.status, x, real(state.R), state.pbest_delta,
        state.windingnumber, state.npredictions, state.iters)
end
