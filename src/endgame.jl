module Endgame

import ..PathTracking
using ..Utilities


abstract type AbstractEndgame end

"""
    EndgameOptions

Options:
* `geomseries_factor=0.3`
* `tol=1e-12`
* `minradius=1e-20`
* `maxnorm=1e6`
"""
struct EndgameOptions
    λ::Float64
    tol::Float64
    minradius::Float64
    maxnorm::Float64
    maxwindingnumber::Float64
end

function EndgameOptions(;geomseries_factor=0.3, tol=1e-10, minradius=1e-20, maxnorm=1e5, maxwindingnumber=12)
    EndgameOptions(geomseries_factor, tol, minradius, maxnorm, maxwindingnumber)
end

include("endgames/cauchy.jl")

mutable struct EndgamerState{T}
    xs::NTuple{3, Vector{T}}
    R::Float64

    p::Vector{T}
    pprev::Vector{T}

    npredictions::Int
    iter::Int

    status::Symbol
    windingnumber_estimate::Int
end

function EndgamerState(x, R::Float64)
    xs = (copy(x), copy(x), copy(x))

    p = copy(x)
    pprev = copy(x)

    npredictions = 0
    iter = 0

    status = :ok
    windingnumber_estimate = 1

    EndgamerState(xs, R, p, pprev, npredictions, iter, status, windingnumber_estimate)
end

function reset!(state::EndgamerState, x, R::Float64)
    state.xs[1] .= x
    state.npredictions = 0
    state.iter = 0
    state.R = R
    state.status = :ok
    state.windingnumber_estimate = 0
end


"""
    Endgamer(endgame, pathtracker, t=0.1, options=EndgameOptions())

"""
struct Endgamer{E<:AbstractEndgame, EC, ES, P<:PathTracking.PathTracker, T}
    alg::E
    alg_cache::EC
    alg_state::ES
    tracker::P
    state::EndgamerState{T}
    options::EndgameOptions
end


function Endgamer(alg::AbstractEndgame, tracker::PathTracking.PathTracker, t=0.1, options=EndgameOptions())
    x = PathTracking.current_x(tracker)
    alg_cache = cache(alg, x)
    alg_state = state(alg, x)

    Endgamer(alg, alg_cache, alg_state, tracker, EndgamerState(x, t), options)
end


"""
    EndgamerResult(endgamer)

## Fields
* `returncode::Symbol`
* `x::Vector{T}`: The solution or last prediction.
* `npredictions`: The number of predictions.
* `iters`: The number of iterations.
"""
struct EndgamerResult{T}
    returncode::Symbol
    x::Vector{T}
    t::Float64
    windingnumber_estimate::Int
    npredictions::Int
    iters::Int
end

function EndgamerResult(endgamer::Endgamer)
    S = endgamer.state

    EndgamerResult(S.status, copy(S.p), S.R, S.windingnumber_estimate, S.npredictions, S.iter)
end


function play(endgamer::Endgamer, x, t)
    play!(endgamer, x, t)
    EndgamerResult(endgamer)
end

"""
    play!(endgamer::Endgamer, x, t)

Run the `Endgamer` `endgamer` starting from `(x, t)`.
"""
function play!(endgamer::Endgamer, x, t)
    reset!(endgamer.state, x, t)
    reset!(endgamer.alg_state, endgamer.alg, x)

    state, options = endgamer.state, endgamer.options

    moveforward!(endgamer)
    moveforward!(endgamer)
    while state.status == :ok
        # if state.npredictions > 0 && state.windingnumber_estimate == 1
        #     # if we have one prediction we also have an estimate for the winding number
        #     # since the heuristics ensure that we are in the endgame convergence zone
        #     try_to_jump_to_target!(endgamer)
        # end
        predict!(endgamer)
        checkconvergence!(endgamer)
        moveforward!(endgamer)
    end
    endgamer
end

"""
    predict!(endgamer::Endgamer)

Try to predict the value of the path `x(t)` and `t=0.0`.
"""
@inline function predict!(endgamer::Endgamer)
    endgamer.state.status == :ok || return nothing

    state = endgamer.state
    p = state.p
    if state.npredictions > 0
        state.pprev .= p
    end
    retcode, windingnumber = predict!(p, endgamer.alg, endgamer.alg_state,
        endgamer.alg_cache, endgamer.tracker, state.xs, state.R, state.npredictions, endgamer.options)
    if retcode == :success
        state.npredictions += 1
        state.windingnumber_estimate = windingnumber
        if infinity_norm(p, 1) > endgamer.options.maxnorm
            endgamer.state.status = :at_infinity
        end
    elseif retcode ≠ :heuristic_failed
        state.status = retcode
    end

    return nothing
end

"""
    checkconvergence!(endgamer)

Try to check whether the endgame is converged.
"""
@inline function checkconvergence!(endgamer::Endgamer)
    endgamer.state.status == :ok || return nothing

    state = endgamer.state
    if state.npredictions > 1
        Δ = infinity_norm(state.p, state.pprev)
        if Δ < endgamer.options.tol
            state.R = 0.0
            state.status = :success
        end
    end
end


"""
    moveforward!(endgamer)

Move the nearer to `0.0` by using the `geomseries_factor` of `EndgameOptions`.
"""
function moveforward!(endgamer::Endgamer)
    endgamer.state.status == :ok || return nothing

    state = endgamer.state

    state.iter += 1

    λ = endgamer.options.λ
    λR = λ * state.R
    if λR < endgamer.options.minradius
        state.status = :minradius_reached
    end

    x = state.xs[3]
    retcode = PathTracking.track!(x, endgamer.tracker, state.xs[1], state.R, λ * state.R)
    if retcode != :success
        state.status = :ill_conditioned_zone
        return nothing
    end
    if infinity_norm(x, 1) > endgamer.options.maxnorm
        state.status = :at_infinity
        state.p .= x
    end

    normalize!(x)
    state.xs = (x, state.xs[1], state.xs[2])

    state.R = λR

    nothing
end

"""
    try_to_jump_to_target!(endgamer)

We try to reach the target without using an endgame. If the path tracking failed
we simply proceed as nothing was. This works well when we actually don't
need an endgame at all.
"""
function try_to_jump_to_target!(endgamer::Endgamer)
    state = endgamer.state

    if state.status != :ok
       return nothing
   end

    retcode = PathTracking.track!(endgamer.tracker, state.xs[1], state.R, 0.0)
    if retcode == :success
        state.p .= PathTracking.current_x(endgamer.tracker)
        state.R = 0.0
        state.status = :success
        state.npredictions += 1
    end
    nothing
end
end
