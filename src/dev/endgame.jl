module Endgame

import ..PathTracking
using ..Utilities


abstract type AbstractEndgame end

struct EndgameOptions
    λ::Float64
    tol::Float64
    minradius::Float64
    maxnorm::Float64
end

function EndgameOptions(;geomseries_factor=0.3, tol=1e-12, minradius=1e-20, maxnorm=1e6)
    EndgameOptions(geomseries_factor, tol, minradius, maxnorm)
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
    state.windingnumber_estimate = 1
end


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

function play!(endgamer::Endgamer, x, t)
    reset!(endgamer.state, x, t)
    reset!(endgamer.alg_state, endgamer.alg, x)

    state, options = endgamer.state, endgamer.options

    moveforward!(endgamer)
    moveforward!(endgamer)
    while state.status == :ok
        if state.npredictions > 0 && state.windingnumber_estimate == 1
            # if we have one prediction we also have an estimate for the winding number
            # since the heuristics ensure that we are in the endgame convergence zone
            try_to_jump_to_target!(endgamer)
        end
        predict!(endgamer)
        checkconvergence!(endgamer)
        moveforward!(endgamer)
    end
    endgamer
end

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
    elseif retcode ≠ :heuristic_failed
        state.status = retcode
    end

    return nothing
end

@inline function checkconvergence!(endgamer::Endgamer)
    endgamer.state.status == :ok || return nothing

    state = endgamer.state
    if state.npredictions > 1
        Δ = infinity_norm(state.p, state.pprev)
        @show Δ, state.p, state.pprev
        if Δ < endgamer.options.tol
            state.status = :success
        end
    end
end


@inline function moveforward!(endgamer::Endgamer)
    endgamer.state.status == :ok || return nothing

    state = endgamer.state

    λ = endgamer.options.λ
    state.iter += 1

    λR = λ * state.R
    if λR < endgamer.options.minradius
        state.status = :minradius_reached
    end

    x = state.xs[3]
    retcode = PathTracking.track!(x, endgamer.tracker, state.xs[1], state.R, λ * state.R)
    normalize!(x)
    state.xs = (x, state.xs[1], state.xs[2])

    if retcode != :success
        state.status = :ill_conditioned_zone
        return nothing
    end

    if infinity_norm(x) > endgamer.options.maxnorm
        state.status = :at_infinity
    end

    state.R = λR
    nothing
end

"""
    try_to_jump_to_target!(endgamer)

We try to reach the target without using an endgame. If the path tracking failed
we simply proceed as nothing was. This works well when we actually don't
need an endgame at all.
"""
@inline function try_to_jump_to_target!(endgamer::Endgamer)
    state = endgamer.state

    if state.status != :ok
       return nothing
   end

    retcode = PathTracking.track!(endgamer.tracker, state.xs[1], state.R, 0.0)
    if retcode == :success
        state.p .= PathTracking.current_x(endgamer.tracker)
        state.status = :success
        state.npredictions += 1
    end
    nothing
end
end
