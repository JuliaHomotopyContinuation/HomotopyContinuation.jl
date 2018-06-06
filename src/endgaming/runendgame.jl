import ..ProjectiveVectors
import ..Homotopies

function setup!(endgamer, x, R)
    PathTracking.setup!(endgamer.tracker, x, R, 0.0im)
    reset!(endgamer.state, x, R)

    # The endgame fails if we already have a solution. Hence we check it now.
    r = PathTracking.residual(endgamer.tracker, x, 0.0)
    if abs(r) < endgamer.options.tol
        endgamer.state.pbest .= x
        endgamer.state.pbest_delta = r
        endgamer.state.status = :success
        endgamer.state.R = 0.0im
    end
    nothing
end

function runendgame(endgamer, x, R)
    setup!(endgamer, x, R)
    runendgame!(endgamer)
    Result(endgamer)
end

function runendgame!(endgamer)
    state, tracker, options, cache = endgamer.state, endgamer.tracker, endgamer.options, endgamer.cache
    while state.status == :ok
        nextsample!(tracker, state, options)
        update_derived_from_samples!(state, options, cache)
        if in_endgame_zone!(state)
            predict_infinity_check!(state, tracker, options, cache)
        end
        checkterminate!(state, endgamer.options)
    end

    endgamer
end


"""
    nextsample!(endgamer)

Obtain a new sample on the geometric series and add it to `endgamer`.
"""
function nextsample!(tracker::PathTracking.PathTracker, state, options)
    state.iters += 1

    R, λ = state.R, options.sampling_factor
    λR = λ * R

    retcode = PathTracking.track!(tracker, state.x, R, λR, precondition=false)
    if retcode != :success
        state.status = :tracker_failed
        return nothing
    end

    state.R = λR
    state.x .= PathTracking.currx(tracker)
    state.nsamples += 1
    for i=1:length(state.x)
        state.samples[i, state.nsamples] = state.x[i]
    end
    nothing
end

function update_derived_from_samples!(state, options, cache)
    # TODO: check whether we still need to collect samples
    for i=1:length(state.x)
        state.logabs_samples[i, state.nsamples] = logabs(state.x[i])
    end

    update_directions!(state, options, cache)
end

function in_endgame_zone!(state)
    if state.in_endgame_zone
        return true
    end

    if state.R < 1e-8
        state.in_endgame_zone = true
        return true
    end

    n, dirs = state.nsamples, state.directions

    if n < 5
        return false
    end

    for i=1:length(state.x)
        di_n = dirs[i, n]
        k = n - 1
        stable = true
        while k ≥ 2 && k ≥ n - 3 && stable
            if abs(dirs[i, k] - di_n) > 1e-2
                stable = false
            end
            k -= 1
        end
        if !stable
            return false
        end
    end
    state.in_endgame_zone = true
end

function analyze_finite_atinfinity(state, options, cache)
    checkfinite(state) && return :finite
    checkatinfinity(state, options) && return :at_infinity
    :unknown
end

checkfinite(state) = checkfinite(state, state.x)
function checkfinite(state, x::ProjectiveVectors.PVector)
    i₀ = ProjectiveVectors.homvar(x)
    d_i₀ = state.directions[i₀, state.nsamples]
    if abs(d_i₀) < 1e-3
        return true
    end

    for i=1:length(state.x)
        i == i₀ && continue
        if abs(state.directions[i, state.nsamples] - d_i₀) > 1e-3
            return false
        end
    end
    true
end

checkatinfinity(state, options) = checkatinfinity(state, options, state.x)
function checkatinfinity(state, options, x::ProjectiveVectors.PVector)
    LOOKBACK = 1
    TOL = 1e-3
    MIN_LOGABS = -4.605170185988091 # = log(0.01)

    i₀ = ProjectiveVectors.homvar(x)
    n = state.nsamples
    dirs = state.directions


    # We want the homogenization variable to be at least somewhat small
    if state.R > 1e-8 && state.logabs_samples[i₀, n] > MIN_LOGABS
        return false
    end
    for i=1:length(x)
        i == i₀ && continue
        di_n = dirs[i, n] - dirs[i₀, n]
        di_n > -0.032 && continue
        # Okay the current direction is negative now we check whether we
        # have at least some stability
        k = n - 1
        stable = true
        while k ≥ 2 && k ≥ n - LOOKBACK && stable
            if abs(dirs[i, k] - dirs[i₀, k] - di_n) > TOL
                stable = false
            end
            k -= 1
        end
        if stable
            return true
        end
    end
    false
end

function predict_infinity_check!(state, tracker, options, cache)
    # @show state.path_type

    path_type = analyze_finite_atinfinity(state, options, cache)
    if path_type == :finite
        # we need to determine the windingnumber
        if state.windingnumber == 0
            determine_windingnumber!(state, tracker, options, cache)
            # If we found that the windingnumber is one then
            # then we have a non-singular solution and can simply
            # track to the solution.
            if state.windingnumber == 1
                tracktoend!(state, tracker)
            end
            return
        else
            predict_cif!(state, tracker, options, cache)
            return
        end
    elseif (path_type == :at_infinity) ||
           path_type == :unknown && checkatinfinity_norm(state, options)
        state.status = :at_infinity
    end

    nothing
end


function tracktoend!(state, tracker)
    retcode = PathTracking.track!(tracker, state.samples[state.nsamples], state.R, 0.0, precondition=false)
    if retcode == :success
        state.R = 0.0
        state.pbest .= PathTracking.currx(tracker)
        state.pbest_delta = PathTracking.currresidual(tracker)
        state.status = :success
    end
    nothing
end


function checkatinfinity_norm(state, options)
    p = state.npredictions > 0 ? state.p : state.x
    infinity_norm(p) > options.maxnorm
end


"""
    checkterminate!(endgamer)

Try to check whether the endgame is converged.
"""
function checkterminate!(state, options)
    state.status == :ok || return nothing
    if state.R < options.minradius
        if state.npredictions > 1 &&
            state.pbest_delta < sqrt(options.tol)
            state.status = :success
        else
            state.status = :minradius_reached
        end
    elseif state.npredictions > 1
        Δ = unsafe_infinity_norm(state.p, state.pprev)
        if Δ < state.pbest_delta
            state.pbest = state.p
            state.pbest_delta = Δ
        end
        if Δ < options.tol
            state.R = 0.0
            state.status = :success
        end
    end
    nothing
end
