function setup!(endgame, x, R)
    setup!(endgame.tracker, x, R, 0.0im)
    reset!(endgame.state, x, R)

    # The endgame fails if we already have a solution. Hence we check it now.
    r = residual(endgame.tracker, x, 0.0)
    if abs(r) < endgame.options.tol
        endgame.state.pbest .= x
        endgame.state.pbest_delta = r
        endgame.state.status = :success
        endgame.state.R = 0.0im
    end
    nothing
end

"""
    runendgame(endgame, x, t)

Start the given endgame to run from `x` at `t` towards `t=0`.

## Math / Algorithm

The Endgame is a hybrid of the Cauchy endgame and the power series endgame.

### Power series representation
Each path ``x(t)`` has an asymptotic expansion as a fractional power series
```math
xᵢ(t) = a t^{\\frac{wᵢ}{m}} (1 + ∑_{j≥1} aᵢⱼt^{\\frac{j}{m}})
```
We want to approximate the ratio ``wᵢ/m``. If ``wᵢ/m > 0`` we have that ``xᵢ(t)`` goes to
``0`` and if ``wᵢ/m < 0`` we have that ``xᵢ(t)`` goes to ``∞``.
Since we have `t^{(k)}=h^kR₀` in the endgame we can approximate it by computing
```julia
\\log(x(t^{(k)})) - \\log(x(t^{(k-1)})) ≈ wᵢ/m \\log(h) + \\mathcal{O}(h^{1/m})
```
See [^1] for some more details.


### Cauchy

The main idea is to use [Cauchy's integral formula](https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula)
to predict the solution of the path ``x(t)``, i.e. ``x(0)``.
At each iteration we are at some point ``(x, t)``. We then track the polygon defined
by ``te^{i2πk/n}`` until we end again at ``x``. Here ``n`` is the number of samples we take
per loop. See [^2] for some more details.


[^1]: Huber, Birkett, and Jan Verschelde. "Polyhedral end games for polynomial continuation." Numerical Algorithms 18.1 (1998): 91-108.
[^2]: Bates, Daniel J., Jonathan D. Hauenstein, and Andrew J. Sommese. "A Parallel Endgame." Contemp. Math 556 (2011): 25-35.
"""
function runendgame(endgame, x, R)
    setup!(endgame, x, R)
    runendgame!(endgame)
    EndgameResult(endgame)
end

function runendgame!(endgame)
    state, tracker, options, cache = endgame.state, endgame.tracker, endgame.options, endgame.cache
    while state.status == :ok
        nextsample!(tracker, state, options)
        update_derived_from_samples!(state, options, cache)
        if in_endgame_zone!(state)
            predict_infinity_check!(state, tracker, options, cache)
        end
        if checkatinfinity_norm(state, options)
            state.status = :at_infinity
        end
        checkterminate!(state, endgame.options)
    end

    endgame
end


"""
    nextsample!(endgame)

Obtain a new sample on the geometric series and add it to `endgame`.
"""
function nextsample!(tracker::PathTracker, state, options)
    state.iters += 1

    R, λ = state.R, options.sampling_factor
    λR = λ * R

    retcode = track!(tracker, state.x, R, λR, setup_patch=false, compute_ẋ=false, checkstartvalue=false)
    if retcode != PathTrackerStatus.success
        state.status = :tracker_failed
        return nothing
    end

    state.R = λR
    state.x .= currx(tracker)
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
    n, dirs = state.nsamples, state.directions

    state.in_endgame_zone && return true

    if state.R < 1e-8
        state.in_endgame_zone = true
        return true
    end

    n < 5 && return false

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

analyze_finite_atinfinity(s, o, c) = analyze_finite_atinfinity(s.x, s, o, c)
function analyze_finite_atinfinity(x::ProjectiveVectors.PVector, state, options, cache)
    options.check_at_infinity || :projective
    checkfinite(state) && return :finite
    checkatinfinity(state, options) && return :at_infinity
    :unknown
end

checkfinite(state) = checkfinite(state, state.x)
function checkfinite(state, x::ProjectiveVectors.PVector{<:Complex})
    i₀, = ProjectiveVectors.homvars(x)
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
function checkatinfinity(state, options, x::ProjectiveVectors.PVector{<:Complex})
    LOOKBACK = 1
    TOL = 1e-3
    MIN_AT_INFINITY_RATIO = 0.032

    i₀, = ProjectiveVectors.homvars(x)
    n = state.nsamples
    dirs = state.directions

    # We want the homogenization variable to be at least somewhat small
    if state.R > 1e-8 && ProjectiveVectors.norm_affine_chart(x) < options.minimal_maxnorm
        return false
    end
    for i=1:length(x)
        i == i₀ && continue
        di_n = dirs[i, n] - dirs[i₀, n]
        di_n > -MIN_AT_INFINITY_RATIO && continue
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

    path_type = analyze_finite_atinfinity(state, options, cache)
    if path_type == :finite || path_type == :projective
        # we need to determine the windingnumber
        if state.windingnumber == 0
            determine_windingnumber!(state, tracker, options, cache)
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

function checkatinfinity_norm(state, options)
    p = state.npredictions > 0 ? state.p : state.x
    options.check_at_infinity && infinity_norm(p) > options.maxnorm
end


"""
    checkterminate!(endgame)

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
