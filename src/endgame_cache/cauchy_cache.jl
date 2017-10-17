mutable struct CauchyEndgameCache{xLowType} <: AbstractEndgameCache
    samples::Vector{xLowType}
end


function alg_cache(alg::CauchyEndgame, tracker)
    xLowType = typeof(tracker.low.x)
    samples = Vector{xLowType}(alg.samples_per_loop)
    CauchyEndgameCache(samples)
end


function predict!(endgamer, cache::CauchyEndgameCache)
    @unpack iter, status, xs, R, predictions, alg, options = endgamer
    λ = options.geometric_series_factor

    # we need at least 4 values for the first heuristic
    if length(xs) < 4
        return nothing
    end

    if status == NotStarted && !firstheuristic(R, λ, xs[end - 3], xs[end - 2], xs[end - 1], xs[end], alg.L)
        return nothing
    end

    # now we need to collect sample points, stored into the cache
    retcode = loop!(endgamer, xs[end], R, cache)
    if retcode == :success
        prediction = predict_with_cif(cache.samples)
        push!(predictions, prediction)
    elseif retcode == :tracker_failed || retcode == :winding_number_too_high
        # TODO: Can we do something here?
        endgamer.status = Failed
        endgamer.failurecode = retcode
    elseif retcode == :heuristic_failed
        # everthing okay
    else
        endgamer.status = Failed
        endgamer.failurecode = retcode
        warn("Unhandled retcode $(retcode) in `predict!`")
    end
end

"""
    firstheuristic(t, x_R, x_λR, x_λ2R, x_λ3R, L = 0.75)

This is the first heuristic when to start the endgame.
It checks whether two consectuive appromixations of the first nonvanishing
monomial of the Puiseux series expansion of x(t) at 0 agree.
For more details see "A Parallel Endgame " by Bates, Hauenstein and Sommese [^1],
page 8 and 9.

[^1]: Bates, Daniel J., Jonathan D. Hauenstein, and Andrew J. Sommese. "A Parallel Endgame." Contemp. Math 556 (2011): 25-35.
"""
function firstheuristic(t, λ, x_R, x_λR, x_λ2R, x_λ3R, L = 0.75)
    # failsafe mechanism
    if t < 1e-8
        return true
    end
    v = rand(Complex128, length(x_R))
    g(x_R, x_λR, x_λ2R) = log(abs2((v⋅x_λR - v⋅x_λ2R) / (v⋅x_R - v⋅x_λR))) / log(λ)

    g_R = g(x_R, x_λR, x_λ2R)
    g_λR = g(x_λR, x_λ2R, x_λ3R)

    if g_R <= 0
        return false
    end

    L < g_R / g_λR < 1 / L
end

"""
    loop(endgame, x, radius, cache::CauchyCache)

Tracks the implicit defined path z(t) around the `n`-gon with vertices
``r⋅exp(i2πk/n)`` where `n=samples_per_loop`.
"""
function loop!(endgamer, x, radius::Real, cache::CauchyEndgameCache)
    @unpack tracker, alg, windingnumber, status, options = endgamer
    @unpack max_winding_number, abstol= options
    @unpack samples_per_loop, loopclosed_tolerance = alg
    @unpack samples = cache

    if length(samples) != windingnumber * samples_per_loop
        resize!(samples, windingnumber * samples_per_loop)
    end

    started = status == Started

    samples[1] = x
    unitroots = UnitRootsIterator(radius, samples_per_loop)

    start = first(unitroots)
    c = 1
    for (k, finish) in enumerate(Iterators.drop(unitroots, 1))
        run!(tracker, samples[k], start, finish)
        retcode, sol = solution(tracker)

        if retcode != :success
            @show tracker.t
            @show retcode
            return :tracker_failed
        end

        # apply the second heuristic to tell us whether it makes sense to start the endgame
        if !started && k == samples_per_loop && windingnumber == 1 &&
            !secondheuristic(radius, samples, abstol, alg.K)
            return :heuristic_failed
        end

        # check whether the loop is closed. This can only happen if we finished one circle.
        if k % samples_per_loop == 0
            if projectivenorm2(x, sol) < loopclosed_tolerance
                break
            else
                c += 1
            end
        end

        if c > max_winding_number
            return :winding_number_too_high
        end

        if c ≤ windingnumber
            samples[k + 1] = sol
        else
            push!(samples, sol)
        end
        start = finish
    end

    endgamer.windingnumber = c
    resize!(samples, c * samples_per_loop)

    :success
end

"""
    secondheuristic(t, samplepoints, β, K=0.5)

This is the second heuristic when to start the endgame.
It enforces that the values around the lopp do not differ radically.
This functions assumes that `samplepoints` has been collected on ``x(Re^{iθ})`` for ``θ∈[0,2π]``.
See "A Parallel Endgame " by Bates, Hauenstein and Sommese [^1], page 9 for more details.

[^1]: Bates, Daniel J., Jonathan D. Hauenstein, and Andrew J. Sommese. "A Parallel Endgame." Contemp. Math 556 (2011): 25-35.
"""
function secondheuristic(radius, samplepoints, β=1e-8, K=0.5)
    # failsafe
    if radius < 1e-14
        true
    end
    # winding number is 1 so this is safe
    norms = norm.(samplepoints)
    m = minimum(norms)
    M = maximum(norms)
    (M - m < β) || (m / M > K)
end

"""
    predict_with_cif(samplepoints)

Predict the value of a function based on the samplepoints using the cauchy integral formula.
This function assumes that the samplepoints are sampled uniformly from a circle around 0,
since then the integral (approximated with the trapezoidal rule) is simply the mean of
all samplepoints.
"""
predict_with_cif(samplepoints) = mean(samplepoints)
