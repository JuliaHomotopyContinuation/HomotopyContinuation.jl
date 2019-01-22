"""
    setup_unitroots!(cache, samples_per_loop)

Fill `cache.unitroots` with the unit roots of order up to
`samples_per_loop - 1`.
"""
function setup_unitroots!(cache, samples_per_loop)
    unitroots = cache.unitroots
    if length(unitroots) ≠ samples_per_loop
        empty!(unitroots)
        push!(unitroots, 1.0+0.0im)
        for k=2:samples_per_loop
            push!(unitroots, cis(2π*(k-1)/samples_per_loop))
        end
    end
    unitroots
end


"""
    determine_windingnumber(tracker, state, options, cache)

Try to determine the windingnumber by tracing around the origin.
Returns a tuple `(retcode, windingnumber, p)` where `windingnumber`
is 0 if the procedure failed, i.e., `retcode != :success`. `p` is a
prediction of `x(0)`.
"""
function determine_windingnumber!(state, tracker, options, cache)
    samples_per_loop = 4
    tol = options.cauchy_loop_closed_tolerance
    unitroots = setup_unitroots!(cache, samples_per_loop)
    x = state.x
    xk1 = cache.xbuffer .= x
    p = cache.pbuffer .= x
    R = state.R
    θk = θk1 = complex(R, 0.0)

    k = c = 1
    # we try to approximate the maximum size of the loop to normalize our error
    dmax = 0.0
    while c ≤ options.maxwindingnumber
        k += 1
        θk = R * unitroots[(k - 1) % samples_per_loop + 1]
        # We go around the unit circle in an `n`-gon
        retcode = track!(tracker, xk1, θk1, θk, setup_patch=false)
        if retcode != PathTrackerStatus.success
            return :loop_failed_tracking_failed
        end
        xk = currx(tracker)
        Δ  = unsafe_infinity_norm(x, xk)
        dmax = max(dmax, Δ)
        if (k - 1) % samples_per_loop == 0
            # @show c, Δ/dmax
            if Δ/dmax < tol
                state.pprev .= state.p
                state.p .= p ./ (c * samples_per_loop)
                state.npredictions += 1
                state.windingnumber = c
                return :success
            end
            c += 1
        end
        p .+= xk
        xk1 .= xk
        θk1 = θk
    end

    return :max_winding_reached
end

"""
    predict_cif!(state, tracker, options, cache)

Try to predict the value of `x(0)` using [Cauchy's integral formula](https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula)
At each iteration we are at some point ``(x, t)``. We then track the polygon defined
by ``te^{i2πk/n}`` until we end again at ``x``. Here ``n`` is the number of samples we take
per loop.
This assumes that the correct windingnumber was determined previously.
Returns a symbol indicating whether the prdiction was successfull.
"""
function predict_cif!(state, tracker, options, cache)
    state.windingnumber ≠ 0 || return :windingnumber_0

    samples_per_loop = options.cauchy_samples_per_loop
    R = state.R
    θk = θk1 = complex(R, 0.0)
    unitroots = setup_unitroots!(cache, samples_per_loop)
    xk1 = cache.xbuffer .= state.x
    p = cache.pbuffer .= xk1

    for k = 2:(state.windingnumber * samples_per_loop)
        θk = R * unitroots[(k - 1) % samples_per_loop + 1]
        # We go around the unit circle in an `n`-gon
        retcode = track!(tracker, xk1, θk1, θk, setup_patch=false)
        if retcode != PathTrackerStatus.success
            return :tracking_failed
        end
        xk = currx(tracker)
        p .+= currx(tracker)
        xk1 .= xk
        θk1 = θk
    end
    state.pprev .= state.p
    state.p .= p ./ (state.windingnumber * samples_per_loop)
    state.npredictions += 1
    :success
end
