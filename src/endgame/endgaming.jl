function setup!(endgamer, x, R::Float64)
    state = endgamer.state

    # We want to fix the patch, but before we need to make sure that we are on the correct
    # local patch. Hence disable fixpatch, do the normal setup, and enable the patch fixing again
    PathTracking.fixpatch!(endgamer.tracker, false)
    PathTracking.setup!(endgamer.tracker, x, R, 0.0)
    PathTracking.fixpatch!(endgamer.tracker, true)

    state.samples[1] .= x
    state.nsamples = 1
    for i = 1:length(state.logabs_samples)
        empty!(state.logabs_samples[i])
        push!(state.logabs_samples[i], logabs(x[i]))
    end

    state.npredictions = 0
    state.iters = 0
    state.R = R
    state.status = :ok

    state.windingnumber_estimate = 1
    state.cons_matching_estimates = 0

    nothing
end


function play(endgamer, x, R)
    setup!(endgamer, x, R)
    play!(endgamer)
end

function play!(endgamer)
    state = endgamer.state
    while state.status == :ok
        nextsample!(endgamer)

        if !confident_in_windingnumber(state, endgamer.options)
            estimate_windingnumber!(endgamer)
        end
        if confident_in_windingnumber(state, endgamer.options) &&
            state.windingnumber_estimate == 1

            try_to_jump_to_target!(endgamer)
        else
            predict!(endgamer)
        end

        checkterminate!(state, endgamer.options)
    end

    # cleanup
    PathTracking.fixpatch!(endgamer.tracker, false)
    endgamer
end

function confident_in_windingnumber(state::State, options::Options)
    if state.cons_matching_estimates < 5
        return false
    end
    true
end

"""
    nextsample!(endgamer)

Obtain a new sample on the geometric series and add it to `endgamer`.
"""
function nextsample!(endgamer::Endgamer)
    state = endgamer.state
    state.iters += 1

    state.status == :ok || return nothing

    R, λ = state.R, endgamer.options.sampling_factor
    λR = λ * R

    retcode = PathTracking.track!(endgamer.tracker, state.samples[state.nsamples], R, λR)
    if retcode != :success
        warn(retcode)
        state.status = :tracker_failed
        return nothing
    end

    state.R = λR
    # now we have to add the new sample point to samples
    update_samples!(state, endgamer.options, PathTracking.currx(endgamer.tracker), λR)

    nothing
end

function update_samples!(state, options, x, λR)
    if state.nsamples == length(state.samples)
        push!(state.samples, copy(x))
        state.nsamples += 1
    else
        state.nsamples += 1
        state.samples[state.nsamples] .= x
    end

    # We only need to add new logabs samples if we are still unsure about the
    # winding number
    if !confident_in_windingnumber(state, options)
        for i=1:length(x)
            push!(state.logabs_samples[i], logabs(x[i]))
        end
    end
    nothing
end

function estimate_windingnumber!(endgamer)
    state, options = endgamer.state, endgamer.options
    nsamples = length(state.logabs_samples[1])
    # We need at least 4 samples to estimate the winding number
    nsamples < 4 && return nothing
    # We use the last samples for the prediction
    range = max(1, nsamples - options.max_extrapolation_samples - 1):nsamples
    extrapolated_w = extrapolate_winding_number(
        rand(state.logabs_samples),
        range,
        options.sampling_factor,
        endgamer.cache.windingnumber_extrapolation_buffer)

    # we round since the winding number is an integer
    w = round(Int, extrapolated_w)
    if w == state.windingnumber_estimate
        state.cons_matching_estimates += 1
    else
        state.cons_matching_estimates = 0
        state.windingnumber_estimate = w
    end

    nothing
end


"""
    extrapolate_winding_number(logabs_samples, range, h, buffer)

Extrapolate an estimate of the winding number from `logabs_samples[range]`.
`logabs_samples` represent the sample points obtained from the geometric
path samples ``log|x(h^k₀)|``.
The extrapolation scheme can be interpreted as a special form of Richardson extrapoliation
and is described in [^1].
The buffer is needed to avoid allocations and it needs to be at least of length
`length(range)-1`.

[^1]: Huber, Birkett, and Jan Verschelde. "Polyhedral end games for polynomial continuation." Numerical Algorithms 18.1 (1998): 91-108.
"""
function extrapolate_winding_number(logabs_samples, range, h, buffer)
    # put consecutive differences into buffer
    cons_diffs = buffer # rename for clarity
    for (i, k) in enumerate(range.start:range.stop-1)
        cons_diffs[i] = logabs_samples[k] - logabs_samples[k + 1]
    end
    # no we need to compute the error expansion
    log_err_exp = cons_diffs
    for i = 1:length(range)-2
        log_err_exp[i] = logabs(cons_diffs[i] - cons_diffs[i + 1])
    end

    logh = fastlog(h)
    extrapolation = log_err_exp # we rename the buffer again

    n = length(range) - 3
    for i = 1:n
        extrapolation[i] = log_err_exp[i+1] - log_err_exp[i]
    end
    w = logh / extrapolation[n]
    for k = 1:n-1
        d = inv(h^(k/w) - 1.0)
        for i=1:n-k
            extrapolation[i] = muladd((extrapolation[i] - extrapolation[i+1]), d, extrapolation[i])
        end
        w = logh / extrapolation[n - k]
    end
    w
end

fastlog(z) = Base.Math.JuliaLibm.log(z)

"""
    predict!(endgamer)
"""
function predict!(endgamer::Endgamer)
    state = endgamer.state

    # We do a prediction if we have at least 2 estimates
    if state.cons_matching_estimates ≥ 2
        state.pprev .= state.p

        # we try to fit a power series
        buffer = endgamer.cache.fitpowerseries_buffer
        range = max(1, state.nsamples - state.cons_matching_estimates + 1):state.nsamples
        resize_buffer!(buffer, length(range))
        fitpowerseries!(state.p, state.samples, range, endgamer.options.sampling_factor,
            state.windingnumber_estimate, buffer)

        state.npredictions += 1
    end
    nothing
end

function resize_buffer!(buffer, target_length)
    n = length(buffer)
    if target_length ≤ n
        return nothing
    end

    for _=1:target_length-n
        push!(buffer, copy(buffer[1]))
    end
    nothing
end


"""
    fitpowerseries!(p, samples, range, h, w, buffer)

Fit a fractional power series to the samples defined in `samples[range]`. The samples are assumed
to be obtained by the geometric series ``s(a h^k)``. `w` is the winding number
of the fractional power series. For details to the implementation take a look
at section 4 and 5 in [^1].

[^1]: Schoenberg, I. J. "On polynomial interpolation at the points of a geometric progression." Proceedings of the Royal Society of Edinburgh Section A: Mathematics 90.3-4 (1981): 195-207.
"""
function fitpowerseries!(p, samples, range, h, w, buffer)
    r = (1 / h)^(1/w)
    n = length(range)
    @assert length(buffer) ≥ n - 1

    d = inv(r - 1)
    @inbounds for (i, k) in enumerate(range.start:range.stop-1)
        @. buffer[i] = (r * samples[k+1] - samples[k]) * d
    end
    r_pow = r
    for k=2:n-1
        r_pow *= r
        d = inv(r_pow - 1)
        @inbounds for i=1:n-k
            @. buffer[i] = (r_pow * buffer[i+1] - buffer[i]) * d
        end
    end
    p .= buffer[1]
end


"""
    checkterminate!(endgamer)

Try to check whether the endgame is converged.
"""
function checkterminate!(state, options)
    state.status == :ok || return nothing
    if state.R < options.minradius
        state.status = :minradius_reached
        return nothing
    end

    if state.npredictions > 1
        Δ = unsafe_infinity_norm(state.p, state.pprev)
        if Δ < options.tol
            state.R = 0.0
            state.status = :success
        elseif infinity_norm(state.p) > options.maxnorm
            state.status = :at_infinity
        end
    end
    return nothing
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

    retcode = PathTracking.track!(endgamer.tracker, state.samples[end], state.R, 0.0)
    if retcode == :success
        state.p .= PathTracking.currx(endgamer.tracker)
        state.R = 0.0
        state.status = :success
        state.npredictions += 1
    end
    nothing
end
