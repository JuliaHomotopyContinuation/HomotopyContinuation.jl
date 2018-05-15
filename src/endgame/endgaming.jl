import ..ProjectiveVectors

function setup!(endgamer, x, R::Float64)
    state = endgamer.state

    PathTracking.setup!(endgamer.tracker, x, R, 0.0)

    state.samples[1] .= x
    state.nsamples = 1
    for i = 1:length(state.logabs_samples)
        empty!(state.logabs_samples[i])
        push!(state.logabs_samples[i], logabs(x[i]))
    end
    state.directions .= 0.0
    state.directions_matching .= 0
    state.n_directions_matching = 0

    state.npredictions = 0
    state.iters = 0
    state.R = R
    state.status = :ok

    state.windingnumber = 1
    state.cons_matching_estimates = 0

    nothing
end


function play(endgamer, x, R)
    setup!(endgamer, x, R)
    play!(endgamer)
    EndgamerResult(endgamer)
end

function play!(endgamer)
    state, options, cache = endgamer.state, endgamer.options, endgamer.cache
    while state.status == :ok
        nextsample!(endgamer.tracker, state, options)

        if state.nsamples > 3 # directions! and windingnumber! need at least 4 samples
            directions!(state, options, cache)
            estimate_windingnumber!(state, options, cache)

            predict!(state, options, cache)
            checkatinfinity!(state, options, cache)
        end

        checkterminate!(state, endgamer.options)
    end

    # cleanup
    endgamer
end

nsamples(state) = length(state.logabs_samples[1])

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
function nextsample!(tracker::PathTracking.PathTracker, state, options)
    state.iters += 1

    R, λ = state.R, options.sampling_factor
    λR = λ * R

    retcode = PathTracking.track!(tracker, state.samples[state.nsamples], R, λR,
        precondition=false, checkstartvalue=false, emit_update=false)
    if retcode != :success
        state.status = :tracker_failed
        return nothing
    end

    state.R = λR
    # now we have to add the new sample point to samples
    update_samples!(state, options, PathTracking.currx(tracker), λR)

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

    for i=1:length(x)
        push!(state.logabs_samples[i], logabs(x[i]))
    end

    nothing
end


"""
    directions!(state, options, cache)

Estimate the lowest order exponent of the fractional power series representation
of x(t).
"""
function directions!(state, options, cache)
    h = options.sampling_factor
    logh = fastlog(h)
    range = samples_range(state, options)
    buffer = cache.direction_buffer
    stable = true
    for (i, samplesᵢ) in enumerate(state.logabs_samples)
        dᵢ = extrapolate_direction(samplesᵢ, range, h, logh, buffer)
        if abs(state.directions[i] - dᵢ) < 0.0005
            state.directions_matching[i] += 1
            state.directions[i] = dᵢ
        else
            state.directions_matching[i] = 0
            state.directions[i] = dᵢ
        end
    end

    if stable
        state.n_directions_matching += 1
    else
        state.n_directions_matching = 0
    end
    # @show state.n_directions_matching
    state.n_directions_matching
end

"""
    extrapolate_direction(samples, range, h, logh, buffer)

Compute the directions by using Richardson extrapolation as explained
in Section 4 of [^1].

[^1]: Huber, Birkett, and Jan Verschelde. "Polyhedral end games for polynomial continuation." Numerical Algorithms 18.1 (1998): 91-108.
"""
function extrapolate_direction(samples, range, h, logh, buffer)
    d = inv(1 - h)
    n = length(range)
    for i=1:n-1
        buffer[i] = samples[range[i+1]] - samples[range[i]]
    end

    @inbounds begin
        for k=2:n-1, i=1:n-k
            buffer[i] = muladd((buffer[i+1] - buffer[i]), d, buffer[i])
        end
    end
    buffer[1] / logh
end

"""
    estimate_windingnumber!(state, options, cache)

This estimates the winding number of the path ``x(t)`` at ``x(0)``.
For this we look at consecutive differences of samples obtained
from the geometric series. The windingnumber is estimated coordinate wise
and since it is the same for each coordinate we take the most common
result.
"""
function estimate_windingnumber!(state, options, cache)
    logh = fastlog(options.sampling_factor)
    map!(cache.windingnumbers, state.logabs_samples) do samples
        windingnumber(samples, logh)
    end
    w = majority(cache.windingnumbers)

    # w < 1 is just wrong
    if w < 1
        state.cons_matching_estimates = 0
        state.windingnumber = 1
    elseif w == state.windingnumber
        state.cons_matching_estimates += 1
    else
        state.cons_matching_estimates = 0
        state.windingnumber = w
    end

    nothing
end

function samples_range(state, options)
    samples_range(length(state.logabs_samples[1]), options)
end
function samples_range(nsamples::Int, options)
    max(1, nsamples - options.max_extrapolation_samples + 1):nsamples
end

"""
    windingnumber(logabs_samples, logh)

Estimate the current winding number from the samples.
`logabs_samples` represent the sample points obtained from the geometric
path samples ``log|x(h^k₀)|``.
We compute an error expansion of the fractional power series of the path
to approximate the winding number. See [^1] for details.

[^1]: Huber, Birkett, and Jan Verschelde. "Polyhedral end games for polynomial continuation." Numerical Algorithms 18.1 (1998): 91-108.
"""
function windingnumber(logabs_samples, logh)
    nsamples = length(logabs_samples)
    Δ = logabs_samples[nsamples-3] - logabs_samples[nsamples-2]
    Δ1 = logabs_samples[nsamples-2] - logabs_samples[nsamples-1]
    Δ2 = logabs_samples[nsamples-1] - logabs_samples[nsamples]

    w = logh / logabs((Δ1 - Δ2) / (Δ - Δ1))
    if !isfinite(w) # we have Δ == Δ1
        return 0
    end
    return round(Int, w)
end


"""
    majority(values)

Returns the value which occurs most often in `values`. If there is a tie
the first one is chosen.
"""
function majority(vec::Vector{Int})
    maxval = 0
    maxnvals = 0
    n = length(vec)
    for i=1:n-1
        nvals = 1
        v = vec[i]
        if v == maxval
            continue
        end
        for j=i+1:n
            if vec[j] == v
                nvals += 1
            end
        end
        if nvals > maxnvals
            maxnvals = nvals
            maxval = v
        end
    end
    maxval
end



"""
    checkatinfinity!(state, options, cache)

We compute an error expansion of the fractional power series of the path
to approximate the "direction" of the path. See [^1] for details.

[^1]: Huber, Birkett, and Jan Verschelde. "Polyhedral end games for polynomial continuation." Numerical Algorithms 18.1 (1998): 91-108.
"""
function checkatinfinity!(state, options, cache)
    i₀ = ProjectiveVectors.homvar(state.samples[1])
    if state.n_directions_matching < 1
        return false
    end
    # return false
    m = state.windingnumber
    r₀ = state.directions[i₀]

    if state.directions_matching[i₀] > 4 &&  r₀ > 0.02 # we assume a winding number of at most 50
        for (i, rᵢ) in enumerate(state.directions)
            i == i₀ && continue
            # rᵢ should be really near zero
            if rᵢ < 0.001 && state.directions_matching[i] > 4
                state.status = :at_infinity
                return true
            end
        end
    else
        p = state.npredictions > 0 ? state.p : state.samples[state.nsamples]
        if infinity_norm(p) > options.maxnorm
            state.status = :at_infinity
            return true
        end
    end

    false
end

#
# It turns out that this extrapolation scheme doesn't seem to work better
# than the naive way.
#
#
# """
#     extrapolate_winding_number(logabs_samples, range, h, buffer)
#
# Extrapolate an estimate of the winding number from `logabs_samples[range]`.
# `logabs_samples` represent the sample points obtained from the geometric
# path samples ``log|x(h^k₀)|``.
# The extrapolation scheme can be interpreted as a special form of Richardson extrapoliation
# and is described in [^1].
# The buffer is needed to avoid allocations and it needs to be at least of length
# `length(range)-1`.
#
# [^1]: Huber, Birkett, and Jan Verschelde. "Polyhedral end games for polynomial continuation." Numerical Algorithms 18.1 (1998): 91-108.
# """
# function extrapolate_winding_number(logabs_samples, range, h, buffer)
#     # put consecutive differences into buffer
#     cons_diffs = buffer # rename for clarity
#     n = length(range)
#     for i=1:n-1
#         cons_diffs[i] = logabs_samples[range[i]] - logabs_samples[range[i+1]]
#     end
#     # no we need to compute the error expansion
#     log_err_exp = cons_diffs
#     for i = 1:length(range)-2
#         log_err_exp[i] = logabs(cons_diffs[i] - cons_diffs[i + 1])
#     end
#
#     logh = fastlog(h)
#     extrapolation = log_err_exp # we rename the buffer again
#
#     n = length(range) - 3
#     for i = 1:n
#         extrapolation[i] = log_err_exp[i+1] - log_err_exp[i]
#     end
#     w = logh / extrapolation[n]
#     for k = 1:n-1
#         d = inv(h^(k/w) - 1.0)
#         for i=1:n-k
#             extrapolation[i] = muladd((extrapolation[i] - extrapolation[i+1]), d, extrapolation[i])
#         end
#         w = logh / extrapolation[n - k]
#     end
#     w
# end



"""
    predict!(endgamer)
"""
function predict!(state, options, cache)
    state.pprev .= state.p
    # we try to fit a power series
    buffer = cache.fitpowerseries_buffer
    range = samples_range(state.nsamples, options)
    fitpowerseries!(state.p, state.samples, range, options.sampling_factor,
        state.windingnumber, buffer)
    state.npredictions += 1

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

    d = inv(r - 1)
    for i = 1:n-1
        @. buffer[i] = (r * samples[range[i+1]] - samples[range[i]]) * d
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
        end
    end
    return nothing
end
