import ..ProjectiveVectors

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
    state.directions .= 0.0

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
    EndgamerResult(endgamer)
end

function play!(endgamer)
    state = endgamer.state
    while state.status == :ok
        nextsample!(endgamer)

        estimate_windingnumber!(endgamer)

        if !checkatinfinity!(endgamer)
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

    for i=1:length(x)
        push!(state.logabs_samples[i], logabs(x[i]))
    end

    nothing
end

function estimate_windingnumber!(endgamer)
    state, options = endgamer.state, endgamer.options
    nsamples = length(state.logabs_samples[1])
    # We need at least 4 samples to estimate the winding number
    nsamples < 4 && return nothing

    logh = fastlog(options.sampling_factor)
    map!(endgamer.cache.windingnumbers, state.logabs_samples) do samples
        windingnumber(samples, logh)
    end
    w = majority(endgamer.cache.windingnumbers)

    if w < 1
        state.cons_matching_estimates = 0
        state.windingnumber_estimate = 1
    elseif w == state.windingnumber_estimate
        state.cons_matching_estimates += 1
    else
        state.cons_matching_estimates = 0
        state.windingnumber_estimate = w
    end

    nothing
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

checkatinfinity!(endgamer) = checkatinfinity!(endgamer, endgamer.state.samples[1])
function checkatinfinity!(endgamer, x::ProjectiveVectors.PVector)
    cache, options, state = endgamer.cache, endgamer.options, endgamer.state

    nsamples = length(state.logabs_samples[1])

    if nsamples < 8 || state.cons_matching_estimates < 3
        return false
    end

    h = options.sampling_factor
    logh = log(h)
    w = 1

    homvar = ProjectiveVectors.homvar(x)

    range = max(1, nsamples - options.max_extrapolation_samples + 1):nsamples
    w₀ = direction(state.logabs_samples[homvar], range, h, logh, cache.direction_buffer)

    maxΔ = abs(w₀ - state.directions[homvar])
    state.directions[homvar] = w₀

    for (i, samples) in enumerate(state.logabs_samples)
        i == homvar && continue
        wᵢ = direction(samples, range, h, logh, cache.direction_buffer)
        maxΔ = max(maxΔ, abs(wᵢ - state.directions[i]))
        state.directions[i] = wᵢ
    end

    if maxΔ > 1e-2
        return false
    end

    for (i, wᵢ) in enumerate(state.directions)
        i == homvar && continue
        w = wᵢ - w₀
        if w * state.windingnumber_estimate < -0.99
            state.status = :at_infinity
            return true
        end
    end

    false
end

function direction(samples, range, h, logh, buffer)
    d = inv(1 - h)
    n = length(range)
    for (i, k) in enumerate(range.start:range.stop-1)
        buffer[i] = samples[k+1] - samples[k]
    end
    for k=2:n-1, i=1:n-k
        buffer[i] = muladd((buffer[i+1] - buffer[i]), d, buffer[i])
    end
    buffer[1] / logh
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

    return round(Int, logh / logabs((Δ1 - Δ2) / (Δ - Δ1)))
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
#     for (i, k) in enumerate(range.start:range.stop-1)
#         cons_diffs[i] = logabs_samples[k] - logabs_samples[k + 1]
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

fastlog(z) = Base.Math.JuliaLibm.log(z)

"""
    predict!(endgamer)
"""
function predict!(endgamer::Endgamer)
    state = endgamer.state

    # We do a prediction if we have at least 2 estimates
    # or at least 10 iterations
    if state.cons_matching_estimates ≥ 2 || state.iters > 10
        state.pprev .= state.p
        # we try to fit a power series
        buffer = endgamer.cache.fitpowerseries_buffer
        range = max(1, state.nsamples - 5):state.nsamples
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
