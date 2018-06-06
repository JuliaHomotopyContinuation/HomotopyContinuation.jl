"""
    update_directions!(state, options, cache)

### Math
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


This function approximates ``wᵢ/m`` using Richardson extrapolation and stores them in
`state.directions`.


[^1]: Huber, Birkett, and Jan Verschelde. "Polyhedral end games for polynomial continuation." Numerical Algorithms 18.1 (1998): 91-108.
"""
function update_directions!(state, options, cache)
    h, logh = options.sampling_factor, cache.log_sampling_factor
    S = state.logabs_samples
    B = cache.direction_buffer
    N = state.nsamples
    m = size(S, 1)

    range = max(1, N - options.max_extrapolation_samples + 1):N
    n = length(range)
    for j=1:n-1, i=1:m
        B[i, j] = S[i, range[j+1]] - S[i, range[j]]
    end

    # This is just a richardson extrapolation
    c_k = 1
    for k=2:n-1
        c_k *= h
        for j=1:n-k, i=1:m
            B[i, j] = (B[i, j+1] - c_k * B[i, j]) / (1 - c_k)
        end
    end

    for i=1:m
        state.directions[i, N] = B[i, 1] / logh
    end
    nothing
end


#
# """
#     predict_powerseries!(state, options, cache)
#
# Fit a fractional power series to the samples defined in `samples[range]`. The samples are assumed
# to be obtained by the geometric series ``s(a h^k)``. `w` is the winding number
# of the fractional power series. For details to the implementation take a look
# at section 4 and 5 in [^1].
#
# [^1]: Schoenberg, I. J. "On polynomial interpolation at the points of a geometric progression." Proceedings of the Royal Society of Edinburgh Section A: Mathematics 90.3-4 (1981): 195-207.
# """
# function predict_powerseries!(state, options, cache)
#     state.pprev .= state.p
#     # we try to fit a power series
#     buffer = cache.fitpowerseries_buffer
#     range = samples_range(state.nsamples, options)
#     state.npredictions += 1
#     fitpowerseries!(state.p, state.samples, range, options.sampling_factor,
#         max(1, state.windingnumber), buffer)
#     nothing
# end
# function fitpowerseries!(p, samples, range, h, w, buffer)
#     r = w == 1 ? inv(h) : (1 / h)^(1/w)
#     n = length(range)
#
#     d = inv(r - 1)
#     for i = 1:n-1
#         @. buffer[i] = (r * samples[range[i+1]] - samples[range[i]]) * d
#     end
#     r_pow = r
#     for k=2:n-1
#         r_pow *= r
#         d = inv(r_pow - 1)
#         @inbounds for i=1:n-k
#             @. buffer[i] = (r_pow * buffer[i+1] - buffer[i]) * d
#         end
#     end
#     p .= buffer[1]
# end
