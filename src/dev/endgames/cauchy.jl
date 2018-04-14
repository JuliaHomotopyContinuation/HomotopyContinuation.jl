"""
    Cauchy(;kwargs...)

The main idea of the Cauchy Endgame is to use [Cauchy's integral formula](https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula)
to predict the solution of the path ``x(t)``, i.e. ``x(0)``.
At each iteration we are at some point ``(x, t)``. We then track the polygon defined
by ``te^{i2πk/n}`` until we end again at ``x``. Here ``n`` is the number of samples we take
per loop.

The following options are available:
* `samples_per_loop=8`: The number of samples we take at one loop.
* `L=0.75` and ``K=0.5``: These are paramters for heuristics. For more details see "A Parallel Endgame " by Bates, Hauenstein and Sommese [^1],
    page 8 and 9.

[^1]: Bates, Daniel J., Jonathan D. Hauenstein, and Andrew J. Sommese. "A Parallel Endgame." Contemp. Math 556 (2011): 25-35.
"""
struct Cauchy <: AbstractEndgame
    samples_per_loop::Int
    cycle_time_cutoff::Float64
    # parameter for the cycle test heuristic
    L::Float64
    num_needed_for_stabilization::Int
    # parameter for the ratio heuristic
    #β::Float64
    K::Float64
end

function Cauchy(;
    samples_per_loop=8,
    cycle_time_cutoff=1e-8,
    L=0.75,
    num_needed_for_stabilization=3,
    K=0.5)
    Cauchy(
        samples_per_loop,
        cycle_time_cutoff,
        L,
        num_needed_for_stabilization,
        #β,
        K)
end

"""
    CauchyCache

Structure to avoid allocations in the cauchy endgame.

## Fields
* `samples`: A buffer holding prior sample points. The length is initialized
with `samples_per_loop + 1`.
* `unitroots::Vector{Complex{Float64}}` Vector of length `samples_per_loop`
holding the values ``exp(i2π k/N)`` for ``k=0,..,N-1`` where ``N`` is `samples_per_loop`.
"""
struct CauchyCache{T}
    samples::Vector{Vector{T}}
    unitroots::Vector{Complex{Float64}}
end

function cache(alg::Cauchy, x)
    N = alg.samples_per_loop
    samples = [copy(x) for _=0:N]
    unitroots = [cis(2π/N*k) for k=0:N-1]
    CauchyCache(samples, unitroots)
end


"""
    CauchyState

"""
struct CauchyState
    k_over_c_values::Vector{Float64}
end

state(alg::Cauchy, x) = CauchyState(Vector{Float64}())
reset!(state::CauchyState, ::Cauchy, x) = empty!(state.k_over_c_values)

function predict!(prediction, alg::Cauchy, state::CauchyState, cache::CauchyCache, tracker, xs, R, npredictions, options::EndgameOptions)
    λ = options.λ
    check_heuristic = npredictions == 0

    if check_heuristic && !cycleheuristic!(state, alg, xs, R, λ)
        return :heuristic_failed, 0
    end
    retcode, windingnumber =
        loop!(cache, alg, tracker, xs[1], R,
              alg.samples_per_loop, options.tol, check_heuristic)

    if retcode != :success
        return retcode, 0
    end

    nsamples = windingnumber * alg.samples_per_loop
    predict_with_cif!(prediction, cache, nsamples)

    return :success, windingnumber
 end

"""
    loop!(cache, alg, tracker, x, R, n, tol, check_heuristic)

Track `x` along a regular `n`-gon around the origin with radius `R` using `tracker`.
The loop is considered closed if the distance is less than `tol`.
Assumption is ``H(x, R)≈0``.
"""
function loop!(cache::CauchyCache, alg::Cauchy, tracker, x, R, samples_per_loop, tol, check_heuristic=true)
    samples_buffer = cache.samples
    unitroots = cache.unitroots

    # c is the winding number
    c = 1
    k = 1

    # as the start of our loop we can use our current point
    Θk₋₁ = R * unitroots[k]
    Θk = Θk₋₁
    xk₋₁ = samples_buffer[k] .= x

    # we need to bring the values on a fixed affine patch
    # otherwise this all *does not work* properly.
    # In order to avoid to deal with at_infinity choose a
    # simple patch which should be good
    maxind, _ = findmax(abs2, xk₋₁) # defined in Utilities
    scale!(xk₋₁, inv(xk₋₁[maxind]))

    while c ≤ 50 # fallback if everything fails...
        k += 1
        # We go around the unit circle in an `n`-gon
        Θk = R * unitroots[(k - 1) % samples_per_loop + 1]
        xk = samples_buffer[k]
        retcode = PathTracking.track!(xk, tracker, xk₋₁, Θk₋₁, Θk)
        tracker.state.iter
        if retcode != :success
            @show retcode
            return (:loop_failed_tracking_failed, 0)
        end
        scale!(xk, inv(xk[maxind]))

        # to avoid unnecessary work we can use a nice heuristic after one loop
        if check_heuristic &&
           k == samples_per_loop &&
           !ratioheuristic(R, samples_buffer, samples_per_loop, tol, alg.K)

           return :(:heuristic_failed, 0)
        end

        if (k - 1) % samples_per_loop == 0
            @show c
            # Check wether the loop is closed
            Δ = infinity_norm(samples_buffer[1], xk)
            if Δ < PathTracking.tol(tracker)
                return (:success, c)
            end

            c += 1
            if length(samples_buffer) < c * samples_per_loop + 1
                extendsamples_buffer!(cache, samples_per_loop)
            end
        end

        Θk₋₁ = Θk
        xk₋₁ = xk
    end

    return
    (:loop_failed_winding_number_too_high, 0)
end

function extendsamples_buffer!(cache, samples_per_loop)
    for _ = 1:samples_per_loop
        push!(cache.samples, copy(cache.samples[1]))
    end
end

"""
    cycleheuristic(t, x_R, x_λR, x_λ2R, x_λ3R, L = 0.75)

This is the first heuristic when to start the endgame.
It checks whether two consectuive appromixations of the first nonvanishing
monomial of the Puiseux series expansion of x(t) at 0 agree.
For more details see "A Parallel Endgame " by Bates, Hauenstein and Sommese [^1],
page 8 and 9 and the Bertini book P. 53

[^1]: Bates, Daniel J., Jonathan D. Hauenstein, and Andrew J. Sommese. "A Parallel Endgame." Contemp. Math 556 (2011): 25-35.
"""
function cycleheuristic!(state, alg::Cauchy, xs, R, λ)
    if R < alg.cycle_time_cutoff
        return true
    end
    compute_k_over_c!(state, alg, xs, λ)

    n = length(state.k_over_c_values)
    N = alg.num_needed_for_stabilization
    if n < N + 1
        return false
    end

    G = state.k_over_c_values
    lower_bound = alg.L
    upper_bound = inv(alg.L)
    for k = n-1:n-N
        if !(lower_bound ≤ G[k] / G[k+1] ≤ upper_bound)
            return false
        end
    end

    return true
end
function compute_k_over_c!(state::CauchyState, alg::Cauchy, xs, λ)
    if length(xs) ≥ 3
        k_over_c = compute_k_over_c(λ, xs[3], xs[2], xs[1])
        push!(state.k_over_c_values, k_over_c)
    end
end
# This is G(t) of page 53 from the Bertini book
function compute_k_over_c(λ, x_R, x_λR, x_λ²R)
    v₁ = rand(Complex{Float64})
    num = v₁ * x_λR[1] - v₁ * x_λ²R[1]
    denom = v₁ * x_R[1] - v₁ * x_λR[1]

    for i=2:length(x_R)
        vᵢ = rand(Complex{Float64})
        num += vᵢ * x_λR[i] - vᵢ * x_λ²R[i]
        denom += vᵢ * x_R[i] - vᵢ * x_λR[i]
    end
    log(abs(num / denom)) / log(λ)
end

"""
    ratioheuristic(radius, samplepoints, nsamplepoints, β, K=0.5)

This is the second heuristic when to start the endgame.
It enforces that the values around the lopp do not differ radically.
This functions assumes that `samplepoints` has been collected on ``x(Re^{iθ})`` for ``θ∈[0,2π]``.
See "A Parallel Endgame " by Bates, Hauenstein and Sommese [^1], page 9 for more details.

[^1]: Bates, Daniel J., Jonathan D. Hauenstein, and Andrew J. Sommese. "A Parallel Endgame." Contemp. Math 556 (2011): 25-35.
"""
function ratioheuristic(radius, samplepoints, nsamplepoints, β, K=0.5)
    # failsafe
    if radius < 1e-14
        true
    end
    # winding number is 1 so this is safe
    m = M = infinity_norm(samplepoints[1])
    for k=2:nsamplepoints
        res = infinity_norm(samplepoints[k])
        m = min(m, res)
        M = max(M, res)
    end
    (M - m < β) || (m / M > K)
end

"""
    predict_with_cif!(x, cache, nsamples)

Predict the value of a function based on the samplepoints using the cauchy integral formula.
This function assumes that the samplepoints are sampled uniformly from a circle around 0,
since then the integral (approximated with the trapezoidal rule) is simply the mean of
all samplepoints.
"""
function predict_with_cif!(x, cache, nsamples)
    samples = cache.samples

    x .= samples[1]
    @inbounds for k=2:(nsamples)
        @. x = x + samples[k]
    end
    scale!(x, inv(nsamples))
end
