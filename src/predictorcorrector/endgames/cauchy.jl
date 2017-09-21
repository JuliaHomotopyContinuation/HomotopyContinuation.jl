export cauchyendgame, ConvergentCluster, CauchyEndgameResult, successfull

"""
    ConvergentCluster(points, convergence_point, t)

Construct a convergent cluster at time `t` with points ``w_1, ..., w_c`` where ``c`` is the
winding number of the `convergence_point`.
See Chapter 15.6 in The Numerical solution of systems of polynomials arising in
engineering and science[^1] for more details.

[^1]: Sommese, Andrew J., and Charles W. Wampler II. The Numerical solution of systems of polynomials arising in engineering and science. World Scientific, 2005.
"""
struct ConvergentCluster{T<:Number}
    points::Vector{Vector{T}}
    convergence_point::Vector{T}
    t::Float64
end

function Base.show(io::IO, cluster::ConvergentCluster)
    println(io, typeof(cluster),":")
    println(io, "* cluster points: ", cluster.points)
    println(io, "* convergence_point: ", cluster.convergence_point)
    println(io, "* t: ", cluster.t)
end

"""
    cauchyendgame(H::AbstractHomotopy, J_H!, Hdt!, x, R, algorithm, kwargs..., pathtrackingkwargs...)

Execute a cauchy endgame for `H` starting at radius `R`. This assumes that ``H(x,t)=0``.
It is based on "A Parallel Endgame " by Bates, Hauenstein and Sommese. [^1]

## Optional arguments
* `prediction_tolerance=1e-4`:
* `geometric_series_factor=0.5`:
* `endgame_tolerance=1e-8`:
* `samples_per_loop=8`:
* `max_winding_number=8`:
* `loopclosed_tolerance=1e-3`
*  optional arguments to [`pathtracking`](@ref)

[^1]: Bates, Daniel J., Jonathan D. Hauenstein, and Andrew J. Sommese. "A Parallel Endgame." Contemp. Math 556 (2011): 25-35.
"""
function cauchyendgame(
    H::AbstractHomotopy{T},
    J_H!::Function,
    Hdt!::Function,
    x::Vector{T},
    t::Float64,
    algorithm::APCA{Val{true}};
    prediction_tolerance=1e-8,
    geometric_series_factor=0.5,
    endgame_tolerance=1e-8,
    samples_per_loop::Int=8,
    max_winding_number::Int=8,
    loopclosed_tolerance=1e-3,
    tolerance_infinity=1e-8,
    pathtrackingkwargs...
) where {T}
    λ = geometric_series_factor
    R = t
    steps = [R]
    xs = [x]
    lastprediction = Nullable(x, false)
    endgame_started = false
    k = 1
    while R > eps(Float64)
        res = trackpath(H, J_H!, Hdt!, last(xs), algorithm, R, λ*R; pathtrackingkwargs...)
        if !issuccessfull(res)
            return CauchyEndgameResult(
                get(lastprediction, last(xs)),
                :IllConditionedZone,
                k, R, λ, xs, steps, Nullable{ConvergentCluster{T}}())
        end
        R = λ*R
        push!(steps, R)
        push!(xs, res.result)

        # we need at least 4 points to start with our first heuristic
        if length(xs) < 4
            continue
        end

        if atinfinity(get(lastprediction, last(xs)), tolerance_infinity)
            return CauchyEndgameResult(
                get(lastprediction, last(xs)),
                :AtInfinity,
                k, R, λ, xs, steps, Nullable{ConvergentCluster{T}}())
        end

        if endgame_started || firstheuristic(R, λ, xs[end-3], xs[end-2], xs[end-1], xs[end])
            retcode, samples =
                loop(H, J_H!, Hdt!, last(xs), R, algorithm,
                    samples_per_loop, max_winding_number, loopclosed_tolerance,
                    endgame_tolerance, apply_heuristic=!endgame_started)
            if retcode == :Success
                endgame_started = true
                prediction = predict_with_cif(samples)
                if !isnull(lastprediction)
                    # check whether to predictions are close enough, so we can say that
                    # the endgame terminates
                    Δprediction = projectivenorm(prediction, get(lastprediction))
                    if Δprediction < prediction_tolerance
                        cluster =
                            ConvergentCluster(samples[1:samples_per_loop:end], prediction, R)
                        return CauchyEndgameResult(prediction, :Success, k, R, λ, xs, steps, Nullable(cluster))
                    end
                end
                lastprediction = Nullable(prediction)
            end
        end

        k += 1
    end

    CauchyEndgameResult(
        get(lastprediction, last(xs)),
        :MachineEpsilon, k, R, λ, xs, steps, Nullable{ConvergentCluster{T}}())
end

"""
    CauchyEndgameResult(result, returncode, iterations, R, λ, trace, nullable_cluster)

Construct a result of a `cauchyendgame`.

## Fields
* `result::Vector{T}`
* `returncode::Symbol`: `:Success`, `:IllConditionedZone` (the path tracking failed)
or `:MachineEpsilon`: (`t` got smaller than machine epsilon)
* `iterations::Int`: How many iterations of the power series was done
* `R`: The radius when the endgame started
* `λ`: The factor of the geometric series ``λ^kR``
* `trace::Vector{Vector{T}}`: The points ``x(λ^kR)`` for `k=0,1,...,iterations`
* `convergent_cluster::Nullable{ConvergentCluster{T}}`: The convergent cluster of the endgame.
"""
struct CauchyEndgameResult{T<:Number}
    result::Vector{T}
    returncode::Symbol
    iterations::Int
    R::Float64
    λ::Float64
    trace::Vector{Vector{T}}
    steps::Vector{Float64}
    convergent_cluster::Nullable{ConvergentCluster{T}}
end

function Base.show(io::IO, res::CauchyEndgameResult)
    println(io, typeof(res),":")
    println(io, "------------------------------")
    println(io, "* result: ", res.result)
    println(io, "* returncode: ", res.returncode)
    println(io, "------------------------------")
    println(io, "* iterations: ", res.iterations)
    println(io, "* R: ", res.R)
    println(io, "* λ: ", res.λ)
    println(io, "* trace: ", length(res.trace), " entries")
    println(io, "* convergent cluster: ", get(map(string, res.convergent_cluster), "---"))
end

issuccessfull(result::CauchyEndgameResult) = result.returncode == :Success

function iscauchykwarg(kwarg)
    symb = first(kwarg)
    symb == :prediction_tolerance ||
    symb == :geometric_series_factor ||
    symb == :endgame_tolerance ||
    symb == :samples_per_loop ||
    symb == :max_winding_number ||
    symb == :loopclosed_tolerance ||
    symb == :tolerance_infinity
end
cauchykwargs(kwargs) = filter(iscauchykwarg, kwargs)

atinfinity(x, tolerance_infinity) = norm(normalize(x)[1]) < tolerance_infinity


"""
    loop(H, J_H!, Hdt!, x, radius, algorithm, samples_per_loop, max_winding_number, loopclosed_tolerance, apply_heuristic=true, pathtrackingkwargs...)

Tracks the implicit defined path z(t) around the `n`-gon with vertices
``r⋅exp(i2πk/n)`` where `n=samples_per_loop`.
"""
function loop(
    H::AbstractHomotopy{T},
    J_H!,
    Hdt!,
    x,
    radius::Real,
    algorithm::APCA,
    samples_per_loop::Int,
    max_winding_number::Int,
    loopclosed_tolerance::Float64,
    endgame_tolerance::Float64;
    apply_heuristic=true,
    pathtrackingkwargs...
) where {T}
    unitroots = UnitRootsIterator(radius, samples_per_loop)
    samples = [x]
    start = first(unitroots)
    for (k, finish) in enumerate(Iterators.drop(unitroots, 1))
        pathresult = trackpath(H, J_H!, Hdt!, samples[end], algorithm, start, finish; pathtrackingkwargs...)
        if pathresult.returncode != :Success
            return (:BadBranch, samples)
        end

        # apply the second heuristic to tell us whether it makes sense to start the endgame
        if apply_heuristic && k == samples_per_loop &&
            !secondheuristic(radius, samples, β=endgame_tolerance)
            return (:HeuristicFailed, samples)
        end

        # check whether the loop is closed. This can only happen if we finished one circle.
        if k % samples_per_loop == 0
            loopnorm = projectivenorm(first(samples), pathresult.result)
            if loopnorm < loopclosed_tolerance
                break
            end
        end

        if k > max_winding_number * samples_per_loop
            return (:WindingNumberTooHigh, samples)
        end

        push!(samples, pathresult.result)
        start = finish
    end
    (:Success, samples)
end

"""
    firstheuristic(t, x_R, x_λR, x_λ2R, x_λ3R; L = 0.75)

This is the first heuristic when to start the endgame.
It checks whether two consectuive appromixations of the first nonvanishing
monomial of the Puiseux series expansion of x(t) at 0 agree.
For more details see "A Parallel Endgame " by Bates, Hauenstein and Sommese [^1],
page 8 and 9.

[^1]: Bates, Daniel J., Jonathan D. Hauenstein, and Andrew J. Sommese. "A Parallel Endgame." Contemp. Math 556 (2011): 25-35.
"""
function firstheuristic(t, λ, x_R, x_λR, x_λ2R, x_λ3R; L = 0.75)
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
    secondheuristic(t, samplepoints, β; K=0.5)

This is the second heuristic when to start the endgame.
It enforces that the values around the lopp do not differ radically.
This functions assumes that `samplepoints` has been collected on ``x(Re^{iθ})`` for ``θ∈[0,2π]``.
See "A Parallel Endgame " by Bates, Hauenstein and Sommese [^1], page 9 for more details.

[^1]: Bates, Daniel J., Jonathan D. Hauenstein, and Andrew J. Sommese. "A Parallel Endgame." Contemp. Math 556 (2011): 25-35.
"""
function secondheuristic(radius, samplepoints; β=1e-8, K=0.5)
    # failsafe
    if radius < 1e-14
        true
    end

    norms = norm.(samplepoints)
    m = minimum(norms)
    M = maximum(norms)
    (M - m < β) || (m / M > K)
end

"""
    projectivenorm(a, b)

Calculate |a-b|, but first bring `a` and `b` on the same patch.
Brings fist `a` and `b` on the same patch by finding diving `a` through it's maximum value
(w.r.t. to the absolute value) with index `i` and diving `b` through `b[i]`.
Then computes the norm of the differences.
"""
function projectivenorm(a::Vector{<:Number}, b::Vector{<:Number})
    _, i = findmax(abs2.(a))
    norm(a ./ a[i] .- b ./ b[i])
end

"""
    predict_with_cif(samplepoints)

Predict the value of a function based on the samplepoints using the cauchy integral formula.
This function assumes that the samplepoints are sampled uniformly from a circle around 0,
since then the integral (approximated with the trapezoidal rule) is simply the mean of
all samplepoints.
"""
predict_with_cif(samplepoints) = mean(samplepoints)

"""
    UnitRootsIterator(r, n)

Construct an infinite iterator which returns the `n`-th scalded unit roots, i.e.
the values ``r⋅exp(i2πk/n)`` for ``k=0,1,...``.
"""
struct UnitRootsIterator
    radius::Float64
    order::Float64
end
UnitRootsIterator(r::Real, order::Real) = UnitRootsIterator(float(r), float(order))

Base.start(::UnitRootsIterator) = 0
Base.next(loop::UnitRootsIterator, k::Int) = (loop.radius *  exp(im * 2π * k / loop.order), k + 1)
Base.done(::UnitRootsIterator, ::Int) = false
Base.iteratorsize(::UnitRootsIterator) = Base.IsInfinite()
Base.iteratoreltype(::UnitRootsIterator) = Base.HasEltype()
Base.eltype(::UnitRootsIterator) = Complex128
