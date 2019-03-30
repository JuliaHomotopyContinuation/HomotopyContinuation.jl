export solve, Result, nresults, nfinite, nsingular, natinfinity, nfailed, nnonsingular, nreal,
    finite, results, mapresults, failed, atinfinity, singular, nonsingular, seed,
    solutions, realsolutions, multiplicities, uniquesolutions, statistics

"""
    solve(F; options...)

Solve the system `F` using a total degree homotopy. `F` can be
- `Vector{<:MultivariatePolynomials.AbstractPolynomial}` (e.g. constructed by `@polyvar`)
- [`AbstractSystem`](@ref) (the system has to represent a **homogeneous** polynomial system.)

### Example
Assume we want to solve the system ``F(x,y) = (x^2+y^2+1, 2x+3y-1)``.
```julia
@polyvar x y
solve([x^2+y^2+1, 2x+3y-1])
```
If you polynomial system is already homogeneous, but you would like to consider it as an affine system
you can do
```julia
@polyvar x y z
solve([x^2+y^2+z^2, 2x+3y-z], homvar=z)
```
This would result in the same result as `solve([x^2+y^2+1, 2x+3y-1])`.

To solve ``F`` by a custom `AbstractSystem` you can do
```julia
@polyvar x y z
# The system `F` has to be homgoenous system
F = SPSystem([x^2+y^2+z^2, 2x+3y-z]) # SPSystem <: AbstractSystem
# To solve the original affine system we have to tell that the homogenization variable has index 3
solve(F, homvar=3)
```
or equivalently (in this case) by
```julia
solve([x^2+y^2+z^2, 2x+3y-z], system=SPSystem)
```

# Start Target Homotopy

    solve(G, F, start_solutions; options...)

Solve the system `F` by tracking the each provided solution of
`G` (as provided by `start_solutions`).

### Example
```julia
@polyvar x y
G = [x^2-1,y-1]
F = [x^2+y^2+z^2, 2x+3y-z]
solve(G, F, [[1, 1], [-1, 1]])
```

# Parameter Homotopy
    solve(F::Vector{<:MultivariatePolynomials.AbstractPolynomial},
        startsolutions; parameters::Vector{<:MP.AbstractVariable}, p₁, p₀, γ₁=nothing, γ₀=nothing)

Solve the parameter homotopy
```math
H(x, t) = F(x, (tγ₁p₁+(1-t)γ₀p₀) / (tγ₁+(1-t)γ₀))
```,
where `p₁` and `p₀` are a vector of parameter values for ``F`` and
`γ₁` and `γ₀` are complex numbers. If `γ₁` or `γ₀` is `nothing`, it is assumed
that `γ₁` and `γ₀` are ``1``.
The input `parameters` specifies the parameter variables of `F`
which should be considered as parameters.
Neccessarily, ``length(parameters) == length(p₁) == length(p₀)``.

    solve(F::Vector{<:MultivariatePolynomials.AbstractPolynomial},
            startsolutions; parameters::Vector{<:MP.AbstractVariable},
            startparameters, targetparameters,
            startgamma=randn(ComplexF64), targetgamma=randn(ComplexF64))

This is a non-unicode variant where `γ₁=start_parameters`, `γ₀=target_parameters`,
    `γ₁=start_gamma`, γ₀=`target_gamma`.

## Example
We want to solve a parameter homotopy ``H(x,t) := F(x; t[1, 0]+(1-t)[2, 4])`` where
```math
F(x; a) := (x₁^2-a₁, x₁x₂-a₁+a₂)
```
and let's say we are only intersted in tracking of ``[1,1]``.
This can be accomplished as follows
```julia
@polyvar x[1:2] a[1:2]
F = [x[1]^2-a[1], x[1]*x[2]-a[1]+a[2]]
startsolutions = [[1, 1]]
solve(F, startsolutions, parameters=a, p₁=p₁, p₀=p₀)
# If you don't like unicode this is also possible
solve(F, startsolutions, parameters=a, startparameters=p₁, targetparameters=p₀)
```

# Abstract Homotopy

    solve(H::AbstractHomotopy, start_solutions; options...)

Solve the homotopy `H` by tracking the each solution of
``H(⋅, t)`` (as provided by `start_solutions`) from ``t=1`` to ``t=0``.
Note that `H` has to be a homotopy between *homogeneous* polynomial systems.
If it should be considered as an affine system indicate which is the index
of the homogenization variable, e.g. `solve(H, startsolutions, homvar=3)`
if the third variable is the homogenization variable.


# Options
General options:

* `system::AbstractSystem`: A constructor to assemble a [`AbstractSystem`](@ref). The default is [`SPSystem`](@ref). This constructor is only applied to the input of `solve`. The constructor is called with `system(polynomials, variables)` where `polynomials` is a vector of `MultivariatePolynomials.AbstractPolynomial`s and `variables` determines the variable ordering. If you experience significant compilation times, consider to change system to `FPSystem`.
* `homotopy::AbstractHomotopy`: A constructor to construct a [`AbstractHomotopy`](@ref). The default is [`StraightLineHomotopy`](@ref). The constructor is called with `homotopy(start, target)` where `start` and `target` are homogeneous [`AbstractSystem`](@ref)s.
* `seed::Int`: The random seed used during the computations.
* `homvar::Union{Int,MultivariatePolynomials.AbstractVariable}`: This considers the *homogeneous* system `F` as an affine system which was homogenized by `homvar`. If `F` is an `AbstractSystem` `homvar` is the index (i.e. `Int`) of the homogenization variable. If `F` is an `AbstractVariables` (e.g. created by `@polyvar x`) `homvar` is the actual variable used in the system `F`.
* `endgame_start=0.1`: The value of `t` for which the endgame is started.
* `report_progress=true`: Whether a progress bar should be printed to `STDOUT`.
* `threading=true`: Enable or disable multi-threading.

Pathtracking specific:
* `corrector::AbstractCorrector`: The corrector used during in the predictor-corrector scheme. The default is [`NewtonCorrector`](@ref).
* `max_corrector_iters=2`: The maximal number of correction steps in a single step.
* `accuracy=1e-7`: The accuracy used to track a value.
* `predictor::AbstractPredictor`: The predictor used during in the predictor-corrector scheme. The default is [`Heun`](@ref).
* `max_refinement_iters=max_corrector_iters`: The maximal number of correction steps used to refine the final value.
* `refinement_accuracy=1e-8`: The precision used to refine the final value.
* `initial_step_size=0.1`: The initial step size for the predictor.
* `min_step_size=1e-14`: The minimal step size. If the size of step is below this the path is considered failed.
* `max_steps=1000`: The maximal number of steps per path.

Endgame specific options
* `cauchy_loop_closed_tolerance=1e-3`: The tolerance for which is used to determine whether a loop is closed. The distance between endpoints is normalized by the maximal difference between any point in the loop and the starting point.
* `cauchy_samples_per_loop=6`: The number of samples used to predict an endpoint. A higher number of samples should result in a better approximation. Note that the error should be roughly ``t^n`` where ``t`` is the current time of the loop and ``n`` is `cauchy_samples_per_loop`.
* `egtol=1e-10`: This is the tolerance necessary to declare the endgame converged.
* `maxnorm=1e5`: If our original problem is affine we declare a path at infinity if the infinity norm with respect to the standard patch is larger than `maxnorm`.
* `maxwindingnumber=15`: The maximal windingnumber we try to find using Cauchys integral formula.
* `max_extrapolation_samples=4`: During the endgame a Richardson extrapolation is used to improve the accuracy of certain approximations. This is the maximal number of samples used for this.
* `minradius=1e-15`: A path is declared false if the endgame didn't finished until then.
* `sampling_factor=0.5`: During the endgame we approach ``0`` by the geometric series ``h^kR₀`` where ``h`` is `sampling_factor` and `R₀` the endgame start provided in `runendgame`.
* `maxiters_per_step=100`: The maximal number of steps bewtween two samples.
"""
function solve end

function solve(args...; threading=true, report_progress=true, kwargs...)
    tracker, start_solutions = pathtracker_startsolutions(args...; kwargs...)
    solve(tracker, start_solutions; threading=threading, report_progress=report_progress)
end

function solve(tracker::PathTracker, start_solutions; threading=true, report_progress=true, details_level=1)
    results = Vector{result_type(tracker)}(undef, length(start_solutions))
    track_paths!(results, tracker, start_solutions, threading, report_progress, details_level)
    path_jumping_check!(results, tracker, details_level)
    Result(results, tracker.problem.seed)
end

function track_paths!(results, tracker, start_solutions, threading, report_progress, details_level)
    n = length(results)

    if report_progress
        progress = ProgressMeter.Progress(n, 0.1, "Tracking $n paths... ")
    else
        progress = nothing
    end

    nthreads = Threads.nthreads()
    if threading && nthreads > 1
        # TODO: We can probably also do this better, but for now we have to collect
        # to support indexing
        S = collect(start_solutions)

        batch_size = 32 * nthreads
        ranges = partition_work(1:min(batch_size, n), nthreads)
        trackers = Threads.resize_nthreads!([tracker])
        batch_tracker = BatchTracker(results, trackers, ranges, S, details_level)

        k = 1
        while k ≤ n
            partition_work!(batch_tracker.ranges, k:min(k+batch_size-1, n), nthreads)
            ccall(:jl_threading_run, Ref{Cvoid}, (Any,), batch_tracker)
            k += batch_size
            update_progress!(progress, results, min(k - 1, n))
        end
    else
        for (k, s) in enumerate(start_solutions)
            results[k] = track(tracker, s, 1.0; path_number=k, details_level=details_level)
            k % 16 == 0 && update_progress!(progress, results, k)
        end
    end
    results
end

function update_progress!(progress, results, N)
    ProgressMeter.update!(progress, N, showvalues=((:tracked, N),))
    nothing
end
update_progress!(::Nothing, results, N) = nothing

mutable struct BatchTracker{Tracker<:PathTracker, V, R} <: Function
    results::Vector{R}
    trackers::Vector{Tracker}
    ranges::Vector{UnitRange{Int}}
    start_solutions::V
    details_level::Int
end

function (batch::BatchTracker)()
    tid = Threads.threadid()
    track_batch!(batch.results, batch.trackers[tid],
                 batch.ranges[tid], batch.start_solutions, batch.details_level)
end
function track_batch!(results, pathtracker, range, starts, details_level)
    for k in range
        results[k] = track(pathtracker, starts[k], 1.0; path_number=k, details_level=details_level)
    end
    results
end

"""
    path_jumping_check!(results, tracker, details_level)

Try to detect path jumping by comparing the winding numbers of finite results.
"""
function path_jumping_check!(results::Vector{<:PathResult}, tracker::PathTracker, details_level::Int)
    finite_results_indices = Int[]
    finite_results = Vector{eltype(results)}()
    for (i, r) in enumerate(results)
        if isfinite(r)
            push!(finite_results, r)
            push!(finite_results_indices, i)
        end
    end
    tol = tracker.core_tracker.options.refinement_accuracy
    clusters = multiplicities(solution, finite_results; tol=tol)
    while true
        for cluster in clusters
            m = length(cluster)
            all_same_winding_number = true
            for i in m
                if unpack(finite_results[i].winding_number, 1) ≠ m
                    all_same_winding_number = false
                    break
                end
            end
            if !all_same_winding_number
                # rerun
                for i in cluster
                    rᵢ = finite_results[i]
                    new_rᵢ = track(tracker, start_solution(rᵢ), 1.0; path_number=rᵢ.path_number,
                                            details_level=details_level,
                                            accuracy=min(1e-8, accuracy(tracker.core_tracker)),
                                            max_corrector_iters=1)
                    finite_results[i] = new_rᵢ
                    results[finite_results_indices[i]] = new_rᵢ

                end
            end
        end
        prev_clusters = clusters
        clusters = multiplicities(solution, finite_results; tol=tol)
        if clusters == prev_clusters
            break
        end
    end

    results
end


"""
    Result{V<:AbstractVector}

The result of `solve`.
"""
struct Result{V}
    pathresults::Vector{PathResult{V}}
    seed::Int
end

Base.length(r::Result) = length(r.pathresults)
Base.getindex(r::Result, I) = getindex(r.pathresults, I)

Base.iterate(r::Result) = iterate(r.pathresults)
Base.iterate(r::Result, state) = iterate(r.pathresults, state)
Base.lastindex(r::Result) = lastindex(r.pathresults)
Base.eltype(r::Type{Result{V}}) where {V} = PathResult{V}

const Results = Union{Result, Vector{<:PathResult}}
const ProjectiveResult = Result{<:PVector}
"""
    nresults(result; onlyreal=false, realtol=1e-6, onlynonsingular=false, singulartol=1e14, onlyfinite=true)

The number of solutions which satisfy the corresponding predicates.

## Example
```julia
result = solve(F)
# Get all non-singular results where all imaginary parts are smaller than 1e-8
nresults(result, onlyreal=true, realtol=1e-8, onlynonsingular=true)
```
"""
function nresults(R::Results; onlyreal=false, realtol=1e-6,
    onlynonsingular=false, onlysingular=false, singulartol=1e14, onlyfinite=true)
    count(R) do r
        (!onlyreal || isreal(r, realtol)) &&
        (!onlynonsingular || isnonsingular(r, singulartol)) &&
        (!onlysingular || issingular(r, singulartol)) &&
        (!onlyfinite || isfinite(r) || isprojective(r))
    end
end

"""
    statistics(R::Result; onlyreal=false, realtol=1e-6,
        onlynonsingular=false, onlysingular=false, singulartol=1e14)

Statistic about the number of (real) singular and non-singular solutions etc. Returns a named tuple with the statistics.

## Example
```julia
julia> statistics(solve([x^2+y^2-2, 2x+3y-1]))
(nonsingular = 2, singular = 0, real_nonsingular = 2, real_singular = 0, real = 2, atinfinity = 0, failed = 0, total = 2)
"""
function statistics(R::Results, onlyreal=false, realtol=1e-6,
    onlynonsingular=false, onlysingular=false, singulartol=1e14)

    failed = atinfinity = nonsingular = singular = real_nonsingular = real_singular = 0

    for r in R
        if isfailed(r)
            failed += 1
        elseif issingular(r, singulartol)
            if isreal(r, realtol)
                real_singular += 1
            end
            singular += 1
        elseif !isprojective(r) && !isfinite(r)
            atinfinity += 1
        else # finite, nonsingular
            if isreal(r, realtol)
                real_nonsingular += 1
            end
            nonsingular += 1
        end
    end
    (nonsingular = nonsingular,
    singular = singular,
    real_nonsingular = real_nonsingular,
    real_singular = real_singular,
    real = real_nonsingular + real_singular,
    atinfinity = atinfinity,
    failed = failed,
    total = length(R))
end

"""
    nfinite(affineresult)

The number of finite solutions.
"""
nfinite(R::Results) = count(isfinite, R)

"""
    nsingular(result; tol=1e10)

The number of singular solutions. A solution is considered singular
if its windingnumber is larger than 1 or the condition number is larger than `tol`.
"""
nsingular(R::Results; tol=1e10) = count(r -> issingular(r, tol), R)

"""
    natinfinity(result)

The number of solutions at infinity.
"""
natinfinity(R::Results) = count(isatinfinity, R)

"""
    nafailed(result)

The number of failed paths.
"""
nfailed(R::Results) = count(isfailed, R)

"""
    nnonsingular(result; tol=1e-10)

The number of non-singular solutions.
"""
nnonsingular(R::Result; tol = 1e10) = count(r -> isnonsingular(r, tol), R)

"""
    nreal(result; tol=1e-6)

The number of real solutions where all imaginary parts of each solution
are smaller than `tol`.
"""
nreal(R::Results; tol = 1e-6) = count(r -> isreal(r, tol), R)


"""
    seed(result)

The random seed used in the computation.
"""
seed(result::Result) = result.seed



# Filtering
"""
    results(result; onlyreal=false, realtol=1e-6, onlynonsingular=false, onlysigular=false, singulartol=1e14, onlyfinite=true)

Return all `PathResult`s for which the given conditions apply.

## Example

```julia
R = solve(F)

# This gives us all PathResults considered non-singular and real (but still as a complex vector).
realsolutions = results(R, onlyreal=true, onlynonsingular=true)
```
"""
results(R::Results; kwargs...) = mapresults(identity, R; kwargs...)
results(f::Function, R::Results; kwargs...) = mapresults(f, R; kwargs...)

"""
    mapresults(f::Function, result; conditions...)

Apply the function `f` to all `PathResult`s for which the given conditions apply. For the possible
conditions see [`results`](@ref).

## Example
```julia
# This gives us all solutions considered real (but still as a complex vector).
realsolutions = mapresults(solution, R, onlyreal=true)
```
"""
function mapresults(f::Function, R::Results;
    onlyreal=false, realtol=1e-6, onlynonsingular=false, onlysingular=false, singulartol=1e14,
    onlyfinite=true)
    [f(r) for r in R if
        (!onlyreal || isreal(r, realtol)) &&
        (!onlynonsingular || isnonsingular(r, singulartol)) &&
        (!onlysingular || issingular(r, singulartol)) &&
        (!onlyfinite || isfinite(r) || isprojective(r))]
end

"""
    solutions(result; conditions...)

Return all solution (as `Vector`s) for which the given conditions apply.
For the possible `conditions` see [`results`](@ref).

## Example
```julia
julia> @polyvar x y
julia> result = solve([(x-2)y, y+x+3]);
julia> solutions(result)
[[2.0+0.0im, -5.0+0.0im], [-3.0+0.0im, 0.0+0.0im]]
```
"""
function solutions(result::Results; kwargs...)
    mapresults(solution, result; kwargs...)
end

"""
    realsolutions(result; tol=1e-6, conditions...)

Return all real solution (as `Vector`s of reals) for which the given conditions apply.
For the possible `conditions` see [`results`](@ref). Note that `onlyreal` is always `true`
and `realtol` is now `tol`.

## Example
```julia
julia> @polyvar x y
julia> result = solve([(x-2)y, y+x+3]);
julia> realsolutions(result)
[[2.0, -5.0], [-3.0, 0.0]]
```
"""
function realsolutions(result::Results; onlyreal=true, tol=1e-6, kwargs...)
    mapresults(r -> real.(solution(r)), result; onlyreal=true, realtol=tol, kwargs...)
end

"""
    nonsingular(result::Results; conditions...)

Return all `PathResult`s for which the solution is non-singular. This is just a shorthand
for `results(R; onlynonsingular=true, conditions...)`. For the possible `conditions` see [`results`](@ref).
"""
nonsingular(R::Results; kwargs...) = results(R; onlynonsingular=true, kwargs...)

"""
    singular(result::Results; conditions...)

Return all `PathResult`s for which the solution is singular. This is just a shorthand
for `results(R; onlysingular=true, conditions...)`. For the possible `conditions` see [`results`](@ref).
"""
function singular(R::Results; singulartol=1e14, tol=singulartol, kwargs...)
    results(R; onlysingular=true, singulartol=tol, kwargs...)
end


"""
    finite(result::AffineResults; conditions...)

Return all `PathResult`s for which the solution is finite. This is just a shorthand
for `results(R; onlyfinite=true, conditions...)`. For the possible `conditions` see [`results`](@ref).
"""
finite(R::Results; kwargs...) = results(R; onlyfinite=true, kwargs...)

"""
    real(result, tol=1e-6)

Get all results where the solutions are real with the given tolerance `tol`.
See [`isreal`](@ref) for details regarding the determination of 'realness'.
"""
Base.real(R::Results; tol=1e-6) = [r for r in R if isreal(r, tol)]

"""
    failed(result)

Get all results where the path tracking failed.
"""
failed(R::Results) = [r for r in R if isfailed(r)]

"""
    atinfinity(result::AffineResult)

Get all results where the solutions is at infinity.
"""
atinfinity(R::Results) = [r for r in R if isatinfinity(r)]

"""
    multiplicities(V::Results; tol=1e-6)

Returns a `Vector` of `Vector{PathResult}`s grouping the `PathResult`s whose solutions appear with multiplicities *greater* 1 in 'V'.
Two solutions are regarded as equal, when their pairwise distance is less than 'tol'.
"""
function multiplicities(results::Results; tol=1e-6)
    map(i -> results[i], multiplicities(solution, results; tol = tol))
end

"""
    uniquesolutions(R::Result; tol=1e-6, multiplicities=false, conditions...)

Return all *unique* solutions. If `multiplicities` is `true`, then
all *unique* solutions with their correspnding multiplicities as pairs `(s, m)`
where `s` is the solution and `m` the multiplicity are returned.
For the possible `conditions` see [`results`](@ref).

## Example
```julia-repl
julia> @polyvar x;
julia> uniquesolutions([(x-3)^3*(x+2)], multiplicities=true)
[([3.0+0.0im], 3), ([-2.0+0.0im], 1)]
julia> uniquesolutions([(x-3)^3*(x+2)])
[[3.0+0.0im], [-2.0+0.0im]]
```
"""
function uniquesolutions(R::Results; tol=1e-6, multiplicities=false, conditions...)
    uniquesolutions(R, Val(multiplicities); tol=tol, conditions...)
end

function uniquesolutions(R::Results, ::Val{Multiplicities}; tol=1e-6, conditions...)  where {Multiplicities}
    sols = solutions(R; conditions...)
    M = multiplicities(sols; tol=tol)
    indicator = trues(length(sols))
    uniques = map(M) do m
        for k in m
            indicator[k] = false
        end
        if Multiplicities
            (sols[m[1]], length(m))
        else
            sols[m[1]]
        end
    end
    for (k, s) in enumerate(sols)
        if indicator[k]
            if Multiplicities
                push!(uniques, (s, 1))
            else
                push!(uniques, s)
            end
        end
    end
    uniques
end


####Show function

function Base.show(io::IO, x::Result)
    s = statistics(x)
    println(io, "Result with $(length(x)) tracked paths")
    println(io, "==================================")
    println(io, "• $(s.nonsingular) non-singular finite $(plural("solution", s.nonsingular)) ($(s.real_nonsingular) real)")
    println(io, "• $(s.singular) singular finite $(plural("solution", s.singular)) ($(s.real_singular) real)")
    println(io, "• $(s.atinfinity) $(plural("solution", s.atinfinity)) at infinity")
    println(io, "• $(s.failed) failed $(plural("path", s.failed))")
    println(io, "• random seed: $(seed(x))")
end

function Base.show(io::IO, x::ProjectiveResult)
    s = statistics(x)
    println(io, "Result with $(length(x)) tracked paths")
    println(io, "==================================")
    println(io, "• $(s.nonsingular) non-singular $(plural("solution", s.nonsingular)) ($(s.real_nonsingular) real)")
    println(io, "• $(s.singular) singular finite $(plural("solution", s.singular)) ($(s.real_singular) real)")
    println(io, "• $(s.failed) failed $(plural("path", s.failed))")
    println(io, "• random seed: $(seed(x))")
end

TreeViews.hastreeview(::Result) = true
TreeViews.hastreeview(::ProjectiveResult) = true
TreeViews.numberofnodes(::Result) = 7
TreeViews.numberofnodes(::ProjectiveResult) = 6
TreeViews.treelabel(io::IO, x::Result, ::MIME"application/prs.juno.inline") =
    print(io, "<span class=\"syntax--support syntax--type syntax--julia\">Result</span>")

function TreeViews.nodelabel(io::IO, x::Result, i::Int, ::MIME"application/prs.juno.inline")
    s = statistics(x)
    if i == 1
        print(io, "Paths tracked")
    elseif i == 2 && s.nonsingular > 0
        print(io, "$(s.nonsingular) finite non-singular ($(s.real_nonsingular) real)")
    elseif i == 3 && s.singular > 0
        print(io, "$(s.singular) finite singular ($(s.real_singular) real)")
    elseif i == 4 && (s.real_nonsingular+s.real_singular) > 0
        print(io, "$(s.real_nonsingular+s.real_singular) finite real")
    elseif i == 5 && s.atinfinity > 0
        print(io, "$(s.atinfinity) atinfinity")
    elseif i == 6 && s.failed > 0
        print(io, "$(s.failed) failed")
    elseif i == 7
        print(io, "Random seed used")
    end
end

function TreeViews.treenode(r::Result, i::Integer)
    s = statistics(r)
    if i == 1
        return length(r)
    elseif i == 2 && s.nonsingular > 0
        return finite(r, onlynonsingular=true)
    elseif i == 3 && s.singular > 0
        return finite(r, onlysingular=true)
    elseif i == 4 && (s.real_nonsingular+s.real_singular) > 0
        return finite(r, onlyreal = true)
    elseif i == 5 && s.atinfinity > 0
        return atinfinity(r)
    elseif i == 6 && s.failed > 0
        return failed(r)
    elseif i == 7
        return seed(r)
    end
    missing
end

function TreeViews.nodelabel(io::IO, x::ProjectiveResult, i::Int, ::MIME"application/prs.juno.inline")
    s = statistics(x)
    if i == 1
        print(io, "Paths tracked")
    elseif i == 2 && s.nonsingular > 0
        print(io, "$(s.nonsingular) non-singular ($(s.real_nonsingular) real)")
    elseif i == 3 && s.singular > 0
        print(io, "$(s.singular) singular ($(s.real_singular) real)")
    elseif i == 4 && (s.real_nonsingular + s.real_singular) > 0
        print(io, "$(s.real_nonsingular + s.real_singular) real solutions")
    elseif i == 5 && s.failed > 0
        print(io, "$(s.failed) failed")
    elseif i == 6
        print(io, "Random seed used")
    end
end

function TreeViews.treenode(r::ProjectiveResult, i::Integer)
    s = statistics(r)
    if i == 1
        return length(r)
    elseif i == 2 && s.nonsingular > 0
        return finite(r, onlynonsingular=true)
    elseif i == 3 && s.singular > 0
        return finite(r, onlysingular=true)
    elseif i == 4 && (s.real_nonsingular + s.real_singular) > 0
        return finite(r, onlyreal=true)
    elseif i == 5 && s.failed > 0
        return failed(r)
    elseif i == 6
        return seed(r)
    end
    missing
end

plural(singularstr, n) = n == 1 ? singularstr : singularstr * "s"
