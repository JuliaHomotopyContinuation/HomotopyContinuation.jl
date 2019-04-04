export solve, Result, nresults, nfinite, nsingular, natinfinity, nfailed, nnonsingular, nreal, ntracked,
    finite, results, mapresults, failed, atinfinity, singular, nonsingular, seed,
    solutions, realsolutions, multiplicities, uniquesolutions, statistics

"""
    solve(args...; options...)::Result

The solve function takes many different arguments and options depending on your specific situation,
but in the it always returns a [`Result`](@ref) containing the result of the computations.
In the following we show the different inputs `solve` takes.

# Total Degree Homotopy

    solve(F; options...)

Solve the system `F` using a total degree homotopy. `F` can be
- `Vector{<:MultivariatePolynomials.AbstractPolynomial}` (e.g. constructed by `@polyvar`)
- A composition of polynomial systems constructed by [`compose`](@ref).
- [`AbstractSystem`](@ref) (the system has to represent a **homogeneous** polynomial system.)

Additionally if `F` has a multi-homogenous structure you can provide variable groups
to use a multi-homogenous totaldegree homotopy.

### Examples
We can solve the system ``F(x,y) = (x^2+y^2+1, 2x+3y-1)`` in the following way:
```julia
julia> @polyvar x y;
julia> solve([x^2+y^2+1, 2x+3y-1])
Result with 2 tracked paths
==================================
• 2 non-singular finite solutions (0 real)
• 0 singular finite solutions (0 real)
• 0 solutions at infinity
• 0 failed paths
• random seed: 448703
```

If your polynomial system is already homogeneous, but you would like to consider it as an affine system
you can do
```julia
@polyvar x y z
solve([x^2+y^2+z^2, 2x+3y-z], homvar=z)
```
This yields the same result as `solve([x^2+y^2+1, 2x+3y-1])`.

By exploiting the multi-homogenous structure of a polynomial system it is possible
to decrease the number of paths necessary to track.
```julia
@polyvar x y
# Use variable groups to only track 2 paths instead of 4
solve([x*y - 6, x^2 - 5], variable_groups=[(x,), (y,)])
```
To check whether a certain variable grouping is beneficial you can use the [`bezout_number`](@ref)
function.


# Start Target Homotopy

    solve(G, F, start_solutions; options...)

This constructs the homotopy ``H(x,t) = tG(x)+(1-t)F(x)`` to compute solutions of the
system `F`. `start_solutions` is a list of solutions of `G` which are tracked to solutions
of `F`.

### Example
```julia
@polyvar x y
G = [x^2-1,y-1]
F = [x^2+y^2+z^2, 2x+3y-z]
solve(G, F, [[1, 1], [-1, 1]])
```

# Parameter Homotopy
    solve(F, startsolutions; parameters, p₁, p₀, γ₁=nothing, γ₀=nothing)

Solve the parameter homotopy
```math
H(x, t) = F(x, (tγ₁p₁+(1-t)γ₀p₀) / (tγ₁+(1-t)γ₀)),
```
where `p₁` and `p₀` are a vector of parameter values for ``F`` and
`γ₁` and `γ₀` are complex numbers. If `γ₁` or `γ₀` is `nothing`, it is assumed
that `γ₁` and `γ₀` are ``1``.
The input `parameters` specifies the parameter variables of `F`
which should be considered as parameters.
Necessarily we have `length(parameters) == length(p₁) == length(p₀)`.

    solve(F, startsolutions; parameters, startparameters, targetparameters, startgamma=nothing, targetgamma=nothing)

This is a non-unicode variant where `γ₁=startparameters`, `γ₀=targetparameters`,
    `γ₁=startgamma`, γ₀=`targetgamma`.

### Example
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
solve(F, startsolutions; parameters=a, p₁=p₁, p₀=p₀)
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

* `seed::Int`: The random seed used during the computations.
* `report_progress=true`: Whether a progress bar should be printed to standard out.
* `threading=true`: Enable or disable multi-threading.
* `path_result_details=:default`: The amount of information computed in each path result. Possible values are `:minimal` (minimal details), `:default` (default) and `:extensive` (all information possible).
* `homvar::Union{Int,MultivariatePolynomials.AbstractVariable}`: This considers the *homogeneous* system `F` as an affine system which was homogenized by `homvar`. If `F` is an `AbstractSystem` `homvar` is the index (i.e. `Int`) of the homogenization variable. If `F` is an `AbstractVariables` (e.g. created by `@polyvar x`) `homvar` is the actual variable used in the system `F`.
* `system::AbstractSystem`: A constructor to assemble a [`AbstractSystem`](@ref). The default is [`SPSystem`](@ref). This constructor is only applied to the input of `solve`. The constructor is called with `system(polynomials, variables)` where `polynomials` is a vector of `MultivariatePolynomials.AbstractPolynomial`s and `variables` determines the variable ordering. If you experience significant compilation times, consider to change system to `FPSystem`.
* `homotopy::AbstractHomotopy`: A constructor to construct a [`AbstractHomotopy`](@ref) for the totaldegree and start target homotopy. The default is [`StraightLineHomotopy`](@ref). The constructor is called with `homotopy(start, target)` where `start` and `target` are homogeneous [`AbstractSystem`](@ref)s.
* `affine_tracking::Bool=false`: Indicate whether path tracking should happen in affine space rather than projective space. Currently this is only supported for parameter homotopies.

Path tracking specific options:

* `corrector::AbstractCorrector`: The corrector used during in the predictor-corrector scheme. The default is [`NewtonCorrector`](@ref).
* `max_corrector_iters=3`: The maximal number of correction steps in a single step.
* `initial_step_size=0.1`: The step size of the first step.
* `max_steps=10_000`: The maximal number of iterations the path tracker has available.
* `min_step_size =1e-14`: The minimal step size.
* `max_step_size =Inf`: The maximal step size.
* `maximal_lost_digits::Real=-(log₁₀(eps) + 3)`: The tracking is terminated if we estimate that we loose more than `maximal_lost_digits` in the linear algebra steps.
* `predictor::AbstractPredictor`: The predictor used during in the predictor-corrector scheme. The default is [`Heun`](@ref)()`.
* `max_refinement_iters=10`: The maximal number of correction steps used to refine the final value.
* `refinement_accuracy=1e-8`: The precision used to refine the final value.
* `accuracy=1e-7`: The precision used to track a value.
* `auto_scaling=true`: This only applies if we track in affine space. Automatically regauges the variables to effectively compute with a relative accuracy instead of an absolute one.

Endgame specific options:

* `at_infinity_check::Bool=true`: Whether the path tracker should stop paths going to infinity early.
* `max_step_size_endgame_start::Float64=1e-6`: The endgame only starts if the step size becomes smaller that the provided value.
* `samples_per_loop::Int=5`: To compute singular solutions Cauchy's integral formula is used. The accuracy of the solutions increases with the number of samples per loop.
* `max_winding_number::Int=12`: The maximal number of loops used in Cauchy's integral formula.
* `max_affine_norm::Float64=1e6`: A fallback heuristic to decide whether a path is going to infinity.
* `min_val_accuracy::Float64=0.001`: A tolerance used to decide whether we are in the endgame zone.
"""
function solve end

const solve_supported_keywords = [:threading, :report_progress, :path_result_details, :save_all_paths]

function solve(args...; kwargs...)
    solve_kwargs, rest = splitkwargs(kwargs, solve_supported_keywords)
    tracker, start_solutions = pathtracker_startsolutions(args...; rest...)
    solve(tracker, start_solutions; solve_kwargs...)
end

function solve(tracker::PathTracker, start_solutions; path_result_details=:default, save_all_paths=false, kwargs...)
    results = track_paths(tracker, start_solutions; path_result_details=path_result_details, save_all_paths=save_all_paths, kwargs...)
    path_jumping_check!(results, tracker, path_result_details; finite_results_only=!save_all_paths)
    Result(results, length(start_solutions), tracker.problem.seed)
end


function track_paths(tracker, start_solutions; threading=true, report_progress=true, path_result_details::Symbol=:default, save_all_paths=false)
    results = Vector{result_type(tracker)}()
    n = length(start_solutions)

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

        chunk_size = 32
        batch_tracker = BatchTracker(tracker, S, chunk_size, path_result_details, save_all_paths)
        k = 1
        while k ≤ n
            prepare_batch!(batch_tracker, k)
            ccall(:jl_threading_run, Ref{Cvoid}, (Any,), batch_tracker)

            for R in batch_tracker.results, i in 1:chunk_size
                R[i] !== nothing && push!(results, R[i])
                R[i] = nothing
            end
            k += batch_tracker.batch_size

            update_progress!(progress, results, min(k - 1, n))
        end
    else
        for (k, s) in enumerate(start_solutions)
            return_code = track!(tracker, s, 1.0)
            if save_all_paths || return_code == PathTrackerStatus.success
                push!(results, PathResult(tracker, s, k; details=path_result_details))
            end
            k % 32 == 0 && update_progress!(progress, results, k)
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
    results::Vector{Vector{Union{Nothing, R}}}
    trackers::Vector{Tracker}
    ranges::Vector{UnitRange{Int}}
    start_solutions::V
    details::Symbol
    all_paths::Bool
    batch_size::Int
end

function BatchTracker(tracker::PathTracker, start_solutions, chunk_size::Int, path_result_details::Symbol, save_all_paths::Bool)
    n = length(start_solutions)
    batch_size = chunk_size * Threads.nthreads()
    batch_results = Threads.resize_nthreads!([Vector{Union{Nothing, result_type(tracker)}}(undef, chunk_size)])
    for R in batch_results
        R .= nothing
    end

    trackers = Threads.resize_nthreads!([tracker])
    ranges = partition_work(1:min(batch_size, n), Threads.nthreads())
    BatchTracker(batch_results, trackers, ranges, start_solutions, path_result_details, save_all_paths, batch_size)
end

function prepare_batch!(batch_tracker::BatchTracker, k::Int)
    n = length(batch_tracker.start_solutions)
    nthreads = length(batch_tracker.trackers)
    partition_work!(batch_tracker.ranges, k:min(k+batch_tracker.batch_size-1, n), nthreads)
end

function (batch::BatchTracker)()
    tid = Threads.threadid()
    track_batch!(batch.results[tid], batch.trackers[tid],
                 batch.ranges[tid], batch.start_solutions, batch.details, batch.all_paths)
end
function track_batch!(results, pathtracker, range, starts, details, all_paths)
    for (i, k) in enumerate(range)
        return_code = track!(pathtracker, starts[k], 1.0)
        if all_paths || return_code == PathTrackerStatus.success
            results[i] = PathResult(pathtracker, starts[k], k; details=details)
        else
            results[i] = nothing
        end
    end
    nothing
end

"""
    path_jumping_check!(results, tracker, details)

Try to detect path jumping by comparing the winding numbers of finite results.
"""
function path_jumping_check!(results::Vector{<:PathResult}, tracker::PathTracker, details::Symbol; finite_results_only=false)
    if finite_results_only
        finite_results_indices = collect(1:length(results))
        finite_results = results
    else
        finite_results_indices = Int[]
        finite_results = Vector{eltype(results)}()
        for (i, r) in enumerate(results)
            if isfinite(r)
                push!(finite_results, r)
                push!(finite_results_indices, i)
            end
        end
    end
    !isempty(finite_results) || return results

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
                                            details=details,
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

The result of `solve`. This is a wrapper around the results of each single path ([`PathResult`](@ref)) and it contains some additional informations like
a random seed to replicate the result.
"""
struct Result{V}
    pathresults::Vector{PathResult{V}}
    tracked_paths::Int
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
    total = R.tracked_paths)
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
    ntracked(R::Result)

Returns the total number of paths tracked.
"""
ntracked(R::Result) = R.tracked_paths

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
    println(io, "Result with $(s.nonsingular + s.singular) solutions")
    println(io, "==================================")
    println(io, "• $(s.nonsingular) non-singular finite $(plural("solution", s.nonsingular)) ($(s.real_nonsingular) real)")
    println(io, "• $(s.singular) singular finite $(plural("solution", s.singular)) ($(s.real_singular) real)")
    s.atinfinity > 0 &&
        println(io, "• $(s.atinfinity) $(plural("solution", s.atinfinity)) at infinity")
    s.failed > 0 &&
        println(io, "• $(s.failed) failed $(plural("path", s.failed))")
    println(io, "• $(ntracked(x)) paths tracked")
    println(io, "• random seed: $(seed(x))")
end

function Base.show(io::IO, x::ProjectiveResult)
    s = statistics(x)
    println(io, "Result with $(s.nonsingular + s.singular) solutions")
    println(io, "==================================")
    println(io, "• $(s.nonsingular) non-singular $(plural("solution", s.nonsingular)) ($(s.real_nonsingular) real)")
    println(io, "• $(s.singular) singular finite $(plural("solution", s.singular)) ($(s.real_singular) real)")
    s.failed > 0 &&
        println(io, "• $(s.failed) failed $(plural("path", s.failed))")
    println(io, "• $(ntracked(x)) paths tracked")
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
        return ntracked(r)
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
