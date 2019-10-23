export Solver, solver_startsolutions, solver, solve

###########
## STATS ##
###########
mutable struct SolveStats
    regular::Int
    regular_real::Int
    singular::Int
    singular_real::Int
end
SolveStats() = SolveStats(0, 0, 0, 0)

function init!(SS::SolveStats)
    SS.regular = SS.regular_real = SS.singular = SS.singular_real = 0
    SS
end

function update!(stats::SolveStats, R::PathResult)
    is_success(R) || return stats

    if is_singular(R)
        stats.singular_real += is_real(R)
        stats.singular += 1
    else
        stats.regular_real += is_real(R)
        stats.regular += 1
    end
    stats
end

function consolidated_stats(stats::Vector{SolveStats})
    (
     regular = sum(s.regular for s in stats),
     regular_real = sum(s.regular_real for s in stats),
     singular = sum(s.singular for s in stats),
     singular_real = sum(s.singular_real for s in stats),
    )
end

######################
## PathJumpingCheck ##
######################
function path_jumping_candidates(
    results::Vector{Union{Nothing,PathResult{V}}},
    tol::Float64,
) where {V}
    index_map = Int[]
    S = Vector{V}()
    for (i, r) in enumerate(results)
        if r !== nothing
            s = r.intermediate_solution
            if s !== nothing
                push!(S, s)
                push!(index_map, i)
            end
        end
    end

    clusters = multiplicities(S; tol = tol)
    indices = Int[]
    for cluster in clusters, i in cluster
        push!(indices, index_map[i])
    end
    indices
end

############
## SOLVER ##
############

"""
    Solver(pathtracker::PathTracker)

A `Solver` is a wrapper around a given `pathtracker` to track multiple paths.
It provides on top of the given `pathtracker` parallelization and an optional path jumping check.
To construct a  `Solver` it is convenient to use the [`solver`](@ref) or
[`solver_startsolutions`](@ref) functions. Given a solver one can use the [`solve`](@ref)
function to solve a system. This struct is constructed for any call to `solve` unless
already explicitly provided.

## Example

Assume we want to solve a polynomial system repeatedly for many different values of the
parameters `p`.
The following example shows how to use a `Solver` to avoid some computational overhead
compared to naively calling [`solve`](@ref).

```julia
@polyvar x y z p[1:3]
F = [
    x + 3 + 2y + 2 * y^2 - p[1],
    (x - 2 + 5y) * z + 4 - p[2] * z,
    (x + 2 + 4y) * z + 5 - p[3] * z,
]
q = randn(ComplexF64, 3)
S = solutions(solve(subs(F, p => q)))
# create some fake parameter values
params = [randn(3) for _ = 1:1000]
# create a `Solver` to reuse for the path tracking
F_solver = solver(F; parameters = p, generic_parameters = q)
# solve the system F for all paramaters p in params
params_solutions = map(params) do p
    solutions(solve(F_solver, S; target_parameters = p))
end
```
"""
struct Solver{PT<:AbstractPathTracker}
    trackers::Vector{PT}
    stats::Vector{SolveStats}
end
Solver(tracker::AbstractPathTracker) = Solver([tracker], [SolveStats()])
function Solver(prob::AbstractProblem, start_solutions; kwargs...)
    Solver(construct_tracker(prob, start_solutions; kwargs...))
end

"""
    solver_startsolutions(args...; kwargs...)

Create a [`Solver`](@ref) and start solutions. Takes almost the same arguments as [`solve`](@ref).
"""
function solver_startsolutions(args...; kwargs...)
    supported, rest = splitkwargs(kwargs, problem_startsolutions_supported_keywords)
    prob, start_solutions = problem_startsolutions(args...; supported...)
    Solver(prob, start_solutions; rest...), start_solutions
end

"""
    solver(args...; kwargs...)

Create a [`Solver`](@ref). Takes almost the same arguments as [`solve`](@ref).
"""
solver(args...; kwargs...) = first(solver_startsolutions(args...; kwargs...))

min_accuracy(T::PathTracker) = T.options.min_accuracy
min_accuracy(T::OverdeterminedTracker) = min_accuracy(T.tracker)
min_accuracy(T::PolyhedralTracker) = min_accuracy(T.generic_tracker)

accuracy(T::PathTracker) = T.default_ct_options.accuracy
accuracy(T::OverdeterminedTracker) = accuracy(T.tracker)
accuracy(T::PolyhedralTracker) = accuracy(T.generic_tracker)

max_corrector_iters(T::PathTracker) = T.default_ct_options.max_corrector_iters
max_corrector_iters(T::OverdeterminedTracker) = max_corrector_iters(T.tracker)
max_corrector_iters(T::PolyhedralTracker) = max_corrector_iters(T.generic_tracker)

const solve_supported_keywords = [
    :path_result_details,
    :path_jumping_check,
    :save_all_paths,
    :show_progress,
    :stop_early_cb,
    :threading,
            # deprecated
    :report_progress,
]

"""
    solve(args...; options...)::Result

The solve function takes many different arguments and options depending on your specific situation,
but in the end it always returns a [`Result`](@ref) containing the result of the computations.
Depending on the prodived arguments different kind of homotopies are constructed. In particular
it is possible to construct the following homotopies:

* Total degree homotopy
* Polyhedral homotopy
* Parameter homotopy
* Multi-homogeneous homotopy
* Start target homotopy

If the input is a *homogeneous* polynomial system, solutions in projective space are computed.
Otherwise affine solutions are computed.

## Options

The `solve` routines takes many different options. In particular all options to
[`CoreTracker`](@ref) and [`PathTracker`](@ref) are allowed.
Additionally the following options are allowed:

* `affine_tracking::Bool=true`: Indicate whether path tracking should happen in affine space.
* `early_stop_cb`: Here it is possible to provide a function (or any callable struct) which
  accepts a `PathResult` `r` as input and returns a `Bool`. If `early_stop_cb(r)` is `true`
  then no further paths are tracked and the computation is finished. This is only called
  for successfull paths unless `save_all_paths` is `true`. This is for example useful
  if you only want to compute one solution of a polynomial system.
  For this `early_stop_cb = _ -> true` would be sufficient.
* `homotopy::AbstractHomotopy`: A constructor to construct a [`AbstractHomotopy`](@ref) for
  the totaldegree and start target homotopy. The default is [`StraightLineHomotopy`](@ref).
  The constructor is called with `homotopy(start, target)` where `start` and `target` are
  homogeneous [`AbstractSystem`](@ref)s.
* `homvar::Union{Int,MultivariatePolynomials.AbstractVariable}`: This considers the
  *homogeneous* system `F` as an affine system which was homogenized by `homvar`.
  If `F` is an `AbstractSystem` `homvar` is the index (i.e. `Int`) of the homogenization
  variable. If `F` is an `AbstractVariables` (e.g. created by `@polyvar x`)
  `homvar` is the actual variable used in the system `F`.
* `path_jumping_check::Bool=true`: Enable a check whether one of
  the paths jumped to another one.
* `path_result_details=:default`: The amount of information computed in each path result.
  Possible values are `:minimal` (minimal details), `:default` (default) and `:extensive`.
* `projective_tracking::Bool=false`: Indicate whether path tracking should happen in
  projective space. The flag `affine_tracking` is dominant.
* `seed`: The random seed used during the computations. The seed is also reported in the
  result. For a given random seed the result is always identical.
* `show_progress` (default `true`): Whether a progress bar should be printed to report the
  progress of the current computation.
* `system::AbstractSystem`: A constructor to assemble a [`AbstractSystem`](@ref).
  The default is [`SPSystem`](@ref). This constructor is only applied to the input of `solve`.
  The constructor is called with `system(polynomials, variables)` where `polynomials` is a
  vector of `MultivariatePolynomials.AbstractPolynomial`s and `variables` determines the
  variable ordering. If you experience significant compilation times,
  consider to change system to `FPSystem`.
* `system_scaling` (default `:equations`) Whether to apply an automatic scaling of
  the equations (:equations), of the equations and variables (`:equations_and_variables`) or
  no scaling at all (`nothing`).
* `threading` (default `true`): Enable or disable multi-threading. The number of threads used
  is controlled by the environment variable `JULIA_NUM_THREADS`.
* `variable_ordering`: Provide a custom ordering of the variables.

# Examples


## Total Degree Homotopy

    solve(F; options...)

Solve the system `F` using a start system computed from the degrees of the entries of `F`.
The number of paths to track is equal to the total degree `d₁⋯dⱼ`, where `dᵢ` is
the degree of the `i`th entry of `F`. `F` can be

- `Vector{<:MultivariatePolynomials.AbstractPolynomial}` (e.g. constructed by using the
  exported `@polyvar`)
- A composition of polynomial systems constructed by [`compose`](@ref)
- Any [`AbstractSystem`](@ref)

We can solve the system ``F(x,y) = (x^2+y^2+1, 2x+3y-1)`` in the following way:

```julia
julia> @polyvar x y;
julia> solve([x^2+y^2+1, 2x+3y-1])
Result with 2 solutions
==================================
• 2 non-singular solutions (0 real)
• 0 singular solutions (0 real)
• 2 paths tracked
• random seed: 661766
```


## Polyhedral Homotopy

    solve(F; start_system = :polyhedral, only_torus=false, options...)

Solve the system `F` using a start system computed from the Newton Polytopes of the
ntries `F`. The number of paths to track is equal to the mixed volume of the
Newton Polytopes of the entries of `F`. The mixed volume is at most the total degree of `F`.
`F` can be
- `Vector{<:MultivariatePolynomials.AbstractPolynomial}` (e.g. constructed by `@polyvar`)
- A composition of polynomial systems constructed by [`compose`](@ref). Note that the
  composition will not preserved.

If `only_torus == true` then only solutions in the algebraic torus
``(ℂ\\setminus \\{0\\})^n`` will be computed.

We can solve the system ``F(x,y) = (x^2+y^2+1, 2x+3y-1)`` in the following way:
```julia
julia> @polyvar x y;
julia> solve([x^2+y^2+1, 2x+3y-1]; start_system = :polyhedral)
Result with 2 solutions
==================================
• 2 non-singular solutions (0 real)
• 0 singular solutions (0 real)
• 2 paths tracked
• random seed: 222880
```


## Parameter Homotopy

    solve(F, startsolutions;
        parameters,
        start_parameters,
        target_parameters,
        start_gamma = nothing,
        target_gamma = nothing,
    )

Solve the parameter homotopy
```math
H(x, t) = F(x, \\frac{tγ₁p₁+(1-t)γ₀p₀}{tγ₁+(1-t)γ₀}),
```
where ``p₁`` (=`start_parameters`) and ``p₀`` (=`target_parameters`) are vectors of
parameter values for ``F`` and ``γ₁`` (=`start_gamma`) and ``γ₀`` (=`target_gamma`)
    are complex numbers.
If `start_parameters` or `target_parameters` is `nothing`, it is assumed that `γ₀=γ₁=1`.
The input `parameters` specifies the variables of `F` which should be considered as parameters.
Necessarily we have `length(parameters) == length(p₁) == length(p₀)`.

    solve(F, startsolutions; parameters, p₁, p₀, γ₁=nothing, γ₀=nothing)

This is a unicode variant where `γ₁=start_parameters`, `γ₀=target_parameters`,
    `γ₁=start_gamma`, γ₀=`target_gamma`.

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
p₁ = [1, 0]
p₀ = [3im, 0.5+2im]
solve(F, startsolutions; parameters=a, start_parameters=p₁, target_parameters=p₀)
# If you like unicode this is also possible
solve(F, startsolutions; parameters=a, p₁=p₁, p₀=p₀)
```

## Start Target Homotopy

    solve(G, F, start_solutions; options...)

This constructs the homotopy ``H(x,t) = tG(x)+(1-t)F(x)`` to compute solutions of the
system `F`. `start_solutions` is a list of solutions of `G` which are tracked to solutions
of `F`.

```julia
@polyvar x y
G = [x^2-1,y-1]
F = [x^2+y^2+z^2, 2x+3y-z]
solve(G, F, [[1, 1], [-1, 1]])
```


## Homogeneous Systems

If `F` has is homogeneous, we return results in projective space

```julia
julia> @polyvar x y z;
julia> solve([x^2+y^2+z^2, 2x+3y-z])
Result{PVector{Complex{Float64},1}} with 2 solutions
====================================================
• 2 non-singular solutions (0 real)
• 0 singular solutions (0 real)
• 2 paths tracked
• random seed: 490575
```

It your polynomial system is not homogeneous, you can homogenize it as follows
```julia
@polyvar x y
g = [x^2+y^2+1, 2x+3y-1]
f = homogenize(g)
```
It is also possible to specify the homogenizing variable.
```julia
@polyvar x y z
g = [x^2+y^2+1, 2x+3y-1]
f = homogenize(g, z)
```

If your polynomial system is already homogeneous, but you would like to consider it as an affine system
you can do
```julia
@polyvar x y z
solve([x^2+y^2+z^2, 2x+3y-z], homvar=z)
```
This yields the same result as `solve([x^2+y^2+1, 2x+3y-1])`.


## Multi-homogeneous Systems

By exploiting the multi-homogeneous structure of a polynomial system it is possible
to decrease the number of paths necessary to track.

```julia
@polyvar x y
# Use variable groups to only track 2 paths instead of 4
solve([x*y - 6, x^2 - 5], variable_groups=[(x,), (y,)])
```
To check whether a certain variable grouping is beneficial you can use the [`bezout_number`](@ref)
function.


## Abstract Homotopy

    solve(H::AbstractHomotopy, start_solutions; options...)

Solve the homotopy `H` by tracking the each solution of
``H(⋅, t)`` (as provided by `start_solutions`) from ``t=1`` to ``t=0``.
"""
function solve(args...; kwargs...)
    solve_kwargs, rest = splitkwargs(kwargs, solve_supported_keywords)
    solver, start_solutions = solver_startsolutions(args...; rest...)
    solve(solver, start_solutions; solve_kwargs...)
end

function solve(
    solver::Solver,
    start_solutions;
    show_progress::Bool = true,
    stop_early_cb = nothing,
    kwargs...,
)
    if show_progress
        n = length(start_solutions)
        progress = ProgressMeter.Progress(
            n;
            dt = 0.1,
            desc = "Tracking $n paths... ",
            clear_output_ijulia = true,
            barlen = 40,
            delay = 0.3,
        )
    else
        progress = nothing
    end
    solve(solver, start_solutions, progress, stop_early_cb; kwargs...)
end

function solve(
    solver::Solver,
    start_solutions,
    progress::Union{Nothing,ProgressMeter.Progress},
    stop_early_cb = nothing;
    path_result_details::Symbol = :default,
    save_all_paths::Bool = false,
    path_jumping_check::Bool = true,
    threading::Bool = true,
    start_parameters = nothing,
    target_parameters = nothing,
)
    @unpack trackers, stats = solver

    Threads.resize_nthreads!(trackers)
    Threads.resize_nthreads!(stats)

    if start_parameters !== nothing
        for tracker in trackers
            start_parameters!(tracker, start_parameters)
        end
    end
    if target_parameters !== nothing
        for tracker in trackers
            target_parameters!(tracker, target_parameters)
        end
    end

    S = collect_startsolutions(start_solutions)
    n = length(S)

    for i = 1:Threads.nthreads()
        prepare!(trackers[i], start_solutions)
        init!(stats[i])
    end

    results = Vector{Union{Nothing,result_type(trackers[1])}}(undef, n)
    results .= nothing

    retracked_paths = 0
    n_blas_threads = single_thread_blas()
    try
        track_parallel!(
            results,
            trackers,
            S,
            1:n,
            stop_early_cb,
            stats,
            progress;
            path_result_details = path_result_details,
            threading = threading,
        )

        if path_jumping_check
            tol = accuracy(trackers[1])
            indices = path_jumping_candidates(results, tol)
            max_correctors = min(max_corrector_iters(trackers[1]) - 1, 2)
            retracked_paths = 0
            while !isempty(indices) && max_correctors > 1
                acc = min(accuracy(trackers[1]), 1e-6)
                track_parallel!(
                    results,
                    trackers,
                    S,
                    indices;
                    threading = threading,
                    max_corrector_iters = min(max_correctors - 1, 2),
                    accuracy = acc,
                )
                retracked_paths += length(indices)
                new_indices = path_jumping_candidates(results, tol)
                if isempty(new_indices) || new_indices ⊆ indices
                    break
                end
                indices = new_indices
            end
        end
    catch e
        if !isa(e, InterruptException)
            rethrow()
        end
    end

    n_blas_threads > 1 && set_num_BLAS_threads(n_blas_threads)
    ntracked = count(!isnothing, results)
    final_results = remove_nothings(results)
    if !save_all_paths
        filter!(
            r -> r.return_code == :success || r.return_code == :terminated_invalid_startvalue,
            final_results,
        )
    end
    Result(
        final_results,
        ntracked,
        seed(trackers[1]);
        retracked_paths = retracked_paths,
        multiplicity_tol = min_accuracy(trackers[1]),
    )
end


function track_parallel!(
    results,
    trackers,
    S,
    range,
    stop_early_cb = nothing,
    stats = nothing,
    progress = nothing;
    threading::Bool = true,
    max_corrector_iters::Union{Nothing,Int} = nothing,
    accuracy::Union{Nothing,Float64} = nothing,
    path_result_details::Symbol = :default,
)
    if threading && Threads.nthreads() > 1
        nthreads = Threads.nthreads()
        ntrackeds = fill(0, nthreads)
        interrupted = Ref(false)
        last_printed = Ref(0)

        @static if VERSION > v"1.3-"
            # create jobs channel
            jobs = Channel{Int}(length(range))
            for i in range
                put!(jobs, i)
            end
            close(jobs)

            # spawn workers which will track the paths
            @sync for _ = 1:nthreads
                Threads.@spawn begin
                    try
                        tid = Threads.threadid()
                        for i in jobs
                            results[i] = track(
                                trackers[tid],
                                S[i];
                                path_number = i,
                                details = path_result_details,
                                max_corrector_iters = max_corrector_iters,
                                accuracy = accuracy,
                            )
                            if is_success(results[i]) && stop_early_cb !== nothing
                                interrupted[] = stop_early_cb(results[i])
                            end

                            stats !== nothing && update!(stats[tid], results[i])
                            ntrackeds[tid] += 1

                            if progress !== nothing
                                ntracked = sum(ntrackeds)
                                if tid == 1 && ntracked - last_printed[] - 32nthreads > 0
                                    last_printed[] = ntracked
                                    update_progress!(
                                        progress,
                                        ntracked,
                                        consolidated_stats(stats),
                                    )
                                end
                            end
                            if interrupted[]
                                break
                            end
                        end
                    catch e
                        if isa(e, InterruptException)
                            interrupted[] = true
                        else
                            rethrow()
                        end
                    end
                end
            end
        else
            Threads.@threads for i in range
                try
                    tid = Threads.threadid()
                    results[i] = track(
                        trackers[tid],
                        S[i];
                        path_number = i,
                        max_corrector_iters = max_corrector_iters,
                        details = path_result_details,
                        accuracy = accuracy,
                    )
                    if is_success(results[i]) && stop_early_cb !== nothing
                        interrupted[] = stop_early_cb(results[i])
                    end

                    stats !== nothing && update!(stats[tid], results[i])
                    ntrackeds[tid] += 1

                    if progress !== nothing
                        ntracked = sum(ntrackeds)
                        if tid == 1 && ntracked - last_printed[] - 32nthreads > 0
                            last_printed[] = ntracked
                            update_progress!(progress, ntracked, consolidated_stats(stats))
                        end
                    end
                    if interrupted[]
                        break
                    end
                catch e
                    if isa(e, InterruptException)
                        interrupted[] = true
                    else
                        rethrow()
                    end
                end
            end
        end

        if progress !== nothing && last_printed != sum(ntrackeds)
            update_progress!(progress, sum(ntrackeds), consolidated_stats(stats))
        end
    else
        ntracked = 0
        for i in range
            results[i] = track(
                trackers[1],
                S[i];
                path_number = i,
                max_corrector_iters = max_corrector_iters,
                details = path_result_details,
                accuracy = accuracy,
            )
            stats !== nothing && update!(stats[1], results[i])
            ntracked += 1
            if ntracked % 32 == 0 && stats !== nothing
                update_progress!(progress, ntracked, stats[1])
            end

            if is_success(results[i]) &&
               stop_early_cb !== nothing && stop_early_cb(results[i])
                break
            end
        end
        if ntracked % 32 != 0 && stats !== nothing
            update_progress!(progress, ntracked, stats[1])
        end
    end
    results
end

prepare!(PT::PathTracker, S) = PT
collect_startsolutions(x::AbstractVector) = x
collect_startsolutions(x) = collect(x)

function update_progress!(progress, ntracked, stats)
    progress === nothing && return nothing

    nsols = stats.regular + stats.singular
    nreal = stats.regular_real + stats.singular_real

    showvalues = (
        ("# paths tracked", ntracked),
        ("# non-singular solutions (real)", "$(stats.regular) ($(stats.regular_real))"),
        ("# singular solutions (real)", "$(stats.singular) ($(stats.singular_real))"),
        ("# total solutions (real)", "$nsols ($nreal)"),
    )

    ProgressMeter.update!(progress, ntracked; showvalues = showvalues)
    nothing
end

function remove_nothings(v::Vector{Union{Nothing,T}}) where {T}
    w = T[]
    for vᵢ in v
        !isnothing(vᵢ) && push!(w, vᵢ)
    end
    w
end
