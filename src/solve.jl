export solve, Result, nresults, nsolutions, nfinite, nsingular, nat_infinity,
    nfailed, nnonsingular, nreal, ntracked, finite, results, mapresults,
    failed, at_infinity, singular, nonsingular, seed,
    solutions, real_solutions, multiplicities, statistics, multiplicities!

"""
    solve(args...; options...)::Result

The solve function takes many different arguments and options depending on your specific situation,
but in the it always returns a [`Result`](@ref) containing the result of the computations.
In the following we show the different inputs `solve` takes and at the end we list all possible options.

# Total Degree Homotopy

    solve(F; options...)

Solve the system `F` using a start system computed from the degrees of the entries of `F`. The number of paths to track is equal to the total degree `d1⋯dn`, where `di` is the degree of the `i`th entry of `F`. `F` can be
- `Vector{<:MultivariatePolynomials.AbstractPolynomial}` (e.g. constructed by `@polyvar`)
- A composition of polynomial systems constructed by [`compose`](@ref).
- [`AbstractSystem`](@ref) (the system has to represent a **homogeneous** polynomial system.)


### Example
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

# Polyhedral Homotopy

    solve(F; start_system = :polyhedral, only_torus=false, options...)

Solve the system `F` using a start system computed from the Newton Polytopes of the entries `F`. The number of paths to track is equal to the mixed volume of the Newton Polytopes of the entries of `F`. The mixed volume is at most the total degree of `F`. `F` can be
- `Vector{<:MultivariatePolynomials.AbstractPolynomial}` (e.g. constructed by `@polyvar`)
- A composition of polynomial systems constructed by [`compose`](@ref).
- [`AbstractSystem`](@ref) (the system has to represent a **homogeneous** polynomial system.)

If `only_torus == true` then only solutions in the algebraic torus ``(ℂ\\setminus \\{0\\})^n`` will be computed.

### Example
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

# Homogeneous Systems

If `F` has is homogeneous, we return results in projective space

### Examples
```julia
julia> @polyvar x y z;
julia> solve([x^2+y^2+z^2, 2x+3y-z])
Result with 2 solutions
==================================
• 2 non-singular solutions (0 real)
• 0 singular solutions (0 real)
• 2 paths tracked
• random seed: 291729
```

If your polynomial system is already homogeneous, but you would like to consider it as an affine system
you can do
```julia
@polyvar x y z
solve([x^2+y^2+z^2, 2x+3y-z], homvar=z)
```
This yields the same result as `solve([x^2+y^2+1, 2x+3y-1])`.

# Multihomogeneous Systems

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

    solve(F, startsolutions; parameters, start_parameters, target_parameters, start_gamma=nothing, target_gamma=nothing)

Solve the parameter homotopy
```math
H(x, t) = F(x, (tγ₁p₁+(1-t)γ₀p₀) / (tγ₁+(1-t)γ₀)),
```
where ``p₁`` (=`start_parameters`) and ``p₀`` (=`target_parameters`) are vectors of parameter values for ``F`` and
``γ₁`` (=`start_gamma`) and ``γ₀`` (=`target_gamma`) are complex numbers.
If `start_parameters` or `target_parameters` is `nothing`, it is assumed that `γ₁` and `γ₀` are ``1``.
The input `parameters` specifies the variables of `F` which should be considered as parameters.
Necessarily we have `length(parameters) == length(p₁) == length(p₀)`.

    solve(F, startsolutions; parameters, p₁, p₀, γ₁=nothing, γ₀=nothing)

This is a unicode variant where `γ₁=start_parameters`, `γ₀=target_parameters`,
    `γ₁=start_gamma`, γ₀=`target_gamma`.

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
p₁ = [1, 0]
p₀ = [3im, 0.5+2im]
solve(F, startsolutions; parameters=a, start_parameters=p₁, target_parameters=p₀)
# If you like unicode this is also possible
solve(F, startsolutions; parameters=a, p₁=p₁, p₀=p₀)
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
### General options:

* `seed::Int`: The random seed used during the computations.
* `show_progress=true`: Whether a progress bar should be printed to standard out.
* `threading=true`: Enable or disable multi-threading.
* `path_result_details=:default`: The amount of information computed in each path result. Possible values are `:minimal` (minimal details), `:default` (default) and `:extensive` (all information possible).
* `homvar::Union{Int,MultivariatePolynomials.AbstractVariable}`: This considers the *homogeneous* system `F` as an affine system which was homogenized by `homvar`. If `F` is an `AbstractSystem` `homvar` is the index (i.e. `Int`) of the homogenization variable. If `F` is an `AbstractVariables` (e.g. created by `@polyvar x`) `homvar` is the actual variable used in the system `F`.
* `system::AbstractSystem`: A constructor to assemble a [`AbstractSystem`](@ref). The default is [`SPSystem`](@ref). This constructor is only applied to the input of `solve`. The constructor is called with `system(polynomials, variables)` where `polynomials` is a vector of `MultivariatePolynomials.AbstractPolynomial`s and `variables` determines the variable ordering. If you experience significant compilation times, consider to change system to `FPSystem`.
* `homotopy::AbstractHomotopy`: A constructor to construct a [`AbstractHomotopy`](@ref) for the totaldegree and start target homotopy. The default is [`StraightLineHomotopy`](@ref). The constructor is called with `homotopy(start, target)` where `start` and `target` are homogeneous [`AbstractSystem`](@ref)s.
* `affine_tracking::Bool=true`: Indicate whether path tracking should happen in affine space.
* `projective_tracking::Bool=false`: Indicate whether path tracking should happen in projective space. The flag `affine_tracking` is dominant.
* `path_jumping_check::Bool=true`: Enable a check whether one of the paths jumped to another one.

### Path tracking specific options:
* `accuracy=1e-7`: The accuracy required during the path tracking.
* `corrector::AbstractCorrector`: The corrector used during in the predictor-corrector scheme. The default is [`NewtonCorrector`](@ref).
* `initial_step_size=0.1`: The size of the first step.
* `max_corrector_iters=2`: The maximal number of correction steps in a single step.
* `max_lost_digits::Real`: The tracking is terminated if we estimate that we loose more than `max_lost_digits` in the linear algebra steps. This threshold depends on the `precision` argument.
* `max_refinement_iters=5`: The maximal number of correction steps used to refine the final value.
* `max_steps=1_000`: The maximal number of iterations the path tracker has available. Note that this changes to `10_000` for parameter homotopies.
* `max_step_size=Inf`: The maximal step size.
* `min_step_size=1e-14`: The minimal step size.
* `precision::PrecisionOption=PRECISION_FIXED_64`: The precision used for evaluating the residual in Newton's method.
* `predictor::AbstractPredictor`: The predictor used during in the predictor-corrector scheme. The default is [`Heun`](@ref)()`.
* `refinement_accuracy=1e-8`: The precision required for the final value.
* `simple_step_size_alg=false`: Use a more simple step size algorithm.
* `steps_jacobian_info_update::Int=1`: Every n-th step a linear system will be solved using a QR factorization to obtain an estimate for the condition number of the Jacobian.
* `terminate_ill_conditioned::Bool=true`: Indicates whether the path tracking should be terminated for ill-conditioned paths. A path is considerd ill-conditioned if the condition number of the Jacobian is larger than ≈1e14 or if it is larger than 1e`max_lost_digits`.

### Endgame specific options:
* `accuracy_eg::Float64=min(accuracy, 1e-5))`: It is possible to change the accuracy during the path tracking. Usually you want lower the accuracy.
* `cond_eg_start::Float64=1e4`: The endgame is only started if the condition of the Jacobian is larger than this threshold.
* `max_winding_number::Int=12`: This limits the maximal number of loops taken in applying Cauchy's formula.
* `min_cond_at_infinity::Float64=1e7`: A path is declared as going to infinity only if it's Jacobian is also larger than this threshold.
* `samples_per_loop::Int=12`: To compute singular solutions Cauchy's integral formula is used. The accuracy of the solutions increases with the number of samples per loop.
* `t_eg_start::Float64=0.1`: The endgame starts only if `t` is smaller than this threshold.
* `tol_val_inf_accurate::Float64=1e-4`: A valuation which would result in a path declared as going to infinity is only accepted if the estimated accuracy of the valuation is less than this threshold.
* `tol_val_finite_accurate::Float64=1e-3`: A valuation which would result in a proper solution is only accepted if the estimated accuracy of the valuation is less than this threshold. This is only affects solutions where the path has at some point near 0 a condition number larger than `cond_eg_start`.
It is recommended to also take a look at the [`PathTracker`](@ref) documentation for some context.

### Overdetermined system specific options:
* `overdetermined_min_accuracy=1e-5`: The minimal accuracy a non-singular solution needs to have to be considered a solution of the original system.
* `overdetermined_min_residual=1e-3`: The minimal residual a singular solution needs to have to be considered a solution of the original system.
"""
function solve end

const solve_supported_keywords = [:threading, :show_progress,
            :path_result_details, :save_all_paths, :path_jumping_check,
            # deprecated
            :report_progress]


function solve(args...; kwargs...)
    solve_kwargs, rest = splitkwargs(kwargs, solve_supported_keywords)
    tracker, start_solutions = pathtracker_startsolutions(args...; rest...)
    solve(tracker, start_solutions; solve_kwargs...)
end

function solve(tracker::Union{<:PathTracker,<:PolyhedralTracker}, start_solutions;
        path_result_details=:default, save_all_paths=false,
        path_jumping_check=true, kwargs...)
    all_results, ntracked = track_paths(tracker, start_solutions; path_result_details=path_result_details, save_all_paths=save_all_paths, kwargs...)
    if path_jumping_check
        path_jumping_check!(all_results, tracker, path_result_details; finite_results_only=!save_all_paths)
    end

    Result(all_results, ntracked, seed(tracker); multiplicity_tol=100*accuracy(tracker))
end

accuracy(T::PathTracker) = T.options.min_accuracy
accuracy(T::PolyhedralTracker) = accuracy(T.generic_tracker)

mutable struct SolveStats
    regular::Int
    regular_real::Int
    singular::Int
    singular_real::Int
end
SolveStats() = SolveStats(0,0,0,0)

function track_paths(tracker, start_solutions;
                threading=true, show_progress=true,
                path_result_details::Symbol=:default, save_all_paths=false,
                # deprecated
                report_progress=nothing)
    results = Vector{result_type(tracker)}()
    n = length(start_solutions)

    @deprecatekwarg report_progress show_progress

    if show_progress
        progress = ProgressMeter.Progress(n; dt=0.1, desc="Tracking $n paths... ",
                                    clear_output_ijulia=true, delay=0.3)
    else
        progress = nothing
    end

    stats = SolveStats()
    ntracked = 0
    try
        nthreads = Threads.nthreads()
        if threading && nthreads > 1
            batch_size = 32 * nthreads
            batches = BatchIterator(start_solutions, batch_size)
            batch_tracker = BatchTracker(tracker, batches, path_result_details, save_all_paths)
            for batch in batches
                prepare_batch!(batch_tracker, batch, ntracked)
                ccall(:jl_threading_run, Ref{Cvoid}, (Any,), batch_tracker)

                for R in batch_tracker.results
                    if R !== nothing
                        push!(results, R)
                        is_success(R) && update!(stats, R)
                    end
                end
                ntracked += length(batch)
                if batch_tracker.interrupted
                    return results, ntracked
                end

                update_progress!(progress, ntracked, stats)
            end
        else
            for (k, s) in enumerate(start_solutions)
                return_code = track!(tracker, s)
                if save_all_paths ||
                   is_success(return_code) ||
                   is_invalid_startvalue(return_code)

                    R = PathResult(tracker, s, k; details=path_result_details)
                    push!(results, R)

                    if is_success(return_code)
                        update!(stats, R)
                    end
                    ntracked = k
                end
                k % 32 == 0 && update_progress!(progress, k, stats)
            end
            # don't print if it already got printed above
            n % 32 != 0 && update_progress!(progress, n, stats)
        end
    catch e
        if isa(e, InterruptException)
            return results, ntracked
        else
            rethrow(e)
        end
    end
    results, n
end

function update!(stats::SolveStats, R::PathResult)
    if is_singular(R)
        stats.singular_real += is_real(R)
        stats.singular += 1
    else
        stats.regular_real += is_real(R)
        stats.regular += 1
    end
end

function update_progress!(progress, ntracked, stats::SolveStats; finished::Bool=false)
    progress === nothing && return nothing

    nsols = stats.regular + stats.singular
    nreal = stats.regular_real + stats.singular_real

    showvalues = (
        ("# paths tracked", ntracked),
        ("# non-singular solutions (real)", "$(stats.regular) ($(stats.regular_real))"),
        ("# singular solutions (real)", "$(stats.singular) ($(stats.singular_real))"),
        ("# total solutions (real)", "$nsols ($nreal)")
    )

    ProgressMeter.update!(progress, ntracked; showvalues=showvalues)
    nothing
end

mutable struct BatchTracker{Tracker<:PathTracker, V, R} <: Function
    trackers::Vector{Tracker}
    results::Vector{Union{Nothing, R}}
    ranges::Vector{UnitRange{Int}}
    batch::Vector{V}
    details::Symbol
    all_paths::Bool
    interrupted::Bool
    ntracked::Int
end

function BatchTracker(tracker::PathTracker, batches::BatchIterator, path_result_details::Symbol, save_all_paths::Bool)
    trackers = Threads.resize_nthreads!([tracker])
    batch_size = length(batches.batch)
    batch_results = Vector{Union{Nothing, result_type(tracker)}}(undef, batch_size)
    batch_results .= nothing
    ranges = partition_work(1:batch_size, length(trackers))
    BatchTracker(trackers, batch_results, ranges, batches.batch,
                 path_result_details, save_all_paths, false, 0)
end

function prepare_batch!(batch_tracker::BatchTracker, batch, ntracked)
    n = length(batch)
    if length(batch_tracker.results) ≠ n
        resize!(batch_tracker.results, n)
    end
    batch_tracker.results .= nothing
    batch_tracker.ntracked = ntracked
    partition_work!(batch_tracker.ranges, 1:n, length(batch_tracker.trackers))
end

function (batch::BatchTracker)()
    tid = Threads.threadid()
    try
        track_batch!(batch.results, batch.trackers[tid], batch.ranges[tid],
                 batch.batch, batch.details, batch.all_paths, batch.ntracked)
     catch e
         if isa(e, InterruptException)
             batch.interrupted = true
         else
             rethrow(e)
         end
     end
end
function track_batch!(results, pathtracker, range, starts, details, all_paths, ntracked)
    for k in range
        return_code = track!(pathtracker, starts[k])
        if all_paths ||
           is_success(return_code) ||
           is_invalid_startvalue(return_code)

            results[k] = PathResult(pathtracker, starts[k], ntracked+k; details=details)
        else
            results[k] = nothing
        end
    end
    nothing
end

"""
    path_jumping_check!(results, tracker, details)

Try to detect path jumping by comparing the winding numbers of finite results.
"""
function path_jumping_check!(results::Vector{<:PathResult}, tracker, details::Symbol; finite_results_only=false)
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

    tol = tracker.core_tracker.options.accuracy
    # find cluster of multiple solutions
    clusters = multiplicities(solution, finite_results; tol=tol)
    while true
        rerun_paths = false
        for cluster in clusters
            m = length(cluster)
            # We have to consider a couple of cases
            # 1) All m solutions in the cluster also have winding number m
            #    -> Everything OK
            # 2) We find two non-singular solutions (winding_number === nothing)
            #    in the same cluster
            #    -> Path jumping happened
            # 3) There can be multiple paths going to the same solution, i.e.,
            #    winding number | m
            #    -> OK
            # 4) It can happen that paths going to a positive dimensional component
            #    have the same endpoint and still winding number = 1.
            #
            # This is quite complicated to tear apart. Therefore we focus for
            # now only on case 2)

            jumping = false
            for i in cluster
                cond = unpack(finite_results[i].condition_jacobian, 1.0)
                w = unpack(finite_results[i].winding_number, 0)
                if w ≤ 1 && cond < 1e10
                    jumping = true
                    break
                end
            end


            if jumping
                rerun_paths = true
                # rerun
                for i in cluster
                    rᵢ = finite_results[i]
                    new_rᵢ = track(tracker, start_solution(rᵢ);
                                            path_number=rᵢ.path_number,
                                            details=details,
                                            accuracy=min(1e-7, accuracy(tracker.core_tracker)),
                                            max_corrector_iters=1)
                    finite_results[i] = new_rᵢ
                    results[finite_results_indices[i]] = new_rᵢ
                end
            end
        end

        rerun_paths || break

        prev_clusters = clusters
        clusters = multiplicities(solution, finite_results; tol=tol)
        if clusters == prev_clusters
            break
        end
    end

    results
end
function path_jumping_check!(results::Vector{<:PathResult}, tracker::PolyhedralTracker, details::Symbol; kwargs...)
    path_jumping_check!(results, tracker.generic_tracker, details; kwargs...)
end
