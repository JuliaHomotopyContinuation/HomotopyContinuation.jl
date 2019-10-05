export Solver, solver_startsolutions, solve!, solve

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
    if is_singular(R)
        stats.singular_real += is_real(R)
        stats.singular += 1
    else
        stats.regular_real += is_real(R)
        stats.regular += 1
    end
    stats
end

######################
## PathJumpingCheck ##
######################
struct PathJumpingCheck{UP<:UniquePoints}
    checkpoint::UP
    duplicate_id::Base.RefValue{Int}
    tol::Base.RefValue{Float64}
    # A vector mapping the start solutions to the index in the checkpoint
    # This is initially 1:1 but can change due to rerunning of paths
    solution_mapping::Vector{Int}
end

function PathJumpingCheck(prob::AbstractProblem, n::Int, tol::Float64)
    checkpoint = UniquePoints(
        tracking_vector_type(prob),
        (x, y) -> distance(x, y, InfNorm());
        check_real = false,
    )
    rtol = Ref(tol)
    duplicate_id = Ref(0)
    solution_mapping = Vector(1:n)
    PathJumpingCheck(checkpoint, duplicate_id, rtol, solution_mapping)
end

function init!(check::PathJumpingCheck, n::Integer)
    check.solution_mapping .= 1:n
    empty!(check.checkpoint)
    check
end

function duplicate_check(x, check::PathJumpingCheck)
    id = add!(check.checkpoint, copy(x), Val(true); tol = check.tol[])
    check.duplicate_id[] = id
    id == NOT_FOUND
end

function track_with_pathjumping_check!(
    results::Vector,
    tracker,
    S::AbstractVector,
    k::Integer,
    check::PathJumpingCheck;
    path_result_details::Symbol = :default,
    save_all_paths::Bool = false,
)
    path_number = k
    return_code = track!(tracker, S[k])

    if is_terminated_callback(return_code)
        # read out the other path to be considered
        j = check.solution_mapping[check.duplicate_id[]]

        # rerun paths and decrease max_correctors first to 2 then to 1.
        max_corrector_iters = min(3, tracker.default_ct_options.max_corrector_iters)
        accuracy = min(tracker.default_ct_options.accuracy, 1e-7)
        while max_corrector_iters > 1
            max_corrector_iters -= 1
            return_code = track!(
                tracker,
                S[k];
                accuracy = accuracy,
                max_corrector_iters = max_corrector_iters,
            )
            # Still duplicate?
            if is_terminated_callback(return_code)
                # If we still have duplicate
                # we assume that we can take the result of path j
                Rⱼ = results[j]
                if Rⱼ === nothing
                    results[k] = nothing
                else
                    Rⱼ.path_number[] = k
                    results[k] = Rⱼ
                end
                # clear the other result
                results[j] = nothing
                # rerun other path
                return_code = track!(
                    tracker,
                    S[j];
                    accuracy = accuracy,
                    max_corrector_iters = max_corrector_iters,
                )
                path_number = j
                if !is_terminated_callback(return_code)
                    check.solution_mapping[k] = j
                    break
                else
                    # TODO: This currently assumes that we can resolve the jump immediately
                    # But it can happen that we have to resolve this recursively
                    # if path j is now on the correct path, but there is a new collision
                    # from another previous path
                end
            else
                check.solution_mapping[k] = length(check.checkpoint)
                break
            end
        end
    else
        check.solution_mapping[k] = length(check.checkpoint)
    end
    return_code, path_number
end


############
## SOLVER ##
############

struct Solver{PT<:AbstractPathTracker,UP<:UniquePoints}
    trackers::PT
    stats::SolveStats
    path_jumping_check::PathJumpingCheck{UP}
end

function Solver(prob::AbstractProblem, start_solutions; kwargs...)
    # tol is only available after we constructed tracker
    path_jumping_check = PathJumpingCheck(prob, length(start_solutions), Inf)
    tracker = construct_tracker(
        prob,
        start_solutions;
        endgame_start_callback = duplicate_check,
        endgame_start_callback_state = path_jumping_check,
        kwargs...,
    )
    path_jumping_check.tol[] = accuracy(tracker)

    Solver(tracker, SolveStats(), path_jumping_check)
end

function Solver(tracker::AbstractPathTracker, checkpoint::UniquePoints)
    Solver(tracker, SolveStats(), checkpoint)
end

function solver_startsolutions(args...; kwargs...)
    supported, rest = splitkwargs(kwargs, problem_startsolutions_supported_keywords)
    prob, start_solutions = problem_startsolutions(args...; supported...)
    Solver(prob, start_solutions; rest...), start_solutions
end

accuracy(T::PathTracker) = T.options.min_accuracy
accuracy(T::OverdeterminedTracker) = accuracy(T.tracker)
accuracy(T::PolyhedralTracker) = accuracy(T.generic_tracker)

############
## solve! ##
############

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
    solve!(solver::Solver, start_solutions;
        show_progress = true,
        path_result_details = :default,
        save_all_paths = false,
        path_jumping_check = true)

Solve the problem encoded in `solver` for the given start solutions `start_solutions`.
"""
function solve!(
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
            delay = 0.3,
        )
    else
        progress = nothing
    end
    solve!(solver, start_solutions, progress, stop_early_cb; kwargs...)
end

function solve!(
    solver::Solver,
    start_solutions,
    progress::Union{Nothing,ProgressMeter.Progress},
    stop_early_cb = nothing;
    path_result_details::Symbol = :default,
    save_all_paths::Bool = false,
    path_jumping_check::Bool = true,
    threading::Bool = true,
)
    @unpack trackers, stats = solver
    tracker = trackers

    S = collect_startsolutions(start_solutions)
    n = length(S)

    init!(stats)
    init!(solver.path_jumping_check, n)
    prepare!(tracker, start_solutions)

    results = Vector{Union{Nothing,result_type(tracker)}}(undef, n)
    results .= nothing

    ntracked = 0
    n_blas_threads = single_thread_blas()

    try
        for k = 1:n
            if path_jumping_check
                return_code, path_number = track_with_pathjumping_check!(
                    results,
                    tracker,
                    S,
                    k,
                    solver.path_jumping_check;
                    path_result_details = path_result_details,
                    save_all_paths = save_all_paths,
                )
            else
                return_code = track!(tracker, S[k])
                path_number = k
            end

            ntracked = k

            if save_all_paths ||
               is_success(return_code) || is_invalid_startvalue(return_code)
                result = PathResult(
                    tracker,
                    S[path_number],
                    path_number;
                    details = path_result_details,
                )
                results[path_number] = result
                if stop_early_cb !== nothing
                    if is_success(result) && stop_early_cb(result)
                        break
                    end
                end
            end

            is_success(return_code) && update!(stats, results[path_number])
            ntracked % 32 == 0 && update_progress!(progress, ntracked, stats)
        end
        # don't print if it already got printed above
        ntracked % 32 != 0 && update_progress!(progress, ntracked, stats)
    catch e
        if !isa(e, InterruptException)
            rethrow()
        end
    end

    n_blas_threads > 1 && set_num_BLAS_threads(n_blas_threads)

    Result(
        remove_nothings(results),
        ntracked,
        seed(tracker);
        multiplicity_tol = 10 * accuracy(tracker),
    )
end

prepare!(PT::PathTracker, S) = PT
collect_startsolutions(x::AbstractVector) = x
collect_startsolutions(x) = collect(x)

function update_progress!(progress, ntracked, stats::SolveStats; finished::Bool = false)
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


"""
    solve(args...; options...)::Result

The solve function takes many different arguments and options depending on your specific situation,
but in the end it always returns a [`Result`](@ref) containing the result of the computations.
Depending on the prodived arguments different kind of homotopies are constructed. In particular
it is possible to construct the following homotopies:

* Total degree homotopy
* Polyhedral homotopy
* Parameter homotopy
* Multi-homogenous homotopy
* Start target homotopy

If the input is a *homogenous* polynomial system, solutions in projective space are computed.
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
* `threading` (default `true`): Enable or disable multi-threading. The number of threads used
  is controlled by the environment variable `JULIA_NUM_THREADS`.

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
H(x, t) = F(x, (tγ₁p₁+(1-t)γ₀p₀) / (tγ₁+(1-t)γ₀)),
```
where ``p₁`` (=`start_parameters`) and ``p₀`` (=`target_parameters`) are vectors of
parameter values for ``F`` and ``γ₁`` (=`start_gamma`) and ``γ₀`` (=`target_gamma`)
    are complex numbers.
If `start_parameters` or `target_parameters` is `nothing`, it is assumed that `γ₁` and `γ₀` are ``1``.
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


## Abstract Homotopy

    solve(H::AbstractHomotopy, start_solutions; options...)

Solve the homotopy `H` by tracking the each solution of
``H(⋅, t)`` (as provided by `start_solutions`) from ``t=1`` to ``t=0``.
Note that `H` has to be a homotopy between *homogeneous* polynomial systems.
If it should be considered as an affine system indicate which is the index
of the homogenization variable, e.g. `solve(H, startsolutions, homvar=3)`
if the third variable is the homogenization variable.


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

If your polynomial system is already homogeneous, but you would like to consider it as an affine system
you can do
```julia
@polyvar x y z
solve([x^2+y^2+z^2, 2x+3y-z], homvar=z)
```
This yields the same result as `solve([x^2+y^2+1, 2x+3y-1])`.


## Multi-homogeneous Systems

By exploiting the multi-homogenous structure of a polynomial system it is possible
to decrease the number of paths necessary to track.

```julia
@polyvar x y
# Use variable groups to only track 2 paths instead of 4
solve([x*y - 6, x^2 - 5], variable_groups=[(x,), (y,)])
```
To check whether a certain variable grouping is beneficial you can use the [`bezout_number`](@ref)
function.
"""
function solve(args...; kwargs...)
    solve_kwargs, rest = splitkwargs(kwargs, solve_supported_keywords)
    solver, start_solutions = solver_startsolutions(args...; rest...)
    solve!(solver, start_solutions; solve_kwargs...)
end
