export solve,
    Solver,
    solver_startsolutions,
    paths_to_track,
    parameter_homotopy,
    linear_subspace_homotopy

struct SolveStats
    regular::Threads.Atomic{Int}
    regular_real::Threads.Atomic{Int}
    singular::Threads.Atomic{Int}
    singular_real::Threads.Atomic{Int}
end
SolveStats() = SolveStats(Threads.Atomic{Int}.((0, 0, 0, 0))...)

function init!(SS::SolveStats)
    SS.regular[] = SS.regular_real[] = SS.singular[] = SS.singular_real[] = 0
    SS
end

function update!(stats::SolveStats, R::PathResult)
    is_success(R) || return stats

    if is_singular(R)
        Threads.atomic_add!(stats.singular_real, Int(is_real(R)))
        Threads.atomic_add!(stats.singular, 1)
    else
        Threads.atomic_add!(stats.regular_real, Int(is_real(R)))
        Threads.atomic_add!(stats.regular, 1)
    end
    stats
end

"""
    Solver(path_tracker; seed = nothing)

A struct containing multiple copies of `path_tracker`. This contains all pre-allocated
data structures to call [`solve`]. The most convenient way to construct a `Solver` is
via [`solver_startsolutions`](@ref).
"""
struct Solver{T<:AbstractPathTracker} <: AbstractSolver
    trackers::Vector{T}
    seed::Union{Nothing,UInt32}
    stats::SolveStats
    start_system::Union{Nothing,Symbol}
end
Solver(
    tracker::AbstractPathTracker;
    seed::Union{Nothing,UInt32} = nothing,
    start_system = nothing,
) = Solver([tracker], seed, SolveStats(), start_system)

Base.show(io::IO, solver::Solver) = print(io, typeof(solver), "()")



"""
    solver_startsolutions(args...; kwargs...)

Takes the same input as [`solve`](@ref) but instead of directly solving the problem
returns a [`Solver`](@ref) struct and the start solutions.

## Example

Calling `solve(args..; kwargs...)` is equivalent to
```julia
solver, starts = solver_startsolutions(args...; kwargs...)
solve(solver, starts)
```
"""
function solver_startsolutions(
    F::AbstractVector{Expression},
    starts = nothing;
    parameters = Variable[],
    variables = setdiff(variables(F), parameters),
    variable_ordering = variables,
    variable_groups = nothing,
    kwargs...,
)
    sys = System(
        F,
        variables = variable_ordering,
        parameters = parameters,
        variable_groups = variable_groups,
    )
    solver_startsolutions(sys, starts; kwargs...)
end
function solver_startsolutions(
    F::AbstractVector{<:MP.AbstractPolynomial},
    starts = nothing;
    parameters = similar(MP.variables(F), 0),
    variables = setdiff(MP.variables(F), parameters),
    variable_ordering = variables,
    variable_groups = nothing,
    target_parameters = nothing,
    kwargs...,
)
    # handle special case that we have no parameters
    # to shift the coefficients of the polynomials to the parameters
    # this was the behaviour of HC.jl v1
    if isnothing(target_parameters) && isempty(parameters)
        sys, target_parameters = ModelKit.system_with_coefficents_as_params(
            F,
            variables = variable_ordering,
            variable_groups = variable_groups,
        )
    else
        sys = System(
            F,
            variables = variable_ordering,
            parameters = parameters,
            variable_groups = variable_groups,
        )
    end
    solver_startsolutions(sys, starts; target_parameters = target_parameters, kwargs...)
end
function solver_startsolutions(
    F::Union{System,AbstractSystem},
    starts = nothing;
    seed = rand(UInt32),
    start_system = isnothing(variable_groups(F)) ? :polyhedral : :total_degree,
    generic_parameters = nothing,
    p₁ = generic_parameters,
    start_parameters = p₁,
    p₀ = generic_parameters,
    target_parameters = p₀,
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    start_subspace = nothing,
    target_subspace = nothing,
    intrinsic = nothing,
    show_progress::Bool = true,
    kwargs...,
)
    !isnothing(seed) && Random.seed!(seed)


    used_start_system = nothing
    if start_subspace !== nothing
        if target_parameters !== nothing
            F = fix_parameters(F, target_parameters; compile = compile)
        end
        H = linear_subspace_homotopy(
            F,
            start_subspace,
            target_subspace;
            intrinsic = intrinsic,
            compile = compile,
        )
        tracker = EndgameTracker(H; kwargs...)
    else
        if target_subspace !== nothing
            F = slice(F, target_subspace, compile = compile)
        end
        if start_parameters !== nothing
            tracker = parameter_homotopy_tracker(
                F;
                start_parameters = start_parameters,
                target_parameters = target_parameters,
                compile = compile,
                kwargs...,
            )
        elseif start_system == :polyhedral
            used_start_system = :polyhedral
            tracker, starts = polyhedral(
                F;
                compile = compile,
                target_parameters = target_parameters,
                show_progress = show_progress,
                kwargs...,
            )
        elseif start_system == :total_degree
            used_start_system = :total_degree
            tracker, starts = total_degree(
                F;
                compile = compile,
                target_parameters = target_parameters,
                kwargs...,
            )
            try
                first(starts)
            catch
                throw("The number of start solutions is zero (total degree is zero).")
            end
        else
            throw(
                KeywordArgumentException(
                    :start_system,
                    start_system,
                    "Possible values are: `:polyhedral` and `:total_degree`.",
                ),
            )
        end
    end

    Solver(tracker; seed = seed, start_system = used_start_system), starts
end

function solver_startsolutions(
    G::Union{System,AbstractSystem},
    F::Union{System,AbstractSystem},
    starts = nothing;
    seed = rand(UInt32),
    tracker_options = TrackerOptions(),
    endgame_options = EndgameOptions(),
    show_progress::Bool = true,
    kwargs...,
)
    !isnothing(seed) && Random.seed!(seed)
    if F isa SlicedSystem && G isa SlicedSystem && system(F) == system(G)
        H = linear_subspace_homotopy(system(F), linear_subspace(F), linear_subspace(G))
    else
        H = start_target_homotopy(G, F; kwargs...)
    end
    tracker =
        EndgameTracker(H; tracker_options = tracker_options, options = endgame_options)

    Solver(tracker; seed = seed), starts
end

function parameter_homotopy_tracker(
    F::Union{System,AbstractSystem};
    tracker_options = TrackerOptions(),
    endgame_options = EndgameOptions(),
    kwargs...,
)
    H = parameter_homotopy(F; kwargs...)
    EndgameTracker(H; tracker_options = tracker_options, options = endgame_options)
end

"""
    parameter_homotopy(F; start_parameters, target_parameters)

Construct a [`ParameterHomotopy`](@ref). If `F` is homogeneous, then a random affine chart
is chosen (via [`AffineChartHomotopy`](@ref)).
"""
function parameter_homotopy(
    F::Union{System,AbstractSystem};
    generic_parameters = randn(ComplexF64, nparameters(F)),
    p₁ = generic_parameters,
    start_parameters = p₁,
    p₀ = generic_parameters,
    target_parameters = p₀,
    tracker_options = TrackerOptions(),
    endgame_options = EndgameOptions(),
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    kwargs...,
)
    unsupported_kwargs(kwargs)
    m, n = size(F)
    H = ParameterHomotopy(fixed(F; compile = compile), start_parameters, target_parameters)
    f = System(F)
    if is_homogeneous(f)
        vargroups = variable_groups(f)
        if vargroups === nothing
            m ≥ (n - 1) || throw(FiniteException(n - 1 - m))
            H = on_affine_chart(H)
        else
            m ≥ (n - length(vargroups)) || throw(FiniteException(n - length(vargroups) - m))
            H = on_affine_chart(H, length.(vargroups) .- 1)
        end
    else
        m ≥ n || throw(FiniteException(n - m))
    end

    H
end

"""
    linear_subspace_homotopy(F, V::LinearSubspace, W::LinearSubspace, intrinsic = nothing)

Constructs an [`IntrinsicSubspaceHomotopy`](@ref) (if `dim(V) <= codim(V)` or
`intrinsic = true`) or [`ExtrinsicSubspaceHomotopy`](@ref).
Compared to the direct constructor, this also takes care of homogeneous systems.
"""
function linear_subspace_homotopy(
    F::Union{System,AbstractSystem},
    V::LinearSubspace,
    W::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    intrinsic = nothing,
    homogeneous::Union{Nothing,Bool} = nothing,
    gamma = cis(2 * pi * randn()),
)

    is_hom = isnothing(homogeneous) ? is_homogeneous(System(F)) : homogeneous

    if dim(V) <= codim(V) || something(intrinsic, false)
        if is_linear(V) && is_linear(W) && is_hom
            IntrinsicSubspaceHomotopy(
                on_affine_chart(F; compile = compile),
                V,
                W;
                gamma = gamma,
            )
        else
            IntrinsicSubspaceHomotopy(F, V, W; compile = compile, gamma = gamma)
        end
    else
        H = ExtrinsicSubspaceHomotopy(F, V, W; compile = compile, gamma = gamma)
        if is_linear(V) && is_linear(W) && is_hom
            on_affine_chart(H)
        else
            H
        end
    end
end

function start_target_homotopy(
    G::Union{System,AbstractSystem},
    F::Union{System,AbstractSystem};
    start_parameters = nothing,
    target_parameters = nothing,
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    γ = 1.0,
    gamma = γ,
    kwargs...,
)
    unsupported_kwargs(kwargs)
    f, g = System(F), System(G)

    size(F) == size(G) || error("The provided systems don't have the same size.")
    is_homogeneous(f) == is_homogeneous(g) ||
        error("The provided systems are not both homogeneous.")
    variable_groups(f) == variable_groups(g) ||
        error("The provided systems don't decalare the same variable groups.")

    m, n = size(F)

    G = fixed(G; compile = compile)
    if !isnothing(start_parameters)
        G = FixedParameterSystem(G, start_parameters)
    end

    F = fixed(F; compile = compile)
    if !isnothing(target_parameters)
        F = FixedParameterSystem(F, target_parameters)
    end

    H = StraightLineHomotopy(G, F; gamma = gamma)
    if is_homogeneous(f)
        vargroups = variable_groups(f)
        if vargroups === nothing
            m ≥ (n - 1) || throw(FiniteException(n - 1 - m))
            H = on_affine_chart(H)
        else
            m ≥ (n - length(vargroups)) || throw(FiniteException(n - length(vargroups) - m))
            H = on_affine_chart(H, length.(vargroups) .- 1)
        end
    else
        m ≥ n || throw(FiniteException(n - m))
    end

    H
end

function solver_startsolutions(
    H::Union{Homotopy,AbstractHomotopy},
    starts = nothing;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    seed = nothing,
    show_progress::Bool = true,
    kwargs...,
)
    !isnothing(seed) && Random.seed!(seed)
    Solver(EndgameTracker(fixed(H; compile = compile)); seed = seed), starts
end

"""
    solve(f; options...)
    solve(f, start_solutions; start_parameters, target_parameters, options...)
    solve(f, start_solutions; start_subspace, target_subspace, options...)
    solve(g, f, start_solutions; options...)
    solve(homotopy, start_solutions; options...)

Solve the given problem. If only a single polynomial system `f` is given, then all
(complex) isolated solutions are computed.
If a system `f` depending on parameters together with start and target parameters is given
then a parameter homotopy is performed.
If two systems `g` and `f` with solutions of `g` are given then the solutions are tracked
during the deformation of `g` to `f`.
Similarly, for a given homotopy `homotopy` ``H(x,t)`` with solutions at ``t=1`` the solutions
at ``t=0`` are computed.
See the documentation for examples.
If the input is a *homogeneous* polynomial system, solutions on a random affine chart of
projective space are computed.

## General Options
The `solve` routines takes the following options:
* `catch_interrupt = true`: If this is `true`, the computation is gracefully stopped and a
  partial result is returned when the computation is interruped.
* `compile = $(COMPILE_DEFAULT[])`: If `true` then a `System` (resp. `Homotopy`) is compiled
  to a straight line program ([`CompiledSystem`](@ref) resp. [`CompiledHomotopy`](@ref))
  for evaluation. This induces a compilation overhead. If `false` then the generated program
  is only interpreted ([`InterpretedSystem`](@ref) resp. [`InterpretedHomotopy`](@ref)).
  This is slower than the compiled version, but does not introduce compilation overhead.
* `endgame_options`: The options and parameters for the endgame.
  See [`EndgameOptions`](@ref).
* `seed`: The random seed used during the computations. The seed is also reported in the
  result. For a given random seed the result is always identical.
* `show_progress= true`: Indicate whether a progress bar should be displayed.
* `stop_early_cb`: Here it is possible to provide a function (or any callable struct) which
  accepts a `PathResult` `r` as input and returns a `Bool`. If `stop_early_cb(r)` is `true`
  then no further paths are tracked and the computation is finished. This is only called
  for successfull paths. This is for example useful if you only want to compute one solution
  of a polynomial system. For this `stop_early_cb = _ -> true` would be sufficient.
* `threading = true`: Enable multi-threading for the computation. The number of available threads is controlled by the environment variable `JULIA_NUM_THREADS`. You can run `Julia` with `n` threads using the command `julia -t n`; e.g., `julia -t 8` for `n=8`. (Some CPUs hang when using multiple threads. To avoid this run Julia with 1 interactive thread for the REPL; e.g., `julia -t 8,1` for `n=8`. Note that some CPUs seem to let `Julia` crash when using that option.)
* `tracker_options`: The options and parameters for the path tracker.
  See [`TrackerOptions`](@ref).
* `path_repair = false`: Opt in to an aggressive first pass followed by conservative
  retracking of failed or conflicting regular paths. Repaired paths are reverse-checked
  for ordinary endgame trackers. Singular multiplicities are handled by the endgame.
  Requires a complete batch (not `iterator_only` or `stop_early_cb`). Endpoint checks
  do not prove completeness or the correspondence of every unretracked path.

## Options depending on input

If only a polynomial system is given:
* `start_system`: Possible values are `:total_degree` and `:polyhedral`. Depending on the
  choice furhter options are possible. See also [`total_degree`](@ref) and
  [`polyhedral`](@ref).

If a system `f` depending on parameters together with start parameters (or start subspace),
start solutions and *multiple* target parameters (or target subspaces) then the following
options are also available:

* `flatten`: Flatten the output of `transform_result`. This is useful for example if
   `transform_result` returns a vector of solutions, and you only want a single vector of
   solutions as the result (instead of a vector of vector of solutions).
* `transform_parameters = identity`: Transform a parameters values `p` before passing it to
  `target_parameters = ...`.
* `transform_result`: A function taking two arguments, the `result` and the
  parameters `p`. By default this returns the tuple `(result, p)`.

## Basic example

```julia-repl
julia> @var x y;

julia> F = System([x^2+y^2+1, 2x+3y-1])
System of length 2
 2 variables: x, y

 1 + x^2 + y^2
 -1 + 2*x + 3*y

julia> solve(F)
Result with 2 solutions
=======================
• 2 non-singular solutions (0 real)
• 0 singular solutions (0 real)
• 2 paths tracked
• random seed: 0x75a6a462
• start_system: :polyhedral
```
"""
function solve(
    args...;
    show_progress::Bool = true,
    threading::Bool = Threads.nthreads() > 1,
    catch_interrupt::Bool = true,
    target_parameters = nothing,
    stop_early_cb = always_false,
    # many parameter options,
    transform_result = nothing,
    transform_parameters = identity,
    flatten = nothing,
    target_subspaces = nothing,
    iterator_only::Bool = false,
    bitmask = nothing,
    path_repair::Bool = false,
    kwargs...,
)
    path_repair &&
        iterator_only &&
        throw(ArgumentError("path_repair requires a complete batch"))
    path_repair &&
        stop_early_cb !== always_false &&
        throw(ArgumentError("path_repair cannot be combined with stop_early_cb"))
    many_parameters = false
    if target_subspaces !== nothing
        many_parameters = true
        solver, starts = solver_startsolutions(
            args...;
            target_subspace = first(target_subspaces),
            show_progress = show_progress,
            kwargs...,
        )
        target_parameters = target_subspaces
    elseif target_parameters !== nothing
        # check if we have many parameters solve
        if !isa(transform_parameters(first(target_parameters)), Number)
            many_parameters = true
            solver, starts = solver_startsolutions(
                args...;
                target_parameters = transform_parameters(first(target_parameters)),
                show_progress = show_progress,
                kwargs...,
            )
        else
            solver, starts = solver_startsolutions(
                args...;
                target_parameters = target_parameters,
                show_progress = show_progress,
                kwargs...,
            )
        end
    else
        solver, starts =
            solver_startsolutions(args...; show_progress = show_progress, kwargs...)
    end
    if many_parameters
        if iterator_only
            map(target_parameters) do p
                solverᵢ = deepcopy(solver)
                target_parameters!(solverᵢ, p)
                ResultIterator(starts, solver; bitmask = bitmask)
            end
        else
            solve(
                solver,
                starts,
                target_parameters;
                show_progress = show_progress,
                threading = threading,
                catch_interrupt = catch_interrupt,
                transform_result = transform_result,
                transform_parameters = transform_parameters,
                flatten = flatten,
                path_repair = path_repair,
            )
        end
    else
        if iterator_only
            ResultIterator(starts, solver; bitmask = bitmask)
        else
            solve(
                solver,
                starts;
                stop_early_cb = stop_early_cb,
                show_progress = show_progress,
                threading = threading,
                catch_interrupt = catch_interrupt,
                path_repair = path_repair,
            )
        end
    end
end

solve(S::Solver, R::Result; kwargs...) =
    solve(S, solutions(R; only_nonsingular = true); kwargs...)
solve(S::Solver, s::AbstractVector{<:Number}; kwargs...) = solve(S, [s]; kwargs...)

function solve(
    S::Solver,
    starts;
    iterator_only::Bool = false,
    bitmask = nothing,
    threading::Bool = Threads.nthreads() > 1,
    path_repair::Bool = false,
    kwargs...,
)
    if iterator_only
        path_repair && throw(ArgumentError("path_repair requires a complete batch"))
        return ResultIterator(starts, S; bitmask = bitmask)
    else
        return solve(
            S,
            collect(starts);
            threading = threading,
            path_repair = path_repair,
            kwargs...,
        )
    end
end

function solve(
    S::Solver,
    starts::AbstractArray;
    stop_early_cb = always_false,
    show_progress::Bool = true,
    threading::Bool = Threads.nthreads() > 1,
    catch_interrupt::Bool = true,
    path_repair::Bool = false,
)
    if path_repair
        stop_early_cb === always_false ||
            throw(ArgumentError("path_repair cannot be combined with stop_early_cb"))
        return repair_solve(S, starts; threading, show_progress, catch_interrupt)
    end
    n = length(starts)
    progress = show_progress ? make_progress(n; delay = 0.3) : nothing
    init!(S.stats)
    if threading
        threaded_solve(
            S,
            starts,
            progress,
            stop_early_cb;
            catch_interrupt = catch_interrupt,
        )
    else
        serial_solve(S, starts, progress, stop_early_cb; catch_interrupt = catch_interrupt)
    end
end
(solver::Solver)(starts; kwargs...) = solve(solver, starts; kwargs...)
track(solver::Solver, s; kwargs...) = track(solver.trackers[1], s; kwargs...)

function make_progress(n::Integer; delay::Float64 = 0.0)
    desc = "Tracking $n paths... "
    barlen = min(ProgressMeter.tty_width(desc, stdout, false), 40)
    progress =
        ProgressMeter.Progress(n; dt = 0.2, desc = desc, barlen = barlen, output = stdout)
    progress.tlast += delay
    progress
end
function update_progress!(progress, stats, ntracked)
    t = time()
    if ntracked == progress.n || t > progress.tlast + progress.dt
        showvalues = make_showvalues(stats, ntracked)
        ProgressMeter.update!(progress, ntracked; showvalues = showvalues)
    end
    nothing
end
@noinline function make_showvalues(stats, ntracked)
    showvalues = (("# paths tracked", ntracked),)
    nsols = stats.regular[] + stats.singular[]
    nreal = stats.regular_real[] + stats.singular_real[]
    (
        ("# paths tracked", ntracked),
        ("# non-singular solutions (real)", "$(stats.regular[]) ($(stats.regular_real[]))"),
        ("# singular endpoints (real)", "$(stats.singular[]) ($(stats.singular_real[]))"),
        ("# total solutions (real)", "$(nsols[]) ($(nreal[]))"),
    )
end
update_progress!(::Nothing, stats, ntracked) = nothing

function serial_solve(
    solver::Solver,
    starts,
    progress = nothing,
    stop_early_cb = always_false;
    catch_interrupt::Bool = true,
)
    path_results = Vector{PathResult}()
    tracker = solver.trackers[1]
    try
        for (k, s) in enumerate(starts)
            r = track(tracker, s; path_number = k)
            push!(path_results, r)
            update!(solver.stats, r)
            update_progress!(progress, solver.stats, k)
            if is_success(r) && stop_early_cb(r)
                break
            end
        end
    catch e
        (catch_interrupt && isa(e, InterruptException)) || rethrow(e)
    end

    Result(path_results; seed = solver.seed, start_system = solver.start_system)
end
function threaded_solve(
    solver::Solver,
    S::AbstractArray,
    progress = nothing,
    stop_early_cb = always_false;
    catch_interrupt::Bool = true,
)

    N = length(S)
    path_results = Vector{PathResult}(undef, N)
    interrupted = Threads.Atomic{Bool}(false)
    finished = Threads.Atomic{Int}(0)
    next_k = Threads.Atomic{Int}(1)  # next index k to process

    # Each thread gets its own tracker
    tracker = solver.trackers[1]
    ntrackers = length(solver.trackers)
    nthr = Threads.nthreads()
    resize!(solver.trackers, nthr)
    for i = (ntrackers+1):nthr
        solver.trackers[i] = deepcopy(tracker)
    end

    # Use a lock for progress
    progress_lock = ReentrantLock()

    # threading block
    try
        Threads.@sync begin
            for lt in solver.trackers
                let tracker = lt
                    Threads.@spawn begin
                        while true
                            # get current index
                            idx = Threads.atomic_add!(next_k, 1)
                            k = idx
                            if k > N || interrupted[]
                                break
                            end
                            # track
                            r = track(tracker, S[k]; path_number = k)
                            path_results[k] = r
                            nfinished = Threads.atomic_add!(finished, 1) + 1
                            # update progress and stats in a locked state
                            lock(progress_lock) do
                                update!(solver.stats, r)
                                update_progress!(progress, solver.stats, nfinished)
                            end
                            if is_success(r) && stop_early_cb(r)
                                interrupted[] = true
                            end
                        end
                    end
                end
            end
        end
    catch e
        if (
            isa(e, InterruptException) ||
            (isa(e, TaskFailedException) && isa(e.task.exception, InterruptException))
        )
            interrupted[] = true
        end
        if !interrupted[] || !catch_interrupt
            rethrow(e)
        end
    end
    # if we got interrupted we need to remove the unassigned fields
    if interrupted[]
        assigned_results = Vector{PathResult}()
        for i in eachindex(path_results)
            if isassigned(path_results, i)
                push!(assigned_results, path_results[i])
            end
        end
        Result(assigned_results; seed = solver.seed, start_system = solver.start_system)
    else
        Result(path_results; seed = solver.seed, start_system = solver.start_system)
    end

end

function start_parameters!(solver::Solver, p)
    for tracker in solver.trackers
        start_parameters!(tracker, p)
    end
    solver
end

function target_parameters!(solver::Solver, p)
    for tracker in solver.trackers
        target_parameters!(tracker, p)
    end
    solver
end
function parameters!(solver::Solver, p, q)
    for tracker in solver.trackers
        parameters!(tracker, p, q)
    end
    solver
end

"""
    paths_to_track(f; optopms..)

Returns the number of paths tracked when calling [`solve`](@ref) with the given arguments.
"""
function paths_to_track(
    f::Union{System,AbstractSystem};
    start_system::Symbol = :polyhedral,
    kwargs...,
)
    paths_to_track(f, Val(start_system); kwargs...)
end

#############################
### Many parameter solver ###
#############################


struct ManySolveStats
    solutions::Threads.Atomic{Int}
end

function solve(
    S::Solver,
    starts,
    target_parameters;
    show_progress::Bool = true,
    threading::Bool = Threads.nthreads() > 1,
    catch_interrupt::Bool = true,
    transform_result = nothing,
    transform_parameters = nothing,
    flatten = nothing,
    path_repair::Bool = false,
)
    transform_result = something(transform_result, tuple) # (solutions ∘ first) ∘ tuple
    transform_parameters = something(transform_parameters, identity)
    flatten = something(flatten, false)

    n = length(target_parameters)

    progress = show_progress ? make_many_progress(n; delay = 0.3) : nothing

    many_solve(
        S,
        starts,
        target_parameters,
        progress,
        transform_result,
        transform_parameters,
        Val(flatten);
        catch_interrupt = catch_interrupt,
        threading = threading,
        path_repair = path_repair,
    )
end


function make_many_progress(n::Integer; delay::Float64 = 0.0)
    desc = "Solving for $n parameters... "
    barlen = min(ProgressMeter.tty_width(desc, stdout, false), 40)
    progress =
        ProgressMeter.Progress(n; dt = 0.3, desc = desc, barlen = barlen, output = stdout)
    progress.tlast += delay
    progress
end
function update_many_progress!(progress, results, k, paths_per_param; flatten::Bool)
    t = time()
    if k == progress.n || t > progress.tlast + progress.dt
        showvalues = make_many_showvalues(results, k, paths_per_param; flatten = flatten)
        ProgressMeter.update!(progress, k; showvalues = showvalues)
    end
    nothing
end
@noinline function make_many_showvalues(results, k, paths_per_param; flatten::Bool)
    if flatten
        [
            ("# parameters solved", k),
            ("# paths tracked", paths_per_param * k),
            ("# results", length(results)),
        ]
    else
        [("# parameters solved", k), ("# paths tracked", paths_per_param * k)]
    end
end
update_many_progress!(::Nothing, results, k, paths_per_param; kwargs...) = nothing

function many_solve(
    solver::Solver,
    starts,
    many_target_parameters,
    progress,
    transform_result,
    transform_parameters,
    ::Val{flatten};
    threading::Bool,
    catch_interrupt::Bool,
    path_repair::Bool = false,
) where {flatten}


    if isa(starts, TotalDegreeStartSolutionsIterator)
        @warn "Solving for many parameters with total degree start system is not recommended. Instead, one should use a two-step approach: first solve a system with generic parameters, and track its solutions to the desired parameters. See https://www.juliahomotopycontinuation.org/guides/many-systems/."
    elseif isa(starts, PolyhedralStartSolutionsIterator)
        @error "Solving for many parameters with polyhedral start system is not implemented and also not recommended. Instead, one should use a two-step approach: solve a system with generic parameters, and track its solutions to the desired parameters. See https://www.juliahomotopycontinuation.org/guides/many-systems/."
    end
    results = nothing
    m = length(starts)
    try
        for (k, q) in enumerate(many_target_parameters)
            target_parameters!(solver, transform_parameters(q))
            if path_repair
                res = solve(
                    solver,
                    collect(starts);
                    threading,
                    path_repair = true,
                    show_progress = false,
                    catch_interrupt = false,
                )
            elseif threading
                res = threaded_solve(solver, collect(starts); catch_interrupt = false)
            else
                res = serial_solve(solver, starts; catch_interrupt = false)
            end

            transformed = transform_result(res, q)
            if results === nothing
                if flatten && !(transformed isa AbstractArray)
                    throw(
                        ArgumentError(
                            "Cannot flatten arguments of type `$(typeof(transformed))`",
                        ),
                    )
                end
                results = flatten ? transformed : [transformed]
            elseif flatten
                append!(results, transformed)
            else
                push!(results, transformed)
            end
            update_many_progress!(progress, results, k, m; flatten = flatten)
        end
    catch e
        if !catch_interrupt || !(
            isa(e, InterruptException) ||
            (isa(e, TaskFailedException) && isa(e.task.exception, InterruptException))
        )
            rethrow(e)
        end
    end

    results === nothing ? Any[] : results
end


###
### Solve Repair 
###

# Optional set-level policy. All tracking still uses the ordinary Solver and
# PathResult machinery. Never modify a caller's tracker configuration.
_repair_trackers(T::EndgameTracker) = (T.tracker,)
_repair_trackers(T::PolyhedralTracker) = (T.toric_tracker, T.generic_tracker.tracker)
_repair_trackers(T::OverdeterminedTracker) = _repair_trackers(T.tracker)
_repair_trackers(T::AbstractPathTracker) =
    throw(ArgumentError("path_repair is not supported by $(typeof(T))"))

function repair_parameters!(S::Solver, stage)
    for T in S.trackers, t in _repair_trackers(T)
        p = deepcopy(t.options.parameters)
        p.β_τ = stage == 0 ? 0.8 : min(p.β_τ, stage == 1 ? 0.25 : 0.1)
        p.β_ω_p = stage == 0 ? min(p.β_ω_p, 2.0) : max(p.β_ω_p, 4.0)
        p.strict_β_τ = min(0.75p.β_τ, 0.4)
        t.options.parameters = p
    end
    S
end

# Candidate generation, not a proof of equality. A projection bounds the
# infinity-norm distance and avoids an all-pairs search on separated endpoints.
function endpoint_conflict_pairs(records)
    ids = findall(is_nonsingular, records)
    isempty(ids) && return Tuple{Int,Int}[]
    n = length(solution(records[first(ids)]))
    n == 0 && return Tuple{Int,Int}[]
    scale = [max(1.0, maximum(abs, solution(records[i]))) for i in ids]
    radius = [max(32eps(), 8accuracy(records[i])) * scale[k] for (k, i) in enumerate(ids)]
    # An uncertain endpoint is retried separately, not allowed to make the
    # projection window infinite for every other point.
    radius .= map(r -> isfinite(r) ? r : 0.0, radius)
    # Use the real or imaginary coordinate with the largest spread. A fixed
    # real projection degenerates on, for example, purely imaginary endpoints.
    coordinate, component, spread = 1, real, -Inf
    for j = 1:n, part in (real, imag)
        lo, hi = extrema(part(solution(records[i])[j]) for i in ids)
        if hi - lo > spread
            coordinate, component, spread = j, part, hi - lo
        end
    end
    projection = [component(solution(records[i])[coordinate]) for i in ids]
    order = sortperm(projection)
    pairs = Tuple{Int,Int}[]
    rmax = maximum(radius)
    for a = 1:length(order)
        i = order[a]
        for b = a+1:length(order)
            j = order[b]
            projection[j] - projection[i] > radius[i] + rmax && break
            xi, xj = solution(records[ids[i]]), solution(records[ids[j]])
            maximum(abs(xi[k] - xj[k]) for k = 1:n) <= radius[i] + radius[j] || continue
            push!(pairs, (ids[i], ids[j]))
        end
    end
    pairs
end

_repair_target(::AbstractPathTracker) = nothing
function _repair_target(T::EndgameTracker)
    H = T.tracker.homotopy
    H isa ParameterHomotopy || return nothing
    H.F isa Union{InterpretedSystem,CompiledSystem,MixedSystem} || return nothing
    H.F, H.q
end

function repair_candidates(records, tracker)
    pairs = endpoint_conflict_pairs(records)
    ids = sort!(unique!([i for pair in pairs for i in pair]))
    target = _repair_target(tracker)
    if !isempty(ids) && target !== nothing
        # Close but genuinely distinct roots must not be classified by a fixed
        # distance threshold. Reuse HC's certification to settle these cases.
        F, p = target
        C = certify(F, [solution(records[i]) for i in ids], p; show_progress = false)
        ndistinct_certified(C) == length(ids) && empty!(pairs)
    end
    selected = Set(i for pair in pairs for i in pair)
    for (i, r) in enumerate(records)
        (
            is_failed(r) ||
            (is_nonsingular(r) && (!isfinite(accuracy(r)) || accuracy(r) > 1e-8))
        ) && push!(selected, i)
    end
    sort!(collect(selected)), pairs
end

function reverse_repair_check(T::EndgameTracker, r, start; separation = Inf)
    start isa AbstractVector{<:Number} || return nothing
    is_nonsingular(r) || return nothing
    back = Tracker(deepcopy(T.tracker.homotopy); options = deepcopy(T.tracker.options))
    check = track(back, solution(r), 0.0, 1.0)
    is_success(check) || return false
    tolerance = max(32eps(), 8check.accuracy) * max(1.0, norm(start, Inf))
    isfinite(tolerance) || return false
    distance(solution(check), start, InfNorm()) <= min(tolerance, separation / 4)
end
reverse_repair_check(::AbstractPathTracker, r, start; separation = Inf) = nothing

# Bound the reverse check by the separation of distinct numeric starts. Equal
# starts describe the same path and are handled separately below.
function repair_start_separation(starts, k)
    start = starts[k]
    start isa AbstractVector{<:Number} || return Inf
    separation = Inf
    for other in starts
        other isa AbstractVector{<:Number} || continue
        length(other) == length(start) || continue
        other == start && continue
        separation = min(separation, distance(other, start, InfNorm()))
    end
    separation
end

repair_accuracy(r) = isfinite(accuracy(r)) ? accuracy(r) : Inf

function repair_solve(S::Solver, starts; threading, show_progress, catch_interrupt)
    fast = repair_parameters!(deepcopy(S), 0)
    result = solve(fast, starts; threading, show_progress, catch_interrupt)
    records = path_results(result)
    # An interrupted first pass has not produced a complete endpoint set.
    if length(records) != length(starts)
        init!(S.stats)
        foreach(r -> update!(S.stats, r), records)
        return result
    end
    try
        for stage = 1:2
            ids, _ = repair_candidates(records, first(S.trackers))
            isempty(ids) && break
            retry = repair_parameters!(deepcopy(S), stage)
            repaired =
                solve(retry, starts[ids]; threading, show_progress = false, catch_interrupt)
            for r in path_results(repaired)
                k = ids[path_number(r)]
                previous = records[k]
                r.path_number = path_number(previous)
                check = reverse_repair_check(
                    first(retry.trackers),
                    r,
                    starts[k];
                    separation = repair_start_separation(starts, k),
                )
                check === false && (r.return_code = :terminated_reverse_check)
                keep_previous =
                    is_success(previous) &&
                    (!is_success(r) || repair_accuracy(previous) < repair_accuracy(r))
                chosen, discarded = keep_previous ? (previous, r) : (r, previous)
                chosen.accepted_steps += discarded.accepted_steps
                chosen.rejected_steps += discarded.rejected_steps
                chosen.steps_eg += discarded.steps_eg
                chosen.extended_precision_used |= discarded.extended_precision_used
                records[k] = chosen
            end
            length(repaired) == length(ids) || break
        end
        unresolved, pairs = repair_candidates(records, first(S.trackers))
        # Proximity alone cannot establish that two paths coincide. Only merge
        # candidates whose equal numeric starts identify the very same path.
        parent = collect(eachindex(records))
        root(i) = parent[i] == i ? i : (parent[i] = root(parent[i]))
        for (i, j) in pairs
            starts[i] isa AbstractVector{<:Number} || continue
            starts[j] isa AbstractVector{<:Number} || continue
            starts[i] == starts[j] || continue
            a, b = root(i), root(j)
            a == b && continue
            if repair_accuracy(records[a]) <= repair_accuracy(records[b])
                parent[b] = a
            else
                parent[a] = b
            end
        end
        for i in eachindex(records)
            root(i) == i || (records[i].return_code = :terminated_path_conflict)
        end
        for i in unresolved
            r = records[i]
            if is_nonsingular(r) && (!isfinite(accuracy(r)) || accuracy(r) > 1e-8)
                r.return_code = :terminated_inaccurate_endpoint
            end
        end
    catch e
        catch_interrupt && e isa InterruptException || rethrow()
    end
    init!(S.stats)
    foreach(r -> update!(S.stats, r), records)
    Result(records; seed = S.seed, start_system = S.start_system)
end
