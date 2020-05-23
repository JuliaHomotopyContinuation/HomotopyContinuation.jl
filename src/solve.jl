export solve, Solver, solver_startsolutions, paths_to_track

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

struct Solver{T<:AbstractPathTracker}
    trackers::Vector{T}
    seed::Union{Nothing,UInt32}
    stats::SolveStats
end
Solver(tracker::AbstractPathTracker, seed::Union{Nothing,UInt32} = nothing) =
    Solver([tracker], seed, SolveStats())

function solver_startsolutions(
    F::AbstractVector{Expression},
    starts = nothing;
    parameters = Variable[],
    variables = setdiff(variables(F), parameters),
    variable_groups = nothing,
    kwargs...,
)
    sys = System(
        F,
        variables = variables,
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
            variables = variables,
            variable_groups = variable_groups,
        )
    else
        sys = System(
            F,
            variables = variables,
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
    p₁ = nothing,
    start_parameters = p₁,
    p₀ = nothing,
    target_parameters = p₀,
    kwargs...,
)
    !isnothing(seed) && Random.seed!(seed)

    if start_parameters !== nothing
        tracker = parameter_homotopy(
            F;
            start_parameters = start_parameters,
            target_parameters = target_parameters,
            kwargs...,
        )
    elseif start_system == :polyhedral
        tracker, starts = polyhedral(F; target_parameters = target_parameters, kwargs...)
    elseif start_system == :total_degree
        tracker, starts = total_degree(F; target_parameters = target_parameters, kwargs...)
    else
        throw(KeywordArgumentException(
            :start_system,
            start_system,
            "Possible values are: `:polyhedral` and `:total_degree`.",
        ))
    end

    Solver(tracker, seed), starts
end

function solver_startsolutions(
    G::Union{System,AbstractSystem},
    F::Union{System,AbstractSystem},
    starts = nothing;
    seed = rand(UInt32),
    tracker_options = TrackerOptions(),
    endgame_options = EndgameOptions(),
    kwargs...,
)
    !isnothing(seed) && Random.seed!(seed)
    H = start_target_homotopy(G, F; kwargs...)
    tracker =
        EndgameTracker(H; tracker_options = tracker_options, options = endgame_options)

    Solver(tracker, seed), starts
end

function parameter_homotopy(
    F::Union{System,AbstractSystem};
    p₁ = nothing,
    start_parameters = p₁,
    p₀ = nothing,
    target_parameters = p₀,
    tracker_options = TrackerOptions(),
    endgame_options = EndgameOptions(),
)
    isnothing(start_parameters) && throw(UndefKeywordError(:start_parameters))
    isnothing(target_parameters) && throw(UndefKeywordError(:target_parameters))
    m, n = size(F)
    H = ParameterHomotopy(F, start_parameters, target_parameters)
    f = System(F)
    if is_homogeneous(f)
        vargroups = variable_groups(f)
        if vargroups === nothing
            m ≥ (n - 1) || throw(FiniteException(n - 1 - m))
            H = on_affine_chart(H)
        else
            m ≥ (n - length(vargroups)) || throw(FiniteException(n - length(vargroups) - m))
            H = on_affine_chart(H, length.(vargroups,) .- 1)
        end
    else
        m ≥ n || throw(FiniteException(n - m))
    end

    EndgameTracker(H; tracker_options = tracker_options, options = endgame_options)
end

function start_target_homotopy(
    G::Union{System,AbstractSystem},
    F::Union{System,AbstractSystem};
    γ = 1.0,
    gamma = γ,
)
    f, g = System(F), System(G)

    size(F) == size(G) || error("The provided systems don't have the same size.")
    is_homogeneous(f) == is_homogeneous(g) ||
        error("The provided systems are not both homogeneous.")
    variable_groups(f) == variable_groups(g) ||
        error("The provided systems don't decalare the same variable groups.")

    m, n = size(F)

    H = StraightLineHomotopy(G, F; gamma = gamma)
    if is_homogeneous(f)
        vargroups = variable_groups(f)
        if vargroups === nothing
            m ≥ (n - 1) || throw(FiniteException(n - 1 - m))
            H = on_affine_chart(H)
        else
            m ≥ (n - length(vargroups)) || throw(FiniteException(n - length(vargroups) - m))
            H = on_affine_chart(H, length.(vargroups,) .- 1)
        end
    else
        m ≥ n || throw(FiniteException(n - m))
    end

    H
end

function solver_startsolutions(
    H::Union{Homotopy,AbstractHomotopy},
    starts = nothing;
    seed = nothing,
    kwargs...,
)
    !isnothing(seed) && Random.seed!(seed)
    Solver(EndgameTracker(H), seed), starts
end

"""
    solve(...)

TODO
"""
function solve(
    args...;
    show_progress::Bool = true,
    threading::Bool = Threads.nthreads() > 1,
    catch_interrupt::Bool = true,
    kwargs...,
)
    solver, starts = solver_startsolutions(args...; kwargs...)
    solve(
        solver,
        starts;
        show_progress = show_progress,
        threading = threading,
        catch_interrupt = catch_interrupt,
    )
end

function solve(
    S::Solver,
    starts;
    show_progress::Bool = true,
    threading::Bool = Threads.nthreads() > 1,
    catch_interrupt::Bool = true,
)
    n = length(starts)
    progress = show_progress ? make_progress(n; delay = 0.3) : nothing
    if threading
        threaded_solve(S, starts, progress; catch_interrupt = catch_interrupt)
    else
        serial_solve(S, starts, progress; catch_interrupt = catch_interrupt)
    end
end
(solver::Solver)(starts; kwargs...) = solve(solver, starts; kwargs...)
track(solver::Solver, s; kwargs...) = track(solver.trackers[1], s; kwargs...)

function make_progress(n::Integer; delay::Float64 = 0.0)
    desc = "Tracking $n paths... "
    barlen = min(ProgressMeter.tty_width(desc), 40)
    progress = ProgressMeter.Progress(n; dt = 0.3, desc = desc, barlen = barlen)
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
        ("# singular solutions (real)", "$(stats.singular[]) ($(stats.singular_real[]))"),
        ("# total solutions (real)", "$(nsols[]) ($(nreal[]))"),
    )
end
update_progress!(::Nothing, stats, ntracked) = nothing

function serial_solve(solver::Solver, starts, progress; catch_interrupt::Bool = true)
    path_results = Vector{PathResult}()
    tracker = solver.trackers[1]
    try
        for (k, s) in enumerate(starts)
            r = track(tracker, s; path_number = k)
            push!(path_results, r)
            update!(solver.stats, r)
            update_progress!(progress, solver.stats, k)
        end
    catch e
        (catch_interrupt && isa(e, InterruptException)) || rethrow(e)
    end

    Result(path_results; seed = solver.seed)
end
function threaded_solve(solver::Solver, starts, progress; catch_interrupt::Bool = true)
    S = collect(starts)
    N = length(S)
    path_results = Vector{PathResult}(undef, N)
    interrupted = false
    started = Threads.Atomic{Int}(0)
    finished = Threads.Atomic{Int}(0)
    try
        Threads.resize_nthreads!(solver.trackers)
        tasks = map(solver.trackers) do tracker
            Threads.@spawn begin
                while (k = Threads.atomic_add!(started, 1) + 1) ≤ N && !interrupted
                    r = track(tracker, S[k]; path_number = k)
                    path_results[k] = r
                    nfinished = Threads.atomic_add!(finished, 1) + 1
                    update!(solver.stats, r)
                    update_progress!(progress, solver.stats, nfinished[])
                end
            end
        end
        for task in tasks
            wait(task)
        end
    catch e
        if (
            isa(e, InterruptException) ||
            (isa(e, TaskFailedException) && isa(e.task.exception, InterruptException))
        )
            interrupted = true
        end
        if !interrupted || !catch_interrupt
            rethrow(e)
        end
    end
    # if we got interrupted we need to remove the unassigned filedds
    if interrupted
        assigned_results = Vector{PathResult}()
        for i = 1:started[]
            if isassigned(path_results, i)
                push!(assigned_results, path_results[i])
            end
        end
        Result(assigned_results; seed = solver.seed)
    else
        Result(path_results; seed = solver.seed)
    end

end

"""
    paths_to_track(
        f::Union{System,AbstractSystem};
        start_system::Symbol = :polyhedral,
        kwargs...)

Returns the number of paths tracked when calling [`solve`](@ref) with the given arguments.
"""
function paths_to_track(
    f::Union{System,AbstractSystem};
    start_system::Symbol = :polyhedral,
    kwargs...,
)
    paths_to_track(f, Val(start_system); kwargs...)
end
