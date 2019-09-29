export Solver, solver_startsolutions, solve!

const solve_supported_keywords = [
    :threading,
    :show_progress,
    :path_result_details,
    :save_all_paths,
    :path_jumping_check,
            # deprecated
    :report_progress,
]


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
                results[k] = results[j]
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
                end
            else
                check.solution_mapping[k] = length(check.checkpoint)
                break
            end
        end
    else
        check.solution_mapping[k] = length(check.checkpoint)
    end
    if save_all_paths || is_success(return_code) || is_invalid_startvalue(return_code)
        results[path_number] = PathResult(
            tracker,
            S[path_number],
            path_number;
            details = path_result_details,
        )
    end
    return_code, path_number
end


############
## SOLVER ##
############

struct Solver{PT,UP<:UniquePoints}
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

function Solver(tracker::PathTracker, checkpoint::UniquePoints)
    Solver(tracker, SolveStats(), checkpoint)
end

function solver_startsolutions(args...; kwargs...)
    supported, rest = splitkwargs(kwargs, problem_startsolutions_supported_keywords)
    prob, start_solutions = problem_startsolutions(args...; supported...)
    Solver(prob, start_solutions; rest...), start_solutions
end

accuracy(T::PathTracker) = T.options.min_accuracy
# accuracy(T::PolyhedralTracker) = accuracy(T.generic_tracker)

function construct_tracker(prob::Problem, startsolutions; kwargs...)
    PathTracker(prob, start_solution_sample(startsolutions); kwargs...)
end

############
## solve! ##
############

function solve!(solver::Solver, start_solutions; show_progress::Bool = true, kwargs...)
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
    solve!(solver, start_solutions, progress; kwargs...)
end

function solve!(
    solver::Solver,
    start_solutions,
    progress::Union{Nothing,ProgressMeter.Progress};
    path_result_details::Symbol = :default,
    save_all_paths::Bool = false,
    path_jumping_check::Bool = true,
)
    @unpack trackers, stats = solver
    tracker = trackers

    S = collect_startsolutions(start_solutions)
    n = length(S)

    init!(stats)
    init!(path_jumping_check, n)

    results = Vector{Union{Nothing,result_type(tracker)}}(undef, n)
    results .= nothing

    ntracked = 0
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
        is_success(return_code) && update!(stats, results[path_number])
        ntracked = k
        k % 32 == 0 && update_progress!(progress, k, stats)
    end
    # don't print if it already got printed above
    n % 32 != 0 && update_progress!(progress, n, stats)

    Result(
        remove_nothings(results),
        n,
        seed(tracker);
        multiplicity_tol = 10 * accuracy(tracker),
    )
end

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
