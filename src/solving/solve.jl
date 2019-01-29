const Solvers = Vector{<:Solver}

include("path_crossing.jl")

internal_solve(solver, start_solutions) = internal_solve(solver, collect(start_solutions))

function internal_solve(solvers, start_solutions::AbstractVector)
    nblas_threads = single_thread_blas()

    endgame_zone_results = track_to_endgamezone(solvers, start_solutions)
    results = endgame(solvers, start_solutions, endgame_zone_results)

    set_num_BLAS_threads(nblas_threads)

    if all(r -> r.solution_type == :affine, results)
        AffineResult(results, seed(solvers))
    elseif all(r -> r.solution_type == :projective, results)
        ProjectiveResult(results, seed(solvers))
    else
        @warn("Something went wrong. There are both affine and projective solutions.")
        ProjectiveResult(results, seed(solvers))
    end
end

"""
    single_thread_blas()

We set the number of BLAS threads to 1 since we multithread by ourself.
But even if we would not, the overhead of the threading is so large
that the single threaded version is around two times faster (at least on Mac Julia v0.6.2).
Returns the previous number of BLAS threads.
"""
function single_thread_blas()
    nblas_threads = get_num_BLAS_threads()
    set_num_BLAS_threads(1)
    nblas_threads
end

seed(solver::Solver) = solver.seed
seed(solvers::Solvers) = seed(solvers[1])
t₁_t_endgame_t₀(solver::Solver) = solver.t₁, solver.options.endgame_start, solver.t₀
t₁_t_endgame_t₀(solvers::Solvers) = t₁_t_endgame_t₀(solvers[1])
report_progress(solver::Solver) = solver.options.report_progress
report_progress(solvers::Solvers) = report_progress(solvers[1])

function track_to_endgamezone(solvers, start_solutions)
    t₁, t_endgame, _ = t₁_t_endgame_t₀(solvers)
    n = length(start_solutions)

    if report_progress(solvers)
        p = ProgressMeter.Progress(n, 0.5, "Tracking $(length(start_solutions)) paths to endgame zone...") # minimum update interval: 1 second

        result = tmap(solvers, 1:n) do solver, tid, k
            if tid == 1
                ProgressMeter.update!(p, max(1, k-1), showvalues=((:tracked, k - 1),))
            end
            trackpath(solver, start_solutions[k], t₁, t_endgame)
        end

        ProgressMeter.update!(p, n, showvalues=((:tracked, n),))
    else
        result = tmap(solvers, 1:n) do solver, tid, k
            trackpath(solver, start_solutions[k], t₁, t_endgame)
        end
    end

    result
end

trackpath(solver::Solver, x₁, t₁, t₀) = track(solver.tracker, embed(solver.prob, x₁), t₁, t₀)

function endgame(solvers, start_solutions, endgame_zone_results)
    _, t_endgame, t₀ = t₁_t_endgame_t₀(solvers)
    n = length(start_solutions)

    if t₀ == t_endgame
        return tmap(solvers, 1:n) do solver, tid, k
            x₁, r = start_solutions[k], endgame_zone_results[k]
            PathResult(solver.prob, k, x₁, r.x, t₀, r, solver.cache.pathresult)
        end
    end

    # report_progress(solvers) && println("Checking for crossed paths...")
    ncrossedpaths, crossing_indices = pathcrossing_check!(endgame_zone_results, solvers, start_solutions)
    if report_progress(solvers)
        p = ProgressMeter.Progress(n, 0.5, "Running endgame for $(length(endgame_zone_results)) paths...")

        result = tmap(solvers, 1:n) do solver, tid, k
            if tid == 1
                ProgressMeter.update!(p, max(1, k-1), showvalues=((:completed, k-1),))
            end
            runendgame(solver, tid, k, start_solutions, endgame_zone_results)
        end
        ProgressMeter.update!(p, n, showvalues=((:completed, n),))
    else
        result = tmap(solvers, 1:n) do solver, tid, k
            runendgame(solver, tid, k, start_solutions, endgame_zone_results)
        end
    end

    result
end

function runendgame(solver, tid, k, start_solutions, endgame_zone_results)
    t₁, t_endgame, t₀ = t₁_t_endgame_t₀(solver)
    x₁, r = start_solutions[k], endgame_zone_results[k]
    if r.returncode == PathTrackerStatus.success
        # Run endgame
        result = runendgame(solver.endgame, r.x, t_endgame)
        # If the tracker failed we are probably to late with the endgame.
        if result.returncode == :tracker_failed
            # Rerun with something more away
            new_t = 0.3*(t₁ - t_endgame)
            pr = trackpath(solver::Solver, x₁, t₁, new_t)
            if pr.returncode == PathTrackerStatus.success
                result = runendgame(solver.endgame, pr.x, new_t)
            end
        end
        return PathResult(solver.prob, k, x₁, r.x, t₀, result, solver.cache.pathresult)
    else
        # If we even didn't come to the endgame zone we start earlier.
        result = runendgame(solver.endgame, embed(solver.prob, x₁), 1.0)
        if result.returncode == :success || result.returncode == :at_infinity
            return PathResult(solver.prob, k, x₁, r.x, t₀, result, solver.cache.pathresult)
        else
            return PathResult(solver.prob, k, x₁, r.x, t₀, r, solver.cache.pathresult)
        end
    end
end
