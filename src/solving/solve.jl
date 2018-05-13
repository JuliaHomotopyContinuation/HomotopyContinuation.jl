include("path_crossing.jl")

solve(solver, start_solutions) = solve(solver, collect(start_solutions))

function solve(solvers, start_solutions::AbstractVector)
    nblas_threads = single_thread_blas()

    if report_progress(solvers)
        println("$(length(start_solutions)) paths to track")
    end
    endgame_zone_results = track_to_endgamezone(solvers, start_solutions)
    results = endgame(solvers, start_solutions, endgame_zone_results)

    BLAS.set_num_threads(nblas_threads)

    results
end

"""
    single_thread_blas()

We set the number of BLAS threads to 1 since we multithread by ourself.
But even if we would not, the overhead of the threading is so large
that the single threaded version is around two times faster (at least on Mac Julia v0.6.2).
Returns the previous number of BLAS threads.
"""
function single_thread_blas()
    nblas_threads = Utilities.get_num_BLAS_threads()
    BLAS.set_num_threads(1)
    nblas_threads
end

t₁_t_endgame_t₀(solver::Solver) = solver.t₁, solver.options.endgame_start, solver.t₀
t₁_t_endgame_t₀(solvers::Solvers) = t₁_t_endgame_t₀(solvers[1])
report_progress(solver::Solver) = solver.options.report_progress
report_progress(solvers::Solvers) = report_progress(solvers[1])

function track_to_endgamezone(solvers, start_solutions)
    t₁, t_endgame, _ = t₁_t_endgame_t₀(solvers)
    n = length(start_solutions)

    if report_progress(solvers)
        p = ProgressMeter.Progress(n, 0.5, "Tracking paths to endgame zone...") # minimum update interval: 1 second

        result = Parallel.tmap(solvers, 1:n) do solver, tid, k
            if tid == 1
                ProgressMeter.update!(p, max(1, k-1), showvalues=((:tracked, k - 1),))
            end
            trackpath(solver, start_solutions[k], t₁, t_endgame)
        end

        ProgressMeter.update!(p, n, showvalues=((:tracked, n),))
    else
        result = Parallel.tmap(solvers, 1:n) do solver, tid, k
            trackpath(solver, start_solutions[k], t₁, t_endgame)
        end
    end

    result
end

@inline trackpath(solver::Solver, x₁, t₁, t₀) = PathTracking.track(solver.tracker, Problems.embed(solver.prob, x₁), t₁, t₀)

function endgame(solvers, start_solutions, endgame_zone_results)
    _, t_endgame, t₀ = t₁_t_endgame_t₀(solvers)
    n = length(start_solutions)

    if t₀ == t_endgame
        return Parallel.tmap(solvers, 1:n) do solver, tid, k
            x₁, r = start_solutions[k], endgame_zone_results[k]
            PathResult(solver.prob, x₁, t₀, r, solver.cache.pathresult)
        end
    end

    report_progress(solvers) && println("Checking for crossed paths...")
    ncrossedpaths, crossing_indices = pathcrossing_check!(endgame_zone_results, solvers, start_solutions)
    if report_progress(solvers)
        p = ProgressMeter.Progress(n, 0.5, "Running endgame...")

        result = Parallel.tmap(solvers, 1:n) do solver, tid, k
            if tid == 1
                ProgressMeter.update!(p, max(1, k-1), showvalues=((:tracked, k-1),))
            end
            runendgame(solver, tid, k, start_solutions, endgame_zone_results, t_endgame, t₀)
        end
        ProgressMeter.update!(p, n, showvalues=((:tracked, n),))
    else
        result = Parallel.tmap(solvers, 1:n) do solver, tid, k
            runendgame(solver, tid, k, start_solutions, endgame_zone_results, t_endgame, t₀)
        end
    end

    result
end

@inline function runendgame(solver, tid, k, start_solutions, endgame_zone_results, t_endgame, t₀)
    x₁, r = start_solutions[k], endgame_zone_results[k]
    if r.returncode == :success
        result = Endgame.play(solver.endgamer, r.x, t_endgame)
        return PathResult(solver.prob, x₁, t₀, result, solver.cache.pathresult)
    else
        return PathResult(solver.prob, x₁, t₀, r, solver.cache.pathresult)
    end
end
