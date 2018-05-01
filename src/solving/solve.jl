include("path_crossing.jl")

solve(solver, start_solutions) = solve(solver, collect(start_solutions))

function solve(solvers, start_solutions::AbstractVector)
    nblas_threads = single_thread_blas()

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

function track_to_endgamezone(solvers, start_solutions)
    t₁, t_endgame, _ = t₁_t_endgame_t₀(solvers)
    Parallel.tmap(solvers, start_solutions) do solver, tid, x₁
        trackpath(solver, x₁, t₁, t_endgame)
    end
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
    ncrossedpaths, crossing_indices =
        pathcrossing_check!(endgame_zone_results, solvers, start_solutions)

    Parallel.tmap(solvers, 1:n) do solver, tid, k
        x₁, r = start_solutions[k], endgame_zone_results[k]
        if r.returncode == :success
            result = Endgame.play(solver.endgamer, r.x, t_endgame)
            PathResult(solver.prob, x₁, t₀, result, solver.cache.pathresult)
        else
            PathResult(solver.prob, x₁, t₀, r, solver.cache.pathresult)
        end
    end
end
