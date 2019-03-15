"""
    pathcrossing_check!(tracked_paths, solvers, start_solutions)

Check for possble path crossings and correct them if possible. Updates the tracked paths
and returns a tuple `(n, indices)` where `n` is the number of initial crossed paths
and `indices` are the indices which could not be resolved.
"""
function pathcrossing_check!(tracked_paths, solvers, start_solutions)
    t₁, t_endgame, _ = t₁_t_endgame_t₀(solvers)
    original_accuracy = accuracy(solvers)
    cross_accuracy = original_accuracy * 10
    # We start with a list of all indices where some crossing happened.
    crossed_path_indices = check_crossed_paths(tracked_paths, cross_accuracy)
    ncrossedpaths = length(crossed_path_indices)
    # No paths crossed -> done :)
    if ncrossedpaths == 0
        return (0, Int[])
    end

    # Now we are trying it with tighter precision to avoid the jumping
    original_accuracy = accuracy(solvers)
    original_corrector_maxiters = max_corrector_iters(solvers)

    set_accuracy!(solvers, min(original_accuracy * 1e-2, 1e-10))

    tforeach(solvers, crossed_path_indices) do solver, tid, k
        x₁ = start_solutions[k]
        tracked_paths[k] = trackpath(solver, x₁, t₁, t_endgame)
    end

    crossed_path_indices = check_crossed_paths(tracked_paths, cross_accuracy)

    if !isempty(crossed_path_indices)
        # Now we are trying it with less corrector steps and even smaller accuracy
        set_accuracy!(solvers, min(original_accuracy * 1e-3, 1e-12))
        set_max_corrector_iters!(solvers, 1)

        tforeach(solvers, crossed_path_indices) do solver, tid, k
            x₁ = start_solutions[k]
            tracked_paths[k] = trackpath(solver, x₁, t₁, t_endgame)
        end
        crossed_path_indices = check_crossed_paths(tracked_paths, cross_accuracy)
    end

    # No we reset the options
    set_accuracy!(solvers, original_accuracy)
    set_max_corrector_iters!(solvers, original_corrector_maxiters)

    return (ncrossedpaths, crossed_path_indices)
end

accuracy(solver::Solver) = accuracy(solver.tracker)
accuracy(solvers::Solvers) = accuracy(solvers[1])
max_corrector_iters(solver::Solver) = max_corrector_iters(solver.tracker)
max_corrector_iters(solvers::Solvers) = max_corrector_iters(solvers[1])

set_accuracy!(solver::Solver, accuracy) = set_accuracy!(solver.tracker, accuracy)
function set_accuracy!(solvers::Solvers, accuracy)
    for solver in solvers
        set_accuracy!(solver, accuracy)
    end
end
set_max_corrector_iters!(solver::Solver, n) = set_max_corrector_iters!(solver.tracker, n)
function set_max_corrector_iters!(solvers::Solvers, n)
    for solver in solvers
        set_max_corrector_iters!(solver, n)
    end
end


"""
    check_crossed_paths(paths, tol)

Check whether two solutions are closer than the given tolerance since this
indicates that path crossing happened.
This assumes that the paths were not tracked until t=0.
"""
function check_crossed_paths(paths::Vector{PR}, tol) where {PR<:PathTrackerResult}
    V = map(p -> LinearAlgebra.normalize(p.x.data), paths)
    vcat(multiplicities(V, fubini_study, tol = tol)...)
end
