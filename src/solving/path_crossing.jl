using ..Utilities
import ..ProjectiveVectors


"""
    pathcrossing_check!(tracked_paths, solvers, start_solutions)

Check for possble path crossings and correct them if possible. Updates the tracked paths
and returns a tuple `(n, indices)` where `n` is the number of initial crossed paths
and `indices` are the indices which could not be resolved.
"""
function pathcrossing_check!(tracked_paths, solvers, start_solutions)
    t₁, t_endgame, _ = t₁_t_endgame_t₀(solvers)

    # We start with a list of all indices where some crossing happened.
    crossed_path_indices = check_crossed_paths(tracked_paths, cross_tol(solvers))

    ncrossedpaths = length(crossed_path_indices)
    # No paths crossed -> done :)
    if ncrossedpaths == 0
        return (0, Int[])
    end

    # Now we are trying it with tighter precision to avoid the jumping
    original_tol = tol(solvers)
    original_corrector_maxiters = corrector_maxiters(solvers)

    set_tol!(solvers, min(original_tol * 1e-2, 1e-10))

    Parallel.tforeach(solvers, crossed_path_indices) do solver, tid, k
        x₁ = start_solutions[k]
        tracked_paths[k] = trackpath(solver, x₁, t₁, t_endgame)
    end

    crossed_path_indices = check_crossed_paths(tracked_paths, cross_tol(solvers))

    if !isempty(crossed_path_indices)
        # Now we are trying it with less corrector steps and even smaller tol
        set_tol!(solvers, min(original_tol * 1e-3, 1e-12))
        set_corrector_maxiters!(solvers, 1)

        Parallel.tforeach(solvers, crossed_path_indices) do solver, tid, k
            x₁ = start_solutions[k]
            tracked_paths[k] = trackpath(solver, x₁, t₁, t_endgame)
        end
        crossed_path_indices = check_crossed_paths(tracked_paths, cross_tol(solvers))
    end

    # No we reset the options
    set_tol!(solvers, original_tol)
    set_corrector_maxiters!(solvers, original_corrector_maxiters)

    return (ncrossedpaths, crossed_path_indices)
end

cross_tol(solver::Solver) = solver.options.pathcrossing_tol
cross_tol(solvers::Solvers) = cross_tol(solvers[1])

tol(solver::Solver) = PathTracking.tol(solver.tracker)
tol(solvers::Solvers) = tol(solvers[1])
corrector_maxiters(solver::Solver) = PathTracking.corrector_maxiters(solver.tracker)
corrector_maxiters(solvers::Solvers) = corrector_maxiters(solvers[1])

set_tol!(solver::Solver, tol) = PathTracking.set_tol!(solver.tracker, tol)
function set_tol!(solvers::Solvers, tol)
    for solver in solvers
        set_tol!(solver, tol)
    end
end
set_corrector_maxiters!(solver::Solver, n) = PathTracking.set_corrector_maxiters!(solver.tracker, n)
function set_corrector_maxiters!(solvers::Solvers, n)
    for solver in solvers
        set_corrector_maxiters!(solver, n)
    end
end


"""
    check_crossed_paths(paths, tol)

Check whether two solutions are closer than the given tolerance since this
indicates that path crossing happened.
This assumes that the paths were not tracked until t=0.
"""
function check_crossed_paths(paths, tol)
    #TODO: This should be multithreaded... AND is really expensive...5
    crossed_path_indices = Int[]
    path_handled = falses(length(paths))
    for i=1:length(paths)-1
        if path_handled[i]
            continue
        end
        x = paths[i].x
        crossing = false
        for j=i+1:length(paths)
            if !path_handled[j] && ProjectiveVectors.infinity_norm(x, paths[j].x) < tol
                push!(crossed_path_indices, j)
                crossing = true
                path_handled[j] = true
            end
        end
        if crossing
            push!(crossed_path_indices, i)
        end
    end

    crossed_path_indices
end
