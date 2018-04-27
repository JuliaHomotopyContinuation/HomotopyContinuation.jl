using ..Utilities
import ..ProjectiveVectors

"""
    pathcrossing_check!(tracked_paths, solver)

Check for possble path crossings and correct them if possible. Updates the tracked paths
and returns a tuple `(n, indices)` where `n` is the number of initial crossed paths
and `indices` are the indices which could not be resolved.
"""
function pathcrossing_check!(tracked_paths, solver, start_solutions)
    t₁ = solver.t₁
    t₀ = solver.options.endgame_start
    tracker = solver.tracker

    cross_tol = solver.options.pathcrossing_tol

    # We start with a list of all indices where some crossing happened.
    crossed_paths_indices = check_crossed_paths(tracked_paths, cross_tol)

    ncrossedpaths = length(crossed_paths_indices)
    # No paths crossed -> done :)
    if ncrossedpaths == 0
        return (0, [])
    end

    # Now we are trying it with tighter precision to avoid the jumping
    original_tol = PathTracking.tol(tracker)
    original_corrector_maxiters = PathTracking.corrector_maxiters(tracker)

    PathTracking.set_tol!(tracker, min(original_tol * 1e-2, 1e-8))
    dmap_indices!(tracked_paths, solver.options, crossed_paths_indices, start_solutions) do x₀
        trackpath(solver, x₀, t₁, t₀)
    end

    crossed_paths_indices = check_crossed_paths(tracked_paths, cross_tol)

    if isempty(crossed_paths_indices)
        PathTracking.set_tol!(tracker, original_tol)

        return (ncrossedpaths, [])
    end

    # Now we are trying it with less corrector steps
    PathTracking.set_tol!(tracker, 1e-13)
    PathTracking.set_corrector_maxiters!(tracker, 1)

    dmap_indices!(tracked_paths, solver.options, crossed_paths_indices, start_solutions) do x₀
        trackpath(solver, x₀, t₁, t₀)
    end

    # No we reset the options
    PathTracking.set_tol!(tracker, original_tol)
    PathTracking.set_corrector_maxiters!(tracker, original_corrector_maxiters)

    crossed_paths_indices = check_crossed_paths(tracked_paths, cross_tol)
    @show crossed_paths_indices

    return (ncrossedpaths, crossed_paths_indices)
end

"""
    check_crossed_paths(paths, tol)

Check whether two solutions are closer than the given tolerance since this
indicates that path crossing happened.
This assumes that the paths were not tracked until t=0.
"""
function check_crossed_paths(paths, tol)
    #TODO: This should be multithreaded...
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
