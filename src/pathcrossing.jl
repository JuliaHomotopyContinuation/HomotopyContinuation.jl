
function pathcrossing_check!(tracked_paths::Vector{PathtrackerResult{T}}, solver, wp=nothing) where T
    @unpack pathtracker, options = solver
    @unpack endgame_start, pathcrossing_tolerance = options
    in_projective = is_projective_tracker(pathtracker)
     # PATH CROSSING CHECK 1
     # get a list of paths where path crossing happened
     crossed_paths_indices =
        check_crossed_paths(tracked_paths, pathcrossing_tolerance, in_projective)
     user_pathtracker_options = deepcopy(pathtracker.options)

     if !isempty(crossed_paths_indices)
          # We try again with a tighter pathtracking tolerance
         pathtracker.options.abstol = min(pathtracker.options.abstol * 1e-2, 1e-8)
         if isa(wp, CachingPool)
             # we need to clear the cachingpool since we modified options
             clear!(wp)
             startvalues = map(i -> tracked_paths[i].startvalue, crossed_paths_indices)
             results =
                pmap(wp, s -> track_to_endgame_zone(s, solver), startvalues; batch_size=solver.options.batch_size)::Vector{PathtrackerResult{T}}

            for (i, r) in zip(crossed_paths_indices, results)
                tracked_paths[i] = r
            end
        else
            for i in crossed_paths_indices
                tracked_paths[i] = track_to_endgame_zone(tracked_paths[i].startvalue, solver)
            end
        end

         crossed_paths_indices =
            check_crossed_paths(tracked_paths, pathcrossing_tolerance, in_projective)
     end
     if !isempty(crossed_paths_indices) && pathtracker.options.corrector_maxiters > 1
          # We try again with less newton correcotr steps
         pathtracker.options.corrector_maxiters = min(pathtracker.options.corrector_maxiters - 1, 2)
         pathtracker.options.abstol *= 1e-2

         if isa(wp, CachingPool)
             # we need to clear the cachingpool since we modified options
             clear!(wp)
             startvalues = map(i -> tracked_paths[i].startvalue, crossed_paths_indices)
             results =
                pmap(wp, s -> track_to_endgame_zone(s, solver), startvalues; batch_size=solver.options.batch_size)::Vector{PathtrackerResult{T}}

            for (i, r) in zip(crossed_paths_indices, results)
                tracked_paths[i] = r
            end
        else
            for i in crossed_paths_indices
                tracked_paths[i] = track_to_endgame_zone(tracked_paths[i].startvalue, solver)
            end

            crossed_paths_indices =
                check_crossed_paths(tracked_paths, pathcrossing_tolerance, in_projective)
        end
     end
     # TODO: SWITCH TO HIGHER PRECISION IF NECESSARY


     # get the defaults back
     solver.pathtracker.options = user_pathtracker_options

     if isa(wp, CachingPool)
         clear!(wp)
     end

     nothing
end

"""
    check_crossed_paths(paths::Vector{PathtrackerResult{T}} [, tolerance=1e-12])

Split the given paths in two arrays. One for all paths were no path-crossing happend and
one for the paths where we assume that path crossing happnened.
This assumes that the paths were not tracked until t=0.
With probability 1 all solutions should be different. Thus, if two solutions are too similar
(given the passed `tolerance`) we assume hat path crossing happened.
"""
function check_crossed_paths(
    paths::Vector{PathtrackerResult{T}}, tol, in_projective)::Vector{Int} where T
    crossed_path_indices = Int[]
    path_handled = falses(length(paths))

    # we will use the squared norm (its cheaper to compute)
    tolerance = tol * tol
    for i=1:length(paths)-1
        if path_handled[i]
            continue
        end
        r = paths[i]
        x0 = r.solution
        crossing = false
        for j=i+1:length(paths)
            if !path_handled[j]
                if in_projective
                    crossed = projectivenorm2(x0, paths[j].solution) < tolerance
                else
                    crossed = norm(x0 - paths[j].solution)^2 < tolerance
                end
                if crossed
                    push!(crossed_path_indices, j)
                    crossing = true
                    path_handled[j] = true
                end
            end
        end
        if crossing
            push!(crossed_path_indices, i)
        end
    end

    crossed_path_indices
end
