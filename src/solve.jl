export solve

"""
    solve(H::AbstractHomotopy, startvalues_s, [algorithm]; kwargs...)

Solve the homotopy `H` via homotopy continuation with the given `startvalues_s` and the given
`algorithm`.


    solve(f::Vector{<:MP.AbstractPolynomial{T}}; homotopytype=GeodesicOnTheSphere, algorithm=SphericalPredictorCorrector(), kwargs...)
    solve(f::MP.AbstractPolynomial{T}; homotopytype=GeodesicOnTheSphere, algorithm=SphericalPredictorCorrector(), kwargs...)

Solves the polynomial system `f` via homotopy continuation. This uses a totaldegree homotopy
of the given `homotopytype` and uses the given `algorithm` for pathtracking.
"""
function solve end

solve(a, b::Vector{<:Number}, c, d, e; kwargs...) = solve(a, [b], c, d, e; kwargs...)
solve(a, b::Vector{<:Number}, c, d; kwargs...) = solve(a, [b], c, d; kwargs...)
solve(a, b::Vector{<:Number}, c; kwargs...) = solve(a, [b], c; kwargs...)
solve(a, b::Vector{<:Number}; kwargs...) = solve(a, [b]; kwargs...)

solve(f::MP.AbstractPolynomial; kwargs...) = solve([f]; kwargs...)
solve(f::MP.AbstractPolynomial, pa::APA; kwargs...) = solve([f], pa; kwargs...)
solve(f::MP.AbstractPolynomial, pa::APA, ea::AEA; kwargs...) = solve([f], pa, ea; kwargs...)
solve(f::MP.AbstractPolynomial, HT; kwargs...) = solve([f], HT; kwargs...)
solve(f::MP.AbstractPolynomial, HT, pa::APA; kwargs...) = solve([f], HT, pa; kwargs...)
solve(f::MP.AbstractPolynomial, HT, pa::APA, ea::AEA; kwargs...) = solve([f], HT, pa, ea; kwargs...)

function solve(F::Vector{<:MP.AbstractPolynomial{T}}; kwargs...) where T
      H, s = totaldegree(StraightLineHomotopy{promote_type(Complex128, T)}, F)
      solve(H, s; kwargs...)
end
function solve(F::Vector{<:MP.AbstractPolynomial{T}}, pa::APA; kwargs...) where T
      H, s = totaldegree(StraightLineHomotopy{promote_type(Complex128, T)}, F)
      solve(H, s, pa; kwargs...)
end
function solve(F::Vector{<:MP.AbstractPolynomial{T}}, pa::APA, ea::AEA; kwargs...) where T
      H, s = totaldegree(StraightLineHomotopy{promote_type(Complex128, T)}, F)
      solve(H, s, pa, ea; kwargs...)
end
function solve(F::Vector{<:MP.AbstractPolynomial{T}}, ::Type{AHT}; kwargs...) where {T, AHT<:AH}
      H, s = totaldegree(promote_type(AHT, promote_type(Complex128, T)), F)
      solve(H, s; kwargs...)
end
function solve(F::Vector{<:MP.AbstractPolynomial{T}}, ::Type{AHT}, pa::APA; kwargs...) where {T, AHT<:AH}
      H, s = totaldegree(promote_type(AHT, promote_type(Complex128, T)), F)
      solve(H, s, pa; kwargs...)
end
function solve(F::Vector{<:MP.AbstractPolynomial{T}}, ::Type{AHT}, pa::APA, ea::AEA; kwargs...) where {T, AHT<:AH}
      H, s = totaldegree(promote_type(AHT, promote_type(Complex128, T)), F)
      solve(H, s, pa, ea; kwargs...)
end

solve(H::AH, s; kwargs...) = solve(Solver(H; kwargs...), s)
solve(H::AH, s, pa::APA; kwargs...) = solve(Solver(H, pa; kwargs...), s)
solve(H::AH, s, pa::APA, ea::AEA; kwargs...) = solve(Solver(H, pa, ea; kwargs...), s)
solve(H::AH, s, pa::APA, ea::AEA, HPT::Type{<:AbstractFloat}; kwargs...) = solve(Solver(H, pa, ea, HPT; kwargs...), s)


# This currently relies on the fact that we can keep all solutions in memory. This could
# not be viable for big systems...
# The main problem for a pure streaming / iterator solution is the pathcrossing check
function solve(solver::Solver{AH}, startvalues) where {T, AH<:AbstractHomotopy{T}}
    @unpack options, pathtracker, endgamer = solver
    @unpack endgame_start, pathcrossing_check = options

    tracked_paths::Vector{PathtrackerResult{T}} = map(startvalues) do startvalue
        track!(pathtracker, startvalue, 1.0, endgame_start)
        PathtrackerResult(pathtracker, false)
    end
    if pathcrossing_check
        pathcrossing_check!(tracked_paths, solver)
    end

    endgame_results = Vector{EndgamerResult{T}}()
    if endgame_start > 0.0
        endgame_results = map(tracked_paths) do result
            if result.retcode == :success
                endgame!(endgamer, result.solution, endgame_start)
                EndgamerResult(endgamer)
            else
                # We just carry over the result from the failed path.
                # This should not happen since we catch things above, but you never know...
                EndgamerResult(endgamer, result)
            end
        end
    else
        # we just carry over the results to make the rest of the code clearer
        endgame_results = pmap(r -> EndgamerResult(endgamer, r), tracked_paths)
    end

    # TODO: We can do a second pathcrossing check here:
    # The cauchy endgame will give us the winding number. This gives the multiplicity.
    # So we should have a match between the winding number and the number of solutions
    # at a given point. Otherwise some pathcrossing happend.

    # Refine solution pass

    results::Vector{PathResult{T}} = map(startvalues, tracked_paths, endgame_results) do s, esr, er
        refine_and_pathresult(s, esr, er, pathtracker, options.abstol, options.refinement_maxiters)
    end

    # Return solution
    Result(results)
end

function pathcrossing_check!(tracked_paths::Vector{PathtrackerResult{T}}, solver) where T
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

         @show pathtracker.options
         for i in crossed_paths_indices
              track!(pathtracker, tracked_paths[i].startvalue, 1.0, endgame_start)
              tracked_paths[i] = PathtrackerResult(pathtracker, false)
         end

         crossed_paths_indices =
            check_crossed_paths(tracked_paths,
                pathcrossing_tolerance, in_projective)
     end
     if !isempty(crossed_paths_indices) && pathtracker.options.corrector_maxiters > 1
          # We try again with less newton correcotr steps
         pathtracker.options.corrector_maxiters = min(pathtracker.options.corrector_maxiters - 1, 2)

         @show pathtracker.options
         for i in crossed_paths_indices
              track!(pathtracker, tracked_paths[i].startvalue, 1.0, endgame_start)
              tracked_paths[i] = PathtrackerResult(pathtracker, false)
         end

         crossed_paths_indices =
            check_crossed_paths(tracked_paths,
                pathcrossing_tolerance, in_projective)
     end
     # TODO: SWITCH TO HIGHER PRECISION IF NECESSARY


     # get the defaults back
     solver.pathtracker.options = user_pathtracker_options

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

function refine_and_pathresult(
    startvalue,
    endgame_start_result::PathtrackerResult{T},
    endgamer_result::EndgamerResult,
    pathtracker,
    abstol,
    refinement_maxiters) where T
    @unpack retcode, solution, windingnumber = endgamer_result

    # we refine the solution if possible
    if retcode == :success
        solution = refinesolution(solution, pathtracker, windingnumber, abstol, refinement_maxiters)
    end

    residual, newton_residual, condition_jacobian = residual_estimates(solution, pathtracker)

    # check whether startvalue was affine and our solution is projective
    N = length(startvalue)
    if length(solution) == N + 1
        # make affine
        # This is a more memory efficient variant from:
        # solution = solution[2:end] / solution[1]
        homog_var = solution[1]
        for i=2:N+1
            solution[i - 1] = solution[i] / homog_var
        end
        resize!(solution, N)

        homogenous_coordinate_magnitude = norm(homog_var)
    else
        homogenous_coordinate_magnitude = 1.0
    end

    PathResult{T}(
        retcode,
        solution,
        residual,
        newton_residual,
        condition_jacobian,
        windingnumber,
        homogenous_coordinate_magnitude,
        copy(startvalue),
        endgame_start_result.iterations,
        endgamer_result.iterations,
        endgamer_result.npredictions
        )
end


function refinesolution(solution, tracker::Pathtracker, windingnumber, abstol, maxiters)
    @unpack H, cfg, cache = tracker.low
    # TODO: we should switch to a higher precision if necessary
    # Since we have the winding number available
    # See the cauchy endgame test, the refinement is nearly useless...
    sol = copy(solution)
    try
        correct!(sol, 0.0, H, cfg, cache, abstol, maxiters)
    catch err
        if !isa(err, Base.LinAlg.SingularException)
            throw(err)
        end
    end
    sol
end

function residual_estimates(solution, tracker::Pathtracker{Low}) where Low
    @unpack H, cfg = tracker.low
    res = evaluate(H, solution, 0.0, cfg)
    jacobian = Homotopy.jacobian(H, solution, 0.0, cfg, true)
    residual = norm(res)
    newton_residual::Float64 = norm(jacobian \ res)
    condition_jacobian::Float64 = cond(jacobian)

    residual, newton_residual, condition_jacobian
end
