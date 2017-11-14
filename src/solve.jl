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

# This are just convenience wrappers to pass arguments to Solver.
# Solver currently takes up to 5 positional arguments.
solve(a, b, c, d, e; kwargs...) = solve(Solver(a, b, c, d, e; kwargs...))
solve(a, b, c, d; kwargs...) = solve(Solver(a, b, c, d; kwargs...))
solve(a, b, c; kwargs...) = solve(Solver(a, b, c; kwargs...))
solve(a, b; kwargs...) = solve(Solver(a, b; kwargs...))
solve(a; kwargs...) = solve(Solver(a; kwargs...))

solve(a, b::Vector{<:Number}, c, d, e; kwargs...) = solve(Solver(a, [b], c, d, e; kwargs...))
solve(a, b::Vector{<:Number}, c, d; kwargs...) = solve(Solver(a, [b], c, d; kwargs...))
solve(a, b::Vector{<:Number}, c; kwargs...) = solve(Solver(a, [b], c; kwargs...))
solve(a, b::Vector{<:Number}; kwargs...) = solve(Solver(a, [b]; kwargs...))


# This currently relies on the fact that we can keep all solutions in memory. This could
# not be viable for big systems...
# The main problem for a pure streaming / iterator solution is the pathcrossing check
function solve(solver::Solver{AH}) where {T, AH<:AbstractHomotopy{T}}
    @unpack options, pathtracker, endgamer, startvalues = solver
    @unpack endgame_start = options

    # TODO: pmap. Will this preserve the order of the arguments? Otherwise we have to
    # return a tuple or something like that
    endgame_start_results::Vector{PathtrackerResult{T}} = pmap(startvalues) do startvalue
        track!(pathtracker, startvalue, 1.0, endgame_start)
        # do we need informations about  condition_jacobian?
        PathtrackerResult(pathtracker, false)
    end


    # TODO: Rerun failed paths with higher precision.

    # TODO: We can add a first pathcrossing check here:
    # Check whether two endgame_start_results solutions are "close".
    # Since the paths should all be unique this should not happen
    # If we find two of them we should rerun them with tighter bounds

    # TODO: pmap. Will this preserve the order of the arguments? Otherwise we have to
    # return a tuple or something like that
    endgame_results = Vector{EndgamerResult{T}}()
    if endgame_start > 0.0
        endgame_results = pmap(endgame_start_results) do result
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
        endgame_results = pmap(endgame_results, r -> EndgamerResult(endgamer, r), endgame_start_results)
    end

    # TODO: We can do a second pathcrossing check here:
    # The cauchy endgame will give us the winding number. This gives the multiplicity.
    # So we should have a match between the winding number and the number of solutions
    # at a given point. Otherwise some pathcrossing happend.

    # Refine solution pass

    results::Vector{PathResult{T}} = map(startvalues, endgame_start_results, endgame_results) do s, esr, er
        refine_and_pathresult(s, esr, er, pathtracker, options.abstol,  options.tol, options.refinement_maxiters)
    end

    # Return solution
    Result(results)
end


function refine_and_pathresult(
    startvalue,
    endgame_start_result::PathtrackerResult{T},
    endgamer_result::EndgamerResult,
    pathtracker,
    abstol,
    tol,
    refinement_maxiters) where T
    @unpack returncode, solution, windingnumber = endgamer_result

    # we refine the solution if possible
    if returncode == :success
        solution = refinesolution(solution, pathtracker, windingnumber, abstol, refinement_maxiters)
    end

    residual, newton_residual, condition_number = residual_estimates(solution, pathtracker)

    # check whether startvalue was affine and our solution is projective
    N = length(startvalue)
    if length(solution) == N + 1
        # make affine

        homog_var = solution[1]
        solution = solution[2:end]
        angle_to_infinity = atan(abs(homog_var)/norm(solution))
        if angle_to_infinity < tol
            returncode = :at_infinity
        end

        scale!(solution, inv(homog_var))

    else
        angle_to_infinity = NaN
    end

    if windingnumber > 1 || condition_number > inv(tol)
        returncode = :Singular
    end

    PathResult{T}(
        returncode,
        solution,
        residual,
        newton_residual,
        condition_number,
        windingnumber,
        angle_to_infinity,
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
    correct!(sol, 0.0, H, cfg, cache, abstol, maxiters)
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
