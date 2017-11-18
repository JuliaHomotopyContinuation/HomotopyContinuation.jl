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
            if result.returncode == :success
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
        endgame_results = map(r -> EndgamerResult(endgamer, r), tracked_paths)
    end

    # TODO: We can do a second pathcrossing check here:
    # The cauchy endgame will give us the winding number. This gives the multiplicity.
    # So we should have a match between the winding number and the number of solutions
    # at a given point. Otherwise some pathcrossing happend.

    # Refine solution pass
    refined_endgame_results = map(r -> refineresult(r, solver), endgame_results)

    pathresults = broadcast(startvalues, tracked_paths, refined_endgame_results) do a, b, c
        PathResult(a, b, c, solver)
    end
    # Return solution
    Result(pathresults)
end

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

function refineresult(r::EndgamerResult, solver)
    if r.returncode == :success
        @unpack abstol, refinement_maxiters = solver.options
        @unpack H, cfg, cache = solver.pathtracker.low
        # TODO: we should switch to a higher precision if necessary
        # Since we have the winding number available
        # See the cauchy endgame test, the refinement is nearly useless...
        sol = copy(r.solution)
        try
            correct!(sol, 0.0, H, cfg, cache, abstol, refinement_maxiters)
        catch err
            if !isa(err, Base.LinAlg.SingularException)
                throw(err)
            end
        end
        return refined_solution(r, sol)
    else
        return r
    end
end

function residual_estimates(solution, tracker::Pathtracker{Low}) where Low
    @unpack H, cfg = tracker.low
    res = evaluate(H, solution, 0.0, cfg)
    jacobian = Homotopy.jacobian(H, solution, 0.0, cfg, true)
    residual = norm(res)
    newton_residual::Float64 = norm(jacobian \ res)


    condition_number::Float64 = Homotopy.κ(H, solution, 0.0, cfg)

    residual, newton_residual, condition_number
end
