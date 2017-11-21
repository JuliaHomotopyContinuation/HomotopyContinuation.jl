export solve

"""
```
    solve(H::AbstractHomotopy, startvalues_s, [algorithm]; kwargs...)
```

Solve the homotopy `H` via homotopy continuation with the given `startvalues_s` and the given
`algorithm`.

```
    solve(f::Vector{<:MP.AbstractPolynomial{T}})
```

Solves the polynomial system `f` via homotopy continuation. This uses a totaldegree homotopy of type StraightLineHomotopy and the SphericalPredictorCorrector pathtracking routine. To specify homotopy and pathtracker, use

    solve(f::Vector{<:MP.AbstractPolynomial{T}}, [homotopy], [algorithm])

Default is homotopy = StraightLineHomotopy and algorithm = SphericalPredictorCorrector. For instance, use

```
solve(f, StraightLineHomotopy, AffinePredictorCorrector())
```

for solving ``f`` with a GeodesicOnTheSphere homotopy and the AffinePredictorCorrector pathtracking routine.
"""
function solve end


function solve(solver::Solver{AH}, startvalues) where {T, AH<:AbstractHomotopy{T}}
    @unpack options, pathtracker, endgamer = solver
    @unpack endgame_start, pathcrossing_check, parallel_type, batch_size = options

    wp = CachingPool(workers())
    # pmap(wp, i -> setup_workers(i, solver), workers())

    if parallel_type == :pmap
        tracked_paths = pmap(wp,
            s -> track_to_endgame_zone(s, solver), startvalues,
            batch_size=batch_size):: Vector{PathtrackerResult{T}}
    else
        tracked_paths = map(
            s -> track_to_endgame_zone(s, solver), startvalues
            )::Vector{PathtrackerResult{T}}
    end

    if pathcrossing_check
        if parallel_type == :pmap
            pathcrossing_check!(tracked_paths, solver, wp)
        else
            pathcrossing_check!(tracked_paths, solver)
        end
    end

    if endgame_start > 0.0
        if parallel_type == :pmap
            endgame_results = pmap(wp,
                r -> track_endgame(r, solver), tracked_paths,
                batch_size=batch_size)::Vector{EndgamerResult{T}}
        else
            endgame_results = map(
                r -> track_endgame(r, solver), tracked_paths
                )::Vector{EndgamerResult{T}}
        end
    else
        # we just carry over the results to make the rest of the code clearer
        endgame_results = map(r -> EndgamerResult(endgamer, r), tracked_paths)::Vector{EndgamerResult{T}}
    end

    # TODO: We can do a second pathcrossing check here:
    # The cauchy endgame will give us the winding number. This gives the multiplicity.
    # So we should have a match between the winding number and the number of solutions
    # at a given point. Otherwise some pathcrossing happend.

    # Refine solution pass
    if parallel_type == :pmap
        refined_endgame_results = pmap(wp,
            r -> refineresult(r, solver), endgame_results,
            batch_size=batch_size)::Vector{EndgamerResult{T}}
    else
        refined_endgame_results = map(
            r -> refineresult(r, solver), endgame_results
            )::Vector{EndgamerResult{T}}
    end

    pathresults::Vector{PathResult{T}} = broadcast(startvalues, tracked_paths, refined_endgame_results) do a, b, c
        PathResult(a, b, c, solver)
    end
    # # Return solution
    Result(pathresults)
end

"""
    setup_workers(s::Solver)

Ensure that everything is properly setup on each worker. Currently `@view`
in SphericalCache makes some problems. Maybe we can get rid of this at a later
point (e.g. if `@view` does not allocate anymore).
"""
function setup_workers(s::Solver)
    setup_workers(s.pathtracker.low.cache)
    setup_workers(s.pathtracker.high.cache)
end

function track_to_endgame_zone(startvalue::AbstractVector{<:Number}, solver)
      setup_workers(solver)
      track!(solver.pathtracker, startvalue, 1.0, solver.options.endgame_start)
      PathtrackerResult(solver.pathtracker, false)
end
track_to_endgame_zone(svs::AbstractVector, solver) = map(r -> track_to_endgame_zone(r, solver), sv)


function track_endgame(tracked_path_result, solver)
    setup_workers(solver)
    if tracked_path_result.returncode == :success
        endgame!(solver.endgamer, tracked_path_result.solution, solver.options.endgame_start)
        EndgamerResult(solver.endgamer)
    else
        # We just carry over the result from the failed path.
        # This should not happen since we catch things above, but you never know...
        EndgamerResult(solver.endgamer, tracked_path_result)
    end
end
track_endgame(rs::AbstractVector, solver) = map(r -> track_endgame(r, solver), rs)

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

solve(H::AH, s::Vector{<:Number}; kwargs...) = solve(Solver(H; kwargs...), [s])
solve(H::AH, s::Vector{<:Number}, pa::APA; kwargs...) = solve(Solver(H, pa; kwargs...), [s])
solve(H::AH, s::Vector{<:Number}, pa::APA, ea::AEA; kwargs...) = solve(Solver(H, pa, ea; kwargs...), [s])
solve(H::AH, s::Vector{<:Number}, pa::APA, ea::AEA, HPT::Type{<:AbstractFloat}; kwargs...) = solve(Solver(H, pa, ea, HPT; kwargs...), [s])

solve(H::AH, s; kwargs...) = solve(Solver(H; kwargs...), s)
solve(H::AH, s, pa::APA; kwargs...) = solve(Solver(H, pa; kwargs...), s)
solve(H::AH, s, pa::APA, ea::AEA; kwargs...) = solve(Solver(H, pa, ea; kwargs...), s)
solve(H::AH, s, pa::APA, ea::AEA, HPT::Type{<:AbstractFloat}; kwargs...) = solve(Solver(H, pa, ea, HPT; kwargs...), s)

function refineresult(r::EndgamerResult, solver)
    setup_workers(solver)
    if r.returncode == :success
        @unpack abstol, refinement_maxiters = solver.options
        @unpack H, cfg, cache = solver.pathtracker.low
        # TODO: we should switch to a higher precision if necessary
        # Since we have the winding number available
        # See the cauchy endgame test, the refinement is nearly useless...
        sol = copy(r.solution)
        try
            correct!(sol, 0.0, H, cfg, abstol, refinement_maxiters, cache)
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
