export Result, solution, returncode, iterations, endgame_iterations, algorithm,
    affine_solution, projective_solution, startvalue, pathtrace, pathsteps,
    convergent_cluster, issuccessfull


"""
    Result(solution, returncode, iterations, affine_solution, projective_solution, startvalue, steps, trace, nullable_cluster)

Construct a result of `solve` with an `AbstractPredictorCorrectorAlgorithm`.

## Fields
* `solution::Vector{T}`: This has the same length as startvalue
* `returncode::Symbol`: :Success, :AtInfinity or any return code from [CauchyEndgameResult](@ref) and [PathResult](@ref).
* `iterations::Int`
* `endgame_iterations::Int`
* `affine_solution::Vector{T}`
* `project_solution::Vector{T}`
* `startvalue::Vector{T}`
* `trace::Vector{Vector{T}}`: Empty if `record_trace=false` (default)
* `convergent_cluster::Nullable{ConvergentCluster{T}}`: The convergent cluster of the endgame.
"""
struct Result{T<:Number, Alg<:AbstractPredictorCorrectorAlgorithm}
    solution::Vector{T}
    returncode::Symbol
    iterations::Int
    endgame_iterations::Int

    algorithm::Alg

    affine_solution::Vector{T}
    projective_solution::Vector{T}

    startvalue::Vector{T}
    #* `steps::Vector{S}`: Empty if `record_steps=false` (default)
    trace::Vector{Vector{T}}
    steps::Vector{Float64}

    convergent_cluster::Nullable{ConvergentCluster{T}}
end


function Result(pathresult::PathResult{T}, startvalue::Vector{T}, alg::APCA{Val{true}}, tolerance_infinity) where T
    result = pathresult.result

    affine_solution = result[2:end] / result[1]
    projective_solution = result
    solution = length(startvalue) == length(affine_solution) ? affine_solution : projective_solution
    returncode = atinfinity(result, tolerance_infinity) ? :AtInfinity : pathresult.returncode

    Result(solution,
           returncode,
           pathresult.iterations,
           0, #endgame_iterations,
           alg,
           affine_solution,
           projective_solution,
           startvalue,
           pathresult.trace,
           pathresult.steps,
           Nullable{ConvergentCluster{T}}())
end

function Result(pathresult::PathResult, endgameresult::CauchyEndgameResult, startvalue::Vector, alg::APCA{Val{true}}, tolerance_infinity::Float64)
    result = endgameresult.result

    affine_solution = result[2:end] / result[1]
    projective_solution = result
    solution = length(startvalue) == length(affine_solution) ? affine_solution : projective_solution

    returncode = atinfinity(result, tolerance_infinity) ? :AtInfinity : endgameresult.returncode

    Result(solution,
           returncode,
           pathresult.iterations,
           endgameresult.iterations,
           alg,
           affine_solution,
           projective_solution,
           startvalue,
           [pathresult.trace; endgameresult.trace],
           [pathresult.steps; endgameresult.steps],
           endgameresult.convergent_cluster)
end

atinfinity(x, tolerance_infinity) = norm(normalize(x)[1]) < tolerance_infinity

function Result(pathresult::PathResult{T}, startvalue, alg::APCA{Val{false}}) where T
    solution = pathresult.result
    affine_solution = solution
    projective_solution = [1; solution]

    Result(solution,
           pathresult.returncode,
           pathresult.iterations,
           0, #endgame_iterations
           alg,
           affine_solution,
           projective_solution,
           startvalue,
           pathresult.trace,
           pathresult.steps,
           Nullable{ConvergentCluster{T}}())
end


solution(r::Result) = r.solution
returncode(r::Result) = r.returncode
iterations(r::Result) = r.iterations
endgame_iterations(r::Result) = r.endgame_iterations
algorithm(r::Result) = r.algorithm
affine_solution(r::Result) = r.affine_solution
projective_solution(r::Result) = r.projective_solution
startvalue(r::Result) = r.startvalue
pathtrace(r::Result) = r.trace
pathsteps(r::Result) = r.steps
convergent_cluster(r::Result) = r.convergent_cluster


issuccessfull(r::Result) = returncode(r) == :Success


#Base.show(io::IO, ::MIME"text/plain", res::Result) = printresult(io, res)
Base.show(io::IO, res::Result) = printresult(io, res)
function printresult(io::IO, res::Result)
    println(io, typeof(res),":")
    println(io, "------------------------------")
    println(io, "* solution: ", solution(res))
    println(io, "* returncode: ", returncode(res))
    println(io, "------------------------------")
    println(io, "* iterations: ", iterations(res))
    println(io, "* endgame_iterations: ", endgame_iterations(res))
    println(io, "* affine_solution: ", affine_solution(res))
    println(io, "* projective_solution: ", projective_solution(res))
    println(io, "* startvalue: ", startvalue(res))
    println(io, "* pathsteps: ", length(pathsteps(res)), " entries")
    println(io, "* pathtrace: ", length(pathtrace(res)), " entries")
    if !isnull(convergent_cluster(res))
        println(io, "* convergent cluster: ", get(map(string, convergent_cluster(res))))
    end
end
