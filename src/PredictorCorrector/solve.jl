function solve(
    H::AbstractHomotopy{T},
    startvalues,
    algorithm::APCA;
    verbose=false,
    kwargs...
) where {T<:Number}

    if eltype(startvalues) != Vector{T}
        error("Expected as start values an iterable with elements of type Vector{$T} " *
            "but got $(eltype(startvalues))")
    end

    H = prepare_homotopy(H, algorithm)
    total_startvalues = length(startvalues)

    if verbose
        println("Total number of paths to track: $total_startvalues")
    end

    solutions = Vector{Result{T}}(total_startvalues)
    for (index, startvalue) in enumerate(startvalues)
        if verbose
            println("Start to track path $index")
        end

        result = solve(H, startvalue, algorithm; homotopy_prepared=true, kwargs...)

        if verbose
            println("Path $index returned: $(result.returncode)")
        end

        solutions[index] = result
    end

    solutions
end

function solve(H::AbstractHomotopy{T}, startvalue::Vector{T}, algorithm::APCA;
    homotopy_prepared=false,
    tolerance_infinity=1e-8,
    trackpathkwargs...
) where {T<:Number}
    if !homotopy_prepared
        H = prepare_homotopy(H, algorithm)
    end
    x = prepare_startvalue(H, startvalue, algorithm)
    pathresult = trackpath(H, x, algorithm, 1.0, 0.0; trackpathkwargs...)
    postprocess(pathresult, startvalue, algorithm, tolerance_infinity)
end

function postprocess(pathresult::PathResult, startvalue, ::APCA{true}, tolerance_infinity)
    result = pathresult.result

    affine_solution = result[2:end] / result[1]
    projective_solution = result
    solution = length(startvalue) == length(affine_solution) ? affine_solution : projective_solution

    atinfinity = norm(result[1]) < tolerance_infinity
    returncode = atinfinity ? :AtInfinity : pathresult.returncode

    Result(solution,
           returncode,
           pathresult.iterations,
           affine_solution,
           projective_solution,
           startvalue,
           #pathresult.steps,
           pathresult.trace)
end
function postprocess(pathresult::PathResult, startvalue, ::APCA{false}, tolerance_infinity)
    solution = pathresult.result
    affine_solution = solution
    projective_solution = [1; solution]

    Result(solution,
           pathresult.returncode,
           pathresult.iterations,
           affine_solution,
           projective_solution,
           startvalue,
           #pathresult.steps,
           pathresult.trace)
end
