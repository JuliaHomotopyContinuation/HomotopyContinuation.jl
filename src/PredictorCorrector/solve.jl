"""
    solve(H::AbstractHomotopy{T}, startvalue::Vector{T}, [algorithm,] endgame_start=0.1, kwargs...)

Track the path ``x(t)`` implicitly defined by ``H(x(t),t)`` from `t=1` to `t=0` where
`x(1)=startvalue`. Returns a [`Result`](@ref) instance.

It uses the defined `algorihm` for the prediction and correction step.
If no `algorithm` is defined the in general best performing algorithm is choosen. Currently
this is `Spherical`.

## Optional arguments
* `endgame_start=0.1`: When to switch to the Cauchy endgame, for a value of t=0.0 no endgame happens.
*  optional arguments to [`pathtracking`](@ref)
*  optional arguments to [`cauchyendgame`](@ref)

    solve(H::AbstractHomotopy, startvalues, algorithm, report_progess=false, kwargs...)

Track multiple paths with the given `startvalues`. Returns a vector of [`Result`](@ref) instances,
one per start value.
If `report_progess=true` the current progess will be logged.
"""
solve(H::AbstractHomotopy, startvalue_s; kwargs...) = solve(H, startvalue_s, Spherical(); kwargs...)
function solve(
    H::AbstractHomotopy{T},
    startvalues,
    algorithm::APCA;
    report_progess=false,
    kwargs...
) where {T<:Number}

    if eltype(startvalues) != Vector{T}
        error("Expected as start values an iterable with elements of type Vector{$T} " *
            "but got $(eltype(startvalues))")
    end

    H = prepare_homotopy(H, algorithm)
    total_startvalues = length(startvalues)

    if report_progess
        println("Total number of paths to track: $total_startvalues")
    end

    solutions = Vector{Result{T}}(total_startvalues)
    for (index, startvalue) in enumerate(startvalues)
        if report_progess
            println("Start to track path $index")
        end

        result = solve(H, startvalue, algorithm; homotopy_prepared=true, kwargs...)

        if report_progess
            println("Path $index returned: $(result.returncode)")
        end

        solutions[index] = result
    end

    solutions
end

function solve(H::AbstractHomotopy{T}, startvalue::Vector{T}, algorithm::APCA{false};
    homotopy_prepared=false,
    trackpathkwargs...
) where {T<:Number}
    if !homotopy_prepared
        H = prepare_homotopy(H, algorithm)
    end
    x = prepare_startvalue(H, startvalue, algorithm)
    pathresult = trackpath(H, x, algorithm, 1.0, 0.0; trackpathkwargs...)
    postprocess(pathresult, startvalue, algorithm)
end

function solve(H::AbstractHomotopy{T}, startvalue::Vector{T}, algorithm::APCA{true};
    homotopy_prepared=false,
    endgame_start=0.1,
    # endgame kwargs copied from cauchyendgame.jl
    prediction_tolerance=1e-4,
    geometric_series_factor=0.5,
    endgame_tolerance=1e-8,
    samples_per_loop::Int=8,
    max_winding_number::Int=8,
    loopclosed_tolerance=1e-3,
    tolerance_infinity=1e-8,
    trackpathkwargs...
) where {T<:Number}
    if !homotopy_prepared
        H = prepare_homotopy(H, algorithm)
    end
    x = prepare_startvalue(H, startvalue, algorithm)
    pathresult = trackpath(H, x, algorithm, 1.0, endgame_start; trackpathkwargs...)
    if nosuccess(pathresult) || (endgame_start <= 0.0)
        return postprocess(pathresult, startvalue, algorithm, tolerance_infinity)
    end
    endgameresult = cauchyendgame(H, pathresult.result, endgame_start, algorithm,
        # this is really not optimal but don't know a better way for now...
        prediction_tolerance=prediction_tolerance,
        geometric_series_factor=geometric_series_factor,
        endgame_tolerance=endgame_tolerance,
        samples_per_loop=samples_per_loop,
        max_winding_number=max_winding_number,
        loopclosed_tolerance=loopclosed_tolerance,
        tolerance_infinity=tolerance_infinity;
        trackpathkwargs...
        )
    postprocess(pathresult, endgameresult, startvalue, algorithm, tolerance_infinity)
end

function postprocess(pathresult::PathResult{T}, startvalue, ::APCA{true}, tolerance_infinity) where T
    result = pathresult.result

    affine_solution = result[2:end] / result[1]
    projective_solution = result
    solution = length(startvalue) == length(affine_solution) ? affine_solution : projective_solution

    returncode = atinfinity(result, tolerance_infinity) ? :AtInfinity : pathresult.returncode

    Result(solution,
           returncode,
           pathresult.iterations,
           0, #endgame_iterations
           affine_solution,
           projective_solution,
           startvalue,
           #pathresult.steps,
           pathresult.trace,
           Nullable{ConvergentCluster{T}}())
end

function postprocess(pathresult::PathResult, endgameresult::CauchyEndgameResult, startvalue, ::APCA{true}, tolerance_infinity)
    result = endgameresult.result

    affine_solution = result[2:end] / result[1]
    projective_solution = result
    solution = length(startvalue) == length(affine_solution) ? affine_solution : projective_solution

    returncode = atinfinity(result, tolerance_infinity) ? :AtInfinity : endgameresult.returncode

    Result(solution,
           returncode,
           pathresult.iterations,
           endgameresult.iterations,
           affine_solution,
           projective_solution,
           startvalue,
           #pathresult.steps,
           [pathresult.trace; endgameresult.trace],
           endgameresult.convergent_cluster)
end

function postprocess(pathresult::PathResult{T}, startvalue, ::APCA{false}) where T
    solution = pathresult.result
    affine_solution = solution
    projective_solution = [1; solution]

    Result(solution,
           pathresult.returncode,
           pathresult.iterations,
           0, #endgame_iterations
           affine_solution,
           projective_solution,
           startvalue,
           #pathresult.steps,
           pathresult.trace,
           Nullable{ConvergentCluster{T}}())
end
