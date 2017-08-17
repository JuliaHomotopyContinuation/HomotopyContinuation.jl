prepare_homotopy(H::AbstractHomotopy{T}, alg::AbstractPredictorCorrectorHomConAlgorithm) where {T<:Complex}= is_projective(alg) ? homogenize(H) : H
function prepare_start_value(H::AbstractHomotopy{T}, start_value::Vector{T}, alg::AbstractPredictorCorrectorHomConAlgorithm) where {T<:Complex}
    N = nvariables(H)
    if N == length(start_value)
        return start_value
    elseif N == length(start_value) + 1 && is_projective(alg)
        return [1.0; start_value]
    else
        if is_projective(alg)
            return error("A start_value has length $(length(start_value)). Excepted length $(N) or $(N-1).")
        else
            return error("A start_value has length $(length(start_value)). Excepted length $(N).")
        end
    end
end


struct Result{T<:Complex,Alg<:AbstractPredictorCorrectorHomConAlgorithm}
    solution::Vector{T}
    retcode::Symbol
    iterations::Int
    affine_result::Bool
    start_value::Vector{T}
    algorithm::Alg
    time_steps::Vector{Float64}
    trace::Vector{Vector{T}}
end


function solve(
    H::AbstractHomotopy{T},
    start_values,
    algorithm::AbstractPredictorCorrectorHomConAlgorithm;
    report_progress=true,
    kwargs...
) where {T<:Number}

    if eltype(start_values) != Vector{T}
        error("Expected as start values an iterable with elements of type Vector{$T} but got $(eltype(start_values))")
    end
    H = prepare_homotopy(H, algorithm)

    total_start_values = length(start_values)

    if report_progress println("Total number of paths to track: $total_start_values") end

    solutions = Vector{Result{T,typeof(algorithm)}}(total_start_values)
    for (index, start_value) in enumerate(start_values)
        start_value = prepare_start_value(H, start_value, algorithm)

        if report_progress
            println("Start to track path $index")
        end
        sol = track_path(H, start_value, algorithm; kwargs...)
        if report_progress println("Path $index returned: $(sol.retcode)") end

        solutions[index] = sol
    end

    solutions
end




"""
    track_path(homotopy, start_value, config)

Tracks the path from the given start_value. Reports additional infos in the result.
"""
function track_path(
    H::AbstractHomotopy{T},
    startvalue::Vector{T},
    algorithm::AbstractPredictorCorrectorHomConAlgorithm;
    tolerance_infinity=1-8,
    affine_result=true,
    trackpathkwargs...
    # #endgame
    # endgame_tolerance=1e-8,
    # endgame_start=0.1,
    # # precision_endgame::S=T
    # endgame_strategy=:none
) where {T<:Number}
    pathresult = trackpath(H, startvalue, algorithm, 1.0, 0.0; trackpathkwargs...)

    t = pathresult.laststep
    retcode = pathresult.returncode;
    x = pathresult.result
    k = pathresult.iterations
    steps = pathresult.steps
    trace = pathresult.trace


    if is_projective(algorithm)
        # check for solutions at infinity
        # assumes that the first variable is the artificial one
        at_infinity = norm(x[1]) < tolerance_infinity
        if affine_result
            sol = x[2:end] / x[1]
        else
            sol = x
        end
    else
        at_infinity = false
        sol = x
    end

    if at_infinity
        retcode = :AtInfinity
    end


    Result(sol, retcode, k, affine_result, startvalue, algorithm, steps, trace)
end
