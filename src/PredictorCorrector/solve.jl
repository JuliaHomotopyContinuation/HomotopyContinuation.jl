prepare_homotopy(H::AbstractHomotopy{T}, alg::AbstractPredictorCorrectorHomConAlgorithm) where {T<:Complex}= is_projective(alg) ? homogenize(H) : H
function prepare_start_value(H::AbstractHomotopy{T}, start_value::Vector{T}, alg::AbstractPredictorCorrectorHomConAlgorithm) where {T<:Complex}
    N = nvars(H)
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

    if report_progress
        println("Total number of paths to track: $total_start_values")
    end

    solutions = Vector{Result{T,typeof(algorithm)}}(total_start_values)
    for (index, start_value) in enumerate(start_values)
        start_value = prepare_start_value(H, start_value, algorithm)

        if report_progress
            println("Start to track path $index")
        end
        sol = track_path(H, start_value, algorithm; kwargs...)
        println("Path $index returned: $(sol.retcode)")

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
    start_value::Vector{T},
    algorithm::AbstractPredictorCorrectorHomConAlgorithm;
    initial_step_length=0.01,
    tolerance=1e-8,
    successfull_steps_until_step_length_increase=3,
    step_length_increase_factor=2.0,
    step_length_decrease_factor=0.5,
    refine_solution_iterations=2,
    max_iterations=10000,
    tolerance_infinity=1e-4,
    affine_result=true,
    track_trace=false,
    iterations_correction_step=3,
    # endgame
    endgame_tolerance=1e-8,
    endgame_start=0.1,
    # precision_endgame::S=T
    endgame_strategy=:none
) where {T<:Complex}
    #some setup
    step_length = initial_step_length
    sucessive_successes = 0
    t = 1.0
    k = 0
    x = copy(start_value)
    u = similar(start_value)

    time_steps::Vector{Float64} = [];
    trace::Vector{Vector{T}} = [];

    # we only need to compute these once
    J_H = jacobian(H)
    ∂H∂t = ∂t(H)

    while t > 0 && k < max_iterations
        if track_trace
            push!(time_steps, t)
            push!(trace, x)
        end

        Δt = min(step_length, t)

        # if t <= endgame_start
        #   tol = tolerance
        # else
        tol = tolerance
        u .= predict(algorithm, H, J_H, ∂H∂t, x, t, Δt)
        converged = correct!(u, algorithm, H, J_H, x, t - Δt, tol, iterations_correction_step)
        if converged
                x .= u
                t -= Δt
                sucessive_successes += 1
                if sucessive_successes == successfull_steps_until_step_length_increase
                    step_length *= step_length_increase_factor
                    sucessive_successes = 0
                end
            else
                sucessive_successes = 0
                step_length *= step_length_decrease_factor
            end
        # end
        k += 1
    end

    # now we refine our solution
    if t ≈ 0
        refinement_converged = correct!(u, algorithm, H, J_H, x, 0.0, endgame_tolerance, iterations_correction_step)
        if refinement_converged
            x .= u
        end
    end

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

    converged = t ≈ 0 && norm(evaluate(H, x, 0.0)) < endgame_tolerance

    if k >= max_iterations
        retcode = :MaxIterations
    elseif !converged
        retcode = :NotConverged
    elseif at_infinity
        retcode = :AtInfinity
    else
        retcode = :Success
    end


    Result(sol, retcode, k, affine_result, start_value, algorithm, time_steps, trace)
end