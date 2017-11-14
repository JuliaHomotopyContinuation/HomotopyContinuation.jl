export PathtrackerResult, solution

struct PathtrackerResult{T}
    returncode::Symbol
    solution::Vector{T}
    startvalue::Vector{T}
    residual::Float64
    iterations::Int
    angle_to_infinity::Float64

    # Extended analysis
    newton_residual::Float64
    condition_jacobian::Float64
end

"""
    solution(pathtracker)

Get `(returncode, solution)` from `pathtracker`. This is more lightwheight than a
`PathtrackerResult`.
"""
@inline function solution(tracker::Pathtracker)
    if tracker.iter ≥ tracker.options.maxiters
        returncode = :max_iterations
    else
        returncode = :success
    end
    if tracker.usehigh
        solution = convert.(Complex{Low}, tracker.high.x)
    else
        solution = copy(tracker.low.x)
    end
    returncode, solution
end

function PathtrackerResult(tracker::Pathtracker{Low}, extended_analysis=true) where Low
    @unpack H, cfg = tracker.low
    if tracker.usehigh
        solution = convert.(Complex{Low}, tracker.high.x)
    else
        solution = copy(tracker.low.x)
    end

    res = evaluate(H, solution, tracker.s)
    residual = norm(res)

    if tracker.iter ≥ tracker.options.maxiters
        returncode = :max_iterations
    else
        returncode = :success
    end

    if extended_analysis
        jacobian = Homotopy.jacobian(H, solution, tracker.s, cfg)

        newton_residual = norm(jacobian \ res)
        condition_jacobian = cond(jacobian)
    else
        newton_residual = NaN
        condition_jacobian = NaN
    end

    if is_projective(tracker.alg)
        angle_to_infinity = convert(Float64, abs(first(solution)))
    else
        angle_to_infinity = 1.0
    end


    PathtrackerResult(
        returncode,
        solution,
        copy(tracker.startvalue),
        residual,
        tracker.iter,
        angle_to_infinity,
        newton_residual,
        condition_jacobian)
end
