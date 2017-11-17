export PathtrackerResult, solution

struct PathtrackerResult{T}
    retcode::Symbol
    solution::Vector{T}
    startvalue::Vector{T}
    residual::Float64
    iterations::Int
    homogenous_coordinate_magnitude::Float64

    # Extended analysis
    newton_residual::Float64
    condition_jacobian::Float64
end

"""
    solution(pathtracker)

Get `(retcode, solution)` from `pathtracker`. This is more lightwheight than a
`PathtrackerResult`.
"""
@inline function solution(tracker::Pathtracker)
    if tracker.iter ≥ tracker.options.maxiters
        retcode = :max_iterations
    elseif tracker.hit_singular_exception
        retcode = :singularity
    else
        retcode = :success
    end
    if tracker.usehigh
        solution = convert.(Complex{Low}, tracker.high.x)
    else
        solution = copy(tracker.low.x)
    end
    retcode, solution
end

function PathtrackerResult(tracker::Pathtracker{Low}, extended_analysis=true) where Low
    @unpack H, cfg = tracker.low
    retcode, sol = solution(tracker)

    res = evaluate(H, sol, tracker.s)
    residual = norm(res)

    if extended_analysis
        jacobian = Homotopy.jacobian(H, sol, tracker.s, cfg)

        newton_residual = norm(jacobian \ res)
        condition_jacobian = cond(jacobian)
    else
        newton_residual = NaN
        condition_jacobian = NaN
    end

    if is_projective(tracker.alg)
        homogenous_coordinate_magnitude = convert(Float64, abs(first(sol)))
    else
        homogenous_coordinate_magnitude = 1.0
    end


    PathtrackerResult(
        retcode,
        sol,
        copy(tracker.startvalue),
        residual,
        tracker.iter,
        homogenous_coordinate_magnitude,
        newton_residual,
        condition_jacobian)
end
