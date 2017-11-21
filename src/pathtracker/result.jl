export PathtrackerResult, solution

"""
    PathtrackerResult(pathtracker, extended_analysis=false)

Reads the result from the current pathtracker state.
A `PathtrackerResult` contains:
* `returncode`: One of `:max_iterations`, `:singularity`, `:invalid_startvalue`, `:success`.
* `solution::Vector{T}`: The solution.
* `residual::Float64`: The value of the infinity norm of `H(solution, 0)`.
* `iterations`: The number of iterations the pathtracker needed.
* `angle_to_infinity`: The angle to infinity is the angle of the solution to the hyperplane where the homogenizing coordinate is ``0``.
If `extended_analysis=true` there is also:
* `newton_residual`: The value of the 2-norm of ``J_H(\\text{solution})^{-1}H(\\text{solution}, 0)``
* `condition_number`: A high condition number indicates singularty. See [`Homotopy.κ`](@ref) for details.
"""
struct PathtrackerResult{T}
    returncode::Symbol
    solution::Vector{T}
    startvalue::Vector{T}
    residual::Float64
    iterations::Int
    angle_to_infinity::Float64

    # Extended analysis
    newton_residual::Float64
    condition_number::Float64
end

"""
    solution(pathtracker)

Get `(returncode, solution)` from `pathtracker`. This is more lightwheight than a
`PathtrackerResult`.
"""
@inline function solution(tracker::Pathtracker{Low}) where {Low}
    if tracker.iter ≥ tracker.options.maxiters
        returncode = :max_iterations
    elseif tracker.status == :singular_exception
        returncode = :singularity
    elseif tracker.status == :invalid_startvalue
        returncode = :invalid_startvalue
    else
        returncode = :success
    end
    if tracker.usehigh
        sol = convert.(Complex{Low}, tracker.high.x)
    else
        sol = copy(tracker.low.x)
    end
    returncode, sol
end

function PathtrackerResult(tracker::Pathtracker{Low}, extended_analysis=true) where Low
    @unpack H, cfg = tracker.low
    returncode, sol = solution(tracker)

    res = evaluate(H, sol, tracker.s)
    residual = convert(Float64, norm(res))

    if extended_analysis
        jacobian = Homotopy.jacobian(H, sol, tracker.s, cfg)

        newton_residual = convert(Float64, norm(jacobian \ res))
        condition_number =  convert(Float64, Homotopy.κ(H, sol, 0.0, cfg))
    else
        newton_residual = NaN
        condition_number = NaN
    end

    if is_projective(tracker.alg)
        angle_to_infinity = convert(Float64, abs(first(sol)))
    else
        angle_to_infinity = 1.0
    end


    PathtrackerResult(
        returncode,
        sol,
        copy(tracker.startvalue),
        residual,
        tracker.iter,
        angle_to_infinity,
        newton_residual,
        condition_number)
end
