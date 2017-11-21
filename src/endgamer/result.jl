export EndgamerResult

"""
    EndgamerResult(endgamer, extended_analysis=false)

Reads the result from the current pathtracker state.
A `EndgamerResult` contains:
* `returncode`: One of `:ill_conditioned_zone`, `:success`, `:windingnumber_too_high`
* `solution::Vector{T}`: The solution.
* `startvalue::Vector{T}`: The solution.
* `residual::Float64`: The value of the infinity norm of `H(solution, 0)`.
* `iterations`: The number of iterations the pathtracker needed.
* `npredictions`: The number of predictions
* `predictions`: All predictions for further introspection (e.g. the path failed)
* `angle_to_infinity`: The angle to infinity is the angle of the solution to the hyperplane where the homogenizing coordinate is ``0``.
* `windingnumber`: The estimated winding number
"""
struct EndgamerResult{T}
    returncode::Symbol
    solution::Vector{T}
    startvalue::Vector{T}
    residual::Float64
    iterations::Int
    npredictions::Int
    predictions::Vector{Vector{T}}
    angle_to_infinity::Float64
    windingnumber::Int
end


function EndgamerResult(endgamer::Endgamer)
    @unpack windingnumber, tracker = endgamer
    @unpack H, cfg = tracker.low


    if endgamer.status == Successfull
        returncode = :success
        solution = copy(endgamer.predictions[end])
    else
        solution = endgamer.xs[end]
        returncode = endgamer.failurecode
    end


    res = evaluate(H, solution, 0.0, cfg)
    residual = convert(Float64, norm(res))

    if is_projective(tracker.alg)
        angle_to_infinity = convert(Float64, abs(first(solution)))
    else
        angle_to_infinity = 1.0
    end

    iterations = endgamer.iter
    npredictions = length(endgamer.predictions)

    EndgamerResult(
        returncode,
        solution,
        copy(endgamer.startvalue),
        residual,
        iterations,
        npredictions,
        deepcopy(endgamer.predictions),
        angle_to_infinity,
        windingnumber)
end

function EndgamerResult(endgamer::Endgamer, result::PathtrackerResult)
    windingnumber = 1
    residual = norm(evaluate(endgamer.tracker.low.H, result.solution, 0.0, endgamer.tracker.low.cfg))
    iterations = 0
    npredictions = 0
    EndgamerResult(
        result.returncode,
        result.solution,
        copy(endgamer.startvalue),
        residual,
        iterations,
        npredictions,
        Vector{typeof(result.solution)}(),
        result.angle_to_infinity,
        windingnumber)
end

function refined_solution(r::EndgamerResult{T}, refined_sol::Vector{T}) where T
    EndgamerResult(
        r.returncode,
        refined_sol,
        r.startvalue,
        r.residual,
        r.iterations,
        r.npredictions,
        deepcopy(r.predictions),
        r.angle_to_infinity,
        r.windingnumber)
end
