export EndgamerResult

struct EndgamerResult{T}
    returncode::Symbol
    solution::Vector{T}
    startvalue::Vector{T}
    residual::Float64
    iterations::Int
    npredictions::Int
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
        r.angle_to_infinity,
        r.windingnumber)
end
