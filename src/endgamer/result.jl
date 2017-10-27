struct EndgamerResult{xType}
    retcode::Symbol
    solution::xType
    residual::Float64
    iterations::Int
    npredictions::Int
    homogenous_coordinate_magnitude::Float64
    windingnumber::Int
end


function EndgamerResult(endgamer::Endgamer)
    @unpack windingnumber, tracker = endgamer
    @unpack H, cfg = tracker.low


    if endgamer.status == Successfull
        retcode = :success
        solution = copy(endgamer.predictions[end])
    else
        solution = endgamer.xs[end]
        retcode = endgamer.failurecode
    end


    res = evaluate(H, solution, 0.0, cfg)
    residual = norm(res)

    if is_projective(tracker.alg)
        homogenous_coordinate_magnitude = convert(Float64, abs(first(solution)))
    else
        homogenous_coordinate_magnitude = 1.0
    end

    iterations = endgamer.iter
    npredictions = length(endgamer.predictions)

    EndgamerResult(
        retcode,
        solution,
        residual,
        iterations,
        npredictions,
        homogenous_coordinate_magnitude,
        windingnumber)
end

function EndgamerResult(endgamer::Endgamer, result::PathtrackerResult)
    windingnumber = 1
    residual = norm(evaluate(endgamer.tracker.low.H, result.solution, 0.0, endgamer.tracker.low.cfg))
    iterations = 0
    npredictions = 0
    EndgamerResult(
        result.retcode,
        result.solution,
        residual,
        iterations,
        npredictions,
        result.homogenous_coordinate_magnitude,
        windingnumber)
end
