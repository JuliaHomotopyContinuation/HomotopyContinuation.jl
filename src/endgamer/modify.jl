function reset!(endgamer::Endgamer, x, R)
    tracker = endgamer.tracker
    predictions = Vector{typeof(tracker.low.x)}()
    xs::Vector{typeof(tracker.low.x)} = [x]

    endgamer.R = R
    endgamer.iter = 0
    endgamer.windingnumber = 1
    endgamer.status = NotStarted
    endgamer.failurecode = :default_failure_code


    empty!(endgamer.xs)
    push!(endgamer.xs, x)

    empty!(endgamer.predictions)
end
