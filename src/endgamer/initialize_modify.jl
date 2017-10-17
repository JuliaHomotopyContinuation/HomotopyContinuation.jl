@inline function initialize(alg::AbstractEndgameAlgorithm, tracker::Pathtracker, x, R::Float64;
    geometric_series_factor=0.5,
    max_winding_number=8,
    abstol=tracker.options.abstol)
    predictions = Vector{typeof(tracker.low.x)}()
    xs::Vector{typeof(tracker.low.x)} = [x]

    iter = 0
    windingnumber = 1
    status = NotStarted
    failurecode = :default_failure_code

    cache = alg_cache(alg, tracker)
    options = EndgamerOptions(geometric_series_factor, abstol, max_winding_number)

    Endgamer(alg, cache, tracker, predictions, xs, iter, R, windingnumber,
        status, failurecode, options)
end

@inline function initialize(alg::AbstractEndgameAlgorithm, tracker::Pathtracker;
    kwargs...)
    initialize(alg, tracker, tracker.low.x, real(tracker.s))
end

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
