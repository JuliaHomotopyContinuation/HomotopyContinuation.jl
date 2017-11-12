export Endgamer, EndgamerStatus

mutable struct EndgamerOptions
    geometric_series_factor::Float64
    abstol::Float64
    max_winding_number::Int
end

@enum EndgamerStatus NotStarted Started Successfull Failed

mutable struct Endgamer{Low, algType<:AbstractEndgameAlgorithm,Cache<:AbstractEndgameCache, Tracker<:Pathtracker{Low}}
    alg::algType
    cache::Cache
    tracker::Tracker

    predictions::Vector{Vector{Complex{Low}}}
    xs::Vector{Vector{Complex{Low}}}

    iter::Int
    R::Float64
    windingnumber::Int

    status::EndgamerStatus
    failurecode::Symbol

    options::EndgamerOptions

    startvalue::Vector{Complex{Low}}
end

function Endgamer(alg::AbstractEndgameAlgorithm, tracker::Pathtracker;
    kwargs...)
    Endgamer(alg, tracker, tracker.low.x, convert(Float64, real(tracker.s)))
end

function Endgamer(alg::AbstractEndgameAlgorithm, tracker::Pathtracker{Low}, x, R::Float64;
    geometric_series_factor=0.5,
    max_winding_number=8,
    abstol=tracker.options.abstol) where Low
    predictions = Vector{Vector{Complex{Low}}}()
    xs::Vector{Vector{Complex{Low}}} = [x]
    startvalue = convert(Vector{Complex{Low}}, x)
    iter = 0
    windingnumber = 1
    status = NotStarted
    failurecode = :default_failure_code

    cache = alg_cache(alg, tracker)
    options = EndgamerOptions(geometric_series_factor, abstol, max_winding_number)

    Endgamer(alg, cache, tracker, predictions, xs, iter, R, windingnumber,
        status, failurecode, options, startvalue)
end

function is_endgamer_kwarg(kwarg)
    kwarg == :geometric_series_factor || 
    kwarg == :max_winding_number ||
    kwarg == :abstol
end
