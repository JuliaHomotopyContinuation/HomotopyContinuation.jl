export Endgamer, EndgamerStatus

mutable struct EndgamerOptions
    geometric_series_factor::Float64
    abstol::Float64
    max_winding_number::Int
end

@enum EndgamerStatus NotStarted Started Successfull Failed

"""
    Endgamer(endgame_algorithm, pathtracker, [x, R]; kwargs...)

Construct an Endgamer object. The Endgamer 'plays' the endgame with the given `endgame_algorithm`
and uses the given `pathtracker` to move forward. The endgame is start at `x` to time `R`
(the endgame radius). In each iteration the endgame moves forward and then performs
one iteration of the endgame algorithm. In each iteration we could get another prediction
and an estimate of the winding number. Convergence is declared if two consecutive predictions
are smaller than a defined tolerance (`endgame_abstol`).

The following options are available:
* `geometric_series_factor=0.5`: The Endgamer moves forward using the geometric series
    ``λ^kR`` where ``λ`` is `geometric_series_factor`.
* `max_winding_number=16` the maximal winding number we allow. If we get a higher winding number
the path is declared failed.
* `endgame_abstol=pathtracker.options.abstol`: The tolerance necessary to declare convergence.

Endgamer supports similar to [`Pathtracker`](@ref) an iterator interface.
"""
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
    Endgamer(alg, tracker, tracker.low.x, convert(Float64, real(tracker.s)); kwargs...)
end

function Endgamer(alg::AbstractEndgameAlgorithm, tracker::Pathtracker{Low}, x, R::Float64;
    geometric_series_factor=0.5,
    max_winding_number=16,
    endgame_abstol=tracker.options.abstol) where Low
    predictions = Vector{Vector{Complex{Low}}}()
    xs::Vector{Vector{Complex{Low}}} = [x]
    startvalue = convert(Vector{Complex{Low}}, x)
    iter = 0
    windingnumber = 1
    status = NotStarted
    failurecode = :default_failure_code

    cache = alg_cache(alg, tracker)
    options = EndgamerOptions(geometric_series_factor, endgame_abstol, max_winding_number)

    Endgamer(alg, cache, tracker, predictions, xs, iter, R, windingnumber,
        status, failurecode, options, startvalue)
end

function is_endgamer_kwarg(kwarg)
    kwarg == :geometric_series_factor || 
    kwarg == :max_winding_number ||
    kwarg == :endgame_abstol
end
