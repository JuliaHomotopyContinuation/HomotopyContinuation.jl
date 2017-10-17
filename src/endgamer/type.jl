mutable struct EndgamerOptions
    geometric_series_factor::Float64
    abstol::Float64
    max_winding_number::Int
end

@enum EndgamerStatus NotStarted Started Successfull Failed

mutable struct Endgamer{algType<:AbstractEndgameAlgorithm,Cache<:AbstractEndgameCache, xLow, Tracker<:Pathtracker}
    alg::algType
    cache::Cache
    tracker::Tracker

    predictions::Vector{xLow}
    xs::Vector{xLow}

    iter::Int
    R::Float64
    windingnumber::Int

    status::EndgamerStatus
    failurecode::Symbol

    options::EndgamerOptions
end
