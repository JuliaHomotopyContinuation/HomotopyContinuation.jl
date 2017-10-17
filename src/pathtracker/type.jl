# This are the options set by the user
mutable struct PathtrackerOptions
    maxiters::Int
    verbose::Bool
    initial_steplength::Float64
    consecutive_successfull_steps_until_steplength_increase::Int
    steplength_increase_factor::Float64
    steplength_decrease_factor::Float64
    abstol::Float64
    corrector_maxiters::Int
end

# This holds all structures which depend on the given precision
mutable struct PathtrackerPrecisionValues{
    T<:Real,
    HT<:AbstractHomotopy{Complex{T}},
    F1,
    F2,
    xType<:AbstractVector{Complex{T}},
    Cache<:AbstractPathtrackerCache{Complex{T}}}
    H::HT
    J_H!::F1
    Hdt!::F2
    x::xType
    xnext::xType
    cache::Cache
end

# ATTENTION:
# If you change something here you probably also want to change `initialize` and `reset!` !
mutable struct Pathtracker{
    LowPrecision<:Real,
    HighPrecision<:Real,
    xLowType<:AbstractVector{Complex{LowPrecision}},
    xHighType<:AbstractVector{Complex{HighPrecision}},
    algType<:AbstractPathtrackingAlgorithm,
    LowHT<:AbstractHomotopy{Complex{LowPrecision}},
    HighHT<:AbstractHomotopy{Complex{HighPrecision}},
    F1, F2, F3, F4,
    LowCache<:AbstractPathtrackerCache{Complex{LowPrecision}},
    HighCache<:AbstractPathtrackerCache{Complex{HighPrecision}},
    O<:PathtrackerOptions}

    alg::algType # the used pathtracking algorithm

    #startvalue::xLowType
    # low precision values
    low::PathtrackerPrecisionValues{LowPrecision, LowHT, F1, F2, xLowType, LowCache}
    # higher precision values
    high::PathtrackerPrecisionValues{HighPrecision, HighHT, F3, F4, xHighType, HighCache}

    usehigh::Bool

    iter::Int
    t::Float64 # from 1 to 0
    steplength::Float64 # steplength for t
    s::Complex{LowPrecision} # value on the path from start to target
    sdiff::Complex{LowPrecision}
    ds::Complex{LowPrecision} # steplength for t
    snext::Complex{LowPrecision}
    step_sucessfull::Bool
    consecutive_successfull_steps::Int
    options::O
end
