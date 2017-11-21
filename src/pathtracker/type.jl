export Pathtracker

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

function PathtrackerOptions(;maxiters::Int=10_000,
    verbose=false,
    initial_steplength::Float64=0.1,
    consecutive_successfull_steps_until_steplength_increase::Int=3,
    steplength_increase_factor::Float64=2.0,
    steplength_decrease_factor::Float64=inv(steplength_increase_factor),
    path_precision::Float64=1e-6,
    corrector_maxiters::Int=3)

    PathtrackerOptions(
        maxiters, verbose, initial_steplength,
        consecutive_successfull_steps_until_steplength_increase,
        steplength_increase_factor, steplength_decrease_factor,
        path_precision, corrector_maxiters)
end

function is_pathtracker_options_kwarg(kwarg)
    kwarg == :maxiters || 
    kwarg == :verbose ||
    kwarg == :initial_steplength ||
    kwarg == :consecutive_successfull_steps_until_steplength_increase ||
    kwarg == :steplength_increase_factor ||
    kwarg == :steplength_decrease_factor ||
    kwarg == :path_precision ||
    kwarg == :corrector_maxiters
end

# This holds all structures which depend on the given precision
mutable struct PathtrackerPrecisionValues{
    T<:Real,
    HT<:AbstractHomotopy{Complex{T}},
    Config<:AbstractHomotopyConfig{Complex{T}},
    #xType<:AbstractVector{Complex{T}},
    Cache<:AbstractPathtrackerCache{Complex{T}}}
    H::HT
    cfg::Config
    x::Vector{Complex{T}}
    xnext::Vector{Complex{T}}
    cache::Cache
end


"""
    Pathtracker(H::AbstractHomotopy{T}, alg, [HT::Type=widen(T)]; kwargs...)

Construct a Pathtracker object. This contains all informations to track a single path
for `H` with the given pathtracking algorithm `alg`. The optional type `HT`
is used if the pathracker decides to switch to a high precision mode.

The following keyword arguments are supported:

* `path_precision=1e-6`: The precision for which a correction step is decleared successfull.
* `corrector_maxiters=3`: The maximal number of correction iterations. A higher value as 3 is not recommended.
* `initial_steplength=0.1`: The initial steplength a preditor-corrector algorithm uses.
* `consecutive_successfull_steps_until_steplength_increase=3`:
    The number of consecutive successfull steps until the step length is increased multiplied
    with the factor `steplength_increase_factor`.
* `steplength_increase_factor=2.0`
* `steplength_decrease_factor=inv(steplength_increase_factor)`: If a correction step fails the step length is multiplied
    with this factor.
* `maxiters=10_000`: The maximum number of iterations.
* `vebose=false`: Print additional informations / warnings during the computation.
"""
# ATTENTION:
# If you change something here you probably also want to change `initialize` and `reset!` !
mutable struct Pathtracker{
    LowPrecision<:Real,
    HighPrecision<:Real,
    algType<:AbstractPathtrackingAlgorithm,
    LowHT<:AbstractHomotopy{Complex{LowPrecision}},
    HighHT<:AbstractHomotopy{Complex{HighPrecision}},
    LowConfig<:AbstractHomotopyConfig{Complex{LowPrecision}},
    HighConfig<:AbstractHomotopyConfig{Complex{HighPrecision}},
    LowCache<:AbstractPathtrackerCache{Complex{LowPrecision}},
    HighCache<:AbstractPathtrackerCache{Complex{HighPrecision}},
    O<:PathtrackerOptions}

    alg::algType # the used pathtracking algorithm

    #startvalue::xLowType
    # low precision values
    low::PathtrackerPrecisionValues{LowPrecision, LowHT, LowConfig, LowCache}
    # higher precision values
    high::PathtrackerPrecisionValues{HighPrecision, HighHT, HighConfig, HighCache}
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
    status::Symbol

    startvalue::Vector{Complex{LowPrecision}}
end

is_pathtracker_kwarg(kwarg) = is_pathtracker_options_kwarg(kwarg)


function Pathtracker(
    H::AbstractHomotopy{Complex{T}},
    alg::algType,
    x0=zeros(Complex{T}, nvariables(H)),
    s_start=1.0,
    s_end=0.0, HT::Type{S}=widen(T); kwargs...) where {algType<:AbstractPathtrackingAlgorithm, T<:AbstractFloat, S<:AbstractFloat}

    options = PathtrackerOptions(;kwargs...)

    x = Vector{Complex{T}}(length(x0))
    x .= x0
    if is_projective(alg)
        H = homogenize(H)
        embed_projective_if_necessary!(x, H)
    end

    usehigh = false
    iter = 0
    t = 1.0
    startvalue = copy(x)
    steplength = options.initial_steplength
    s = convert(Complex{T}, s_start)
    sdiff = convert(Complex{T}, s_end) - s
    ds = steplength * sdiff
    snext = s
    step_sucessfull = false
    consecutive_successfull_steps = 0

    # Assemble the precision dependent types
    cfg = Homotopies.config(H)
    xnext = similar(x)
    cache = alg_cache(alg, H, x)
    low = PathtrackerPrecisionValues{T, typeof(H), typeof(cfg), typeof(cache)}(H, cfg, x, xnext, cache)

    highH = convert(promote_type(typeof(H), Complex{HT}), H)
    highcfg = Homotopies.config(highH)
    highx0 = convert(Vector{Complex{HT}}, x)
    highxnext = similar(highx0)
    highcache = alg_cache(alg, highH, highx0)
    high = PathtrackerPrecisionValues{
        HT, typeof(highH), typeof(highcfg), typeof(highcache)
        }(highH, highcfg, highx0, highxnext, highcache)

    status = :default

    Pathtracker{T, HT, algType,
        typeof(low.H), typeof(high.H),
        typeof(low.cfg), typeof(high.cfg),
        typeof(low.cache), typeof(high.cache),
        typeof(options)
        }(alg, low, high, usehigh,
        iter, t, steplength, s, sdiff, ds, snext,
        step_sucessfull, consecutive_successfull_steps, options,
        status, copy(x))
end

function Pathtracker(
    H::AbstractHomotopy{Complex{T}},
    alg::algType,
    HT::Type{S}; kwargs...) where {algType<:AbstractPathtrackingAlgorithm, T<:AbstractFloat, S<:AbstractFloat}
    Pathtracker(H, alg, zeros(Complex{T}, nvariables(H)), 1.0, 0.0, HT; kwargs...)
end
