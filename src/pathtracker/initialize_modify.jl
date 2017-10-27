
@inline function initialize(
    alg::algType,
    H::AbstractHomotopy{Complex{T}},
    x0::Vector{Complex{T}} = zeros(Complex{T}, nvariables(H)),
    s_start=1.0,
    s_end=0.0;
    highprecisiontype = _widen(T),
    maxiters::Int=10_000,
    verbose::Bool=false,
    consecutive_successfull_steps_until_steplength_increase::Int=3,
    steplength_increase_factor::Float64=2.0,
    steplength_decrease_factor::Float64=inv(steplength_increase_factor),
    abstol::Float64=1e-6,
    corrector_maxiters::Int=3,
    initial_steplength::Float64=0.1) where {algType<:AbstractPathtrackingAlgorithm, T<:Real}

    x = copy(x0)
    if is_projective(alg)
        H = homogenize(H)
        embed_projective_if_necessary!(x, H)
    end

    options = PathtrackerOptions(
        maxiters,
        verbose,
        initial_steplength,
        consecutive_successfull_steps_until_steplength_increase,
        steplength_increase_factor,
        steplength_decrease_factor,
        abstol,
        corrector_maxiters
        )

    usehigh = false
    iter = 0
    t = 1.0
    startvalue = copy(x)
    steplength = initial_steplength
    s = convert(Complex{T}, s_start)
    sdiff = convert(Complex{T}, s_end) - s
    ds = steplength * sdiff
    snext = s
    step_sucessfull = false
    consecutive_successfull_steps = 0

    # Assemble the precision dependent types
    cfg = Homotopy.config(H)
    xnext = similar(x)
    cache = alg_cache(alg, H, x)
    low = PathtrackerPrecisionValues{
        T, typeof(H), typeof(cfg), Vector{Complex{T}}, typeof(cache)
        }(H, cfg, x, xnext, cache)

    highH = convert(promote_type(typeof(H), Complex{highprecisiontype}), H)
    highcfg = Homotopy.config(highH)
    highx0 = convert(Vector{Complex{highprecisiontype}}, x)
    highxnext = similar(highx0)
    highcache = alg_cache(alg, highH, highx0)
    high = PathtrackerPrecisionValues{
        highprecisiontype, typeof(highH), typeof(highcfg),
        Vector{Complex{highprecisiontype}}, typeof(highcache)
        }(highH, highcfg, highx0, highxnext, highcache)

    Pathtracker{T, highprecisiontype, Vector{Complex{T}}, Vector{Complex{highprecisiontype}}, algType,
        typeof(low.H), typeof(high.H),
        typeof(low.cfg), typeof(high.cfg),
        typeof(low.cache), typeof(high.cache),
        typeof(options)
        }(alg, low, high, usehigh,
        iter, t, steplength, s, sdiff, ds, snext,
        step_sucessfull, consecutive_successfull_steps, options)
end

"""
    embed_projective_if_necessary(x, H)

Embeds a vector into the projective space if necessary, i.e. if it's length is one less
than the number of variables of `H`. `H` is assumed to be homogenized. After the (eventual)
embedding the value is normalized.
"""
function embed_projective_if_necessary!(x::AbstractVector{T}, H::AbstractHomotopy{T}) where T
    N = Homotopy.nvariables(H)
    n = length(x)
    if N - 1 == n
        unshift!(x, one(T))
    elseif N != n
        error("A start value has length $n. Excepted length $N or $(N-1).")
    end
    x
end

function reset!(tracker::Pathtracker{Low}, x0, s_start, s_end) where Low
    # TODO: we should carry over the original options
    setprecisionvalues!(tracker, x0)

    #tracker.startvalue = copy(tracker.low.x)

    tracker.iter = 0
    tracker.t = 1.0
    tracker.steplength = tracker.options.initial_steplength
    tracker.s = convert(Complex{Low}, s_start)
    tracker.sdiff = convert(Complex{Low}, s_end) - tracker.s
    tracker.ds = tracker.steplength * tracker.sdiff
    tracker.snext = tracker.s
    tracker.step_sucessfull = false
    tracker.consecutive_successfull_steps = 0
    tracker
end

function setprecisionvalues!(tracker::Pathtracker{Low, High}, x0::AbstractVector{Complex{Low}}) where {Low, High}
    N = length(tracker.low.x)
    if length(x0) == N - 1
        tracker.low.x[2:end] .= x0
        tracker.low.x[1] = one(Low)
    elseif length(x0) == N
        tracker.low.x .= x0
    else
        error("A start value has length $n. Excepted length $N or $(N-1).")
    end
    tracker.low.xnext .= tracker.low.x
    tracker.high.x .= tracker.low.x #this automatically converts to High
    tracker.high.xnext .= tracker.high.x
    tracker.usehigh = false
end

function setprecisionvalues!(tracker::Pathtracker{Low, High}, x0::AbstractVector{Complex{High}}) where {Low, High}
    N = length(tracker.high.x)
    if length(x0) == N - 1
        tracker.high.x[2:end] .= x0
        tracker.high.x[1] = one(High)
    elseif length(x0) == N
        tracker.high.x .= x0
    else
        error("A start value has length $n. Excepted length $N or $(N-1).")
    end
    tracker.high.xnext .= tracker.high.x
    tracker.low.x .= tracker.high.x #this automatically converts to Low
    tracker.low.xnext .= tracker.low.x

    tracker.usehigh = true
end
