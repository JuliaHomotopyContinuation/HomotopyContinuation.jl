export AbstractPathTracker,
    Tracker,
    TrackerResult,
    TrackerOptions,
    TrackerParameters,
    TrackerCode,
    DEFAULT_TRACKER_PARAMETERS,
    FAST_TRACKER_PARAMETERS,
    CONSERVATIVE_TRACKER_PARAMETERS,
    track,
    init!,
    track!,
    step!,
    start_parameters!,
    target_parameters!,
    status,
    state,
    is_success,
    is_terminated,
    is_invalid_startvalue,
    is_tracking,
    solution,
    steps,
    accepted_steps,
    rejected_steps,
    iterator


include("predictor.jl")
include("newton_corrector.jl")
include("tracker_options.jl")
include("tracker_result.jl")


Base.@kwdef mutable struct TrackerStatistics
    accepted_steps::Int = 0
    rejected_steps::Int = 0
    last_steps_accepted::Int = 0
    last_steps_failed::Int = 0
    ext_accepted_steps::Int = 0
    ext_rejected_steps::Int = 0
    failed_attempts_to_go_to_zero::Int = 0
end
Base.show(io::IO, S::TrackerStatistics) = print_fieldnames(io, S)
function init!(S::TrackerStatistics; keep_steps::Bool = false)
    if !keep_steps
        S.accepted_steps = 0
        S.rejected_steps = 0
        S.ext_accepted_steps = 0
        S.ext_rejected_steps = 0
        S.failed_attempts_to_go_to_zero = 0
    end

    S.last_steps_accepted = 0
    S.last_steps_failed = 0
end

steps(S::TrackerStatistics) = S.accepted_steps + S.rejected_steps
ext_steps(S::TrackerStatistics) = S.ext_accepted_steps + S.ext_rejected_steps

###
### STATE
###

mutable struct MovingAverageVector
    v::Vector{Float64}
    n::Int
    i::Int
    sum::Float64
end
function MovingAverageVector(n::Int)
    v = ones(Float64, n)
    MovingAverageVector(v, n, 1, n)
end

function Base.push!(M::MovingAverageVector, x::Float64)
    M.sum -= M.v[M.i]
    M.v[M.i] = x
    M.sum += x
    M.i = (M.i % M.n) + 1
end

mean_value(M::MovingAverageVector) = M.sum / M.n
max_value(M::MovingAverageVector) = maximum(M.v)

function init!(M::MovingAverageVector)
    fill!(M.v, 1.0)
    M.i = 1
    M.sum = M.n
end


Base.@kwdef mutable struct TrackerState
    code::TrackerCode.codes = TrackerCode.tracking
    x::Vector{ComplexF64} # current x
    x̂::Vector{ComplexF64} # last prediction
    x̄::Vector{ComplexF64} # candidate for new x
    # internal step size
    segment_stepper::SegmentStepper = SegmentStepper(1.0, 0.0)
    Δs_prev::Float64 = NaN # previous step size 
    # path tracking algorithm
    accuracy::Float64 = NaN # norm(x - x(t))
    ω::Float64 = NaN # liptschitz constant estimate, see arxiv:1902.02968
    ω_prev::Float64 = NaN
    μ::Float64 = NaN # limit accuracy
    τ::Float64 = NaN # trust region size
    norm_Δx₀::Float64 = NaN# debug info only
    norm_Δx₀_prev::Float64 = NaN# debug info only
    extended_prec::Bool = false
    used_extended_prec::Bool = false
    refined_extended_prec::Bool = false
    keep_extended_prec::Bool = false
    norm::WeightedNorm{InfNorm}
    use_strict_β_τ::Bool = false
    β_trust_region::Float64 = 0.8
    act_err_pred_err_ratio::MovingAverageVector = MovingAverageVector(6)
    statistics::TrackerStatistics = TrackerStatistics()
end

function TrackerState(H; norm::WeightedNorm{InfNorm})
    m, n = size(H)
    x = zeros(ComplexF64, m)
    x̂ = zero(x)
    x̄ = zero(x)

    TrackerState(x = x, x̂ = x̂, x̄ = x̄, norm = norm)
end

Base.show(io::IO, state::TrackerState) = print_fieldnames(io, state)
function Base.getproperty(state::TrackerState, sym::Symbol)
    if sym === :t
        return getfield(state, :segment_stepper).t
    elseif sym == :Δt
        return getfield(state, :segment_stepper).Δt
    elseif sym == :t′
        return getfield(state, :segment_stepper).t′
    elseif sym == :last_step_failed
        return getfield(state.statistics, :last_steps_failed) > 0
    else # fallback to getfield
        return getfield(state, sym)
    end
end

function init!(
    state::TrackerState,
    H::AbstractHomotopy,
    x::AbstractVector,
    t₁,
    t₀;
    keep_steps::Bool = false,
)
    set_solution!(state.x, H, x, t₁)
    state.code = TrackerCode.tracking
    init!(state.segment_stepper, t₁, t₀)
    state.Δs_prev = NaN
    state.accuracy = NaN
    state.ω = NaN
    state.ω_prev = NaN
    state.μ = NaN
    state.τ = NaN
    state.norm_Δx₀ = state.norm_Δx₀_prev = NaN
    state.extended_prec = false
    state.used_extended_prec = false
    state.refined_extended_prec = false
    state.keep_extended_prec = false
    state.use_strict_β_τ = false
    state.β_trust_region = 0.8
    init!(state.statistics; keep_steps = keep_steps)
    init!(state.act_err_pred_err_ratio)
    state
end

function current_tol(state::TrackerState)
    μ_ = 2.842170943040401e-14 # 128eps()
    ω = isnan(state.ω) ? 1.0 : state.ω
    if (state.extended_prec)
        return √(μ_ / ω)
    else
        return 1e-5 / ω
    end
end


###
### TRACKER
###

"""
    Tracker(H::AbstractHomotopy;
            options = TrackerOptions(),
            weighted_norm_options = WeightedNormOptions())

Construct a tracker for the given homotopy `H`. The algorithm computes along the path ``x(t)``
the local derivatives up to order 4.
For `options` see also [`TrackerOptions`](@ref).
The algorithm uses as a weighted infinity norm to measure distances.
See also [`WeightedNorm`](@ref).

[^Tim20]: Timme, S. "Mixed Precision Path Tracking for Polynomial Homotopy Continuation". arXiv:1902.02968 (2020)


## Example

We want to solve the system
```julia
@var x y t
F = System([x^2 + y^2 - 3, 2x^2 + 0.5x*y + 3y^2 - 2])
```
using a total degree homotopy and `Tracker`.
```julia
# construct start system and homotopy
G = System(im * [x^2 - 1, y^2 - 1])
H = StraightLineHomotopy(G, F)
start_solutions = [[1,1], [-1,1], [1,-1], [-1,-1]]
# construct tracker
tracker = Tracker(H)
# track each start solution separetely
results = track.(tracker, start_solutions)
println("# successfull: ", count(is_success, results))
```
We see that we tracked all 4 paths successfully.
```
# successfull: 4
```
"""
struct Tracker{H<:AbstractHomotopy}
    homotopy::H
    ws::LinearSolveWorkspace
    predictor::Predictor
    corrector::NewtonCorrector
    # these are mutable
    state::TrackerState
    options::TrackerOptions
end

Tracker(H::ModelKit.Homotopy; compile::Union{Bool,Symbol} = COMPILE_DEFAULT[], kwargs...) =
    Tracker(fixed(H; compile = compile); kwargs...)

function Tracker(
    H::AbstractHomotopy;
    weighted_norm_options::WeightedNormOptions = WeightedNormOptions(),
    options = TrackerOptions(),
)
    if !isa(options, TrackerOptions)
        options = TrackerOptions(; options...)
    else
        options = deepcopy(options)
    end
    m, n = size(H)
    norm = WeightedNorm(ones(n), InfNorm(), weighted_norm_options)

    ws = LinearSolveWorkspace(m, n, norm)
    state = TrackerState(H; norm = norm)
    predictor = Predictor(H, ws, norm; max_order = options.parameters.predictor_local_order)
    corrector = NewtonCorrector(ws, norm)

    Tracker(H, ws, predictor, corrector, state, options)
end

Base.show(io::IO, C::Tracker) = print(io, typeof(C), "()")
Base.broadcastable(C::Tracker) = Ref(C)

"""
    state(tracker::Tracker)

Return the state of the tracker.
"""
state(tracker::Tracker) = tracker.state

"""
    status(tracker::Tracker)

Get the current [`TrackerCode`](@ref) of `tracker`.
"""

"""
    solution(tracker::Tracker)

Get the current solution.
"""
solution(T::Tracker) = get_solution(T.homotopy, T.state.x, T.state.t)

status(tracker::Tracker) = tracker.state.code

is_done(tracker::Tracker) = is_done(tracker.state.segment_stepper)

#############
# TRACKING ##
#############

"""
    track!(tracker::Tracker, x, t₁ = 1.0, t₀ = 0.0; debug::Bool = false)

The same as [`track`](@ref) but only returns the final [`TrackerCode`](@ref).

    track!(tracker::Tracker, t₀; debug::Bool = false)

Track with `tracker` the current solution to `t₀`. This keeps the current state.
"""
function track!(
    tracker::Tracker,
    x::AbstractVector,
    t₁ = 1.0,
    t₀ = 0.0;
    keep_steps::Bool = false,
)
    init!(tracker, x, t₁, t₀; keep_steps = keep_steps)

    while is_tracking(tracker.state.code)
        step!(tracker)
    end

    tracker.state.code
end

function track!(tracker::Tracker, r::TrackerResult, t₁ = 1.0, t₀ = 0.0; debug::Bool = false)
    track!(
        tracker,
        solution(r),
        t₁,
        t₀;
        ω = r.ω,
        μ = r.μ,
        τ = r.τ,
        extended_precision = r.extended_precision,
    )
end
function track!(tracker::Tracker, t₀; max_initial_step_size::Float64 = Inf)
    init!(tracker, t₀; max_initial_step_size = max_initial_step_size)

    while is_tracking(tracker.state.code)
        step!(tracker)
    end

    tracker.state.code
end


##########
## INIT ##
##########

"""
    init!(tracker::Tracker, x₁, t₁, t₀)

Setup `tracker` to track `x₁` from `t₁` to `t₀`.

    init!(tracker::Tracker, t₀)

Setup `tracker` to continue tracking the current solution to `t₀`.
This keeps the current state.
"""
init!(tracker::Tracker, r::TrackerResult, t₁::Number = 1.0, t₀::Number = 0.0) =
    init!(tracker, solution(r), t₁, t₀; ω = r.ω, μ = r.μ)

function init!(
    tracker::Tracker,
    x₁::AbstractVector,
    t₁::Number = 1.0,
    t₀::Number = 0.0;
    max_initial_step_size::Float64 = Inf,
    keep_steps::Bool = false,
)
    @unpack state, predictor, corrector, homotopy, options = tracker
    @unpack x, x̄, norm = state

    init!(state, homotopy, x₁, t₁, t₀; keep_steps = keep_steps)
    init!(norm, x)
    init!(tracker.ws)
    is_valid = check_valid_startsolution!(tracker)
    if !is_valid
        return
    end


    # initialize the predictor
    init!(tracker.predictor)

    localize!(predictor, homotopy, state.x, state.t)
    state.τ = trust_region(predictor)

    # Deal with the case that the initial solution is already a solution of the homotopy.
    # It can still happen that the solution is singular at t=0 so we cannot really on the corrector step working

    # This condition is only the case if all derivatives are 0
    if !isfinite(state.τ) && predictor.order[] == 1
        state.code = TrackerCode.success
        return
    end


    # compute initial step size

    Δs′ = initial_step_size(state, predictor, tracker.options)
    Δs = clamp(Δs′, tracker.options.min_step_size, max_initial_step_size)
    propose_step!(state.segment_stepper, Δs)


    is_tracking(state.code)
end


function check_valid_startsolution!(tracker)
    @unpack state, corrector, homotopy, ws = tracker
    @unpack x, x̄, t = state
    for extended_precision in (false, true)
        res = newton!(
            x̄,
            corrector,
            homotopy,
            x,
            t;
            max_iters = 3,
            tol = 1e-10,
            extended_precision = extended_precision,
        )
        if is_converged(res)
            state.accuracy = res.accuracy
            state.μ = max(res.μ_low, eps())
            x .= x̄
            return true
        end
    end

    state.code = if rcond!(ws) < 1e-12
        TrackerCode.terminated_invalid_startvalue_singular_jacobian
    else
        TrackerCode.terminated_invalid_startvalue
    end

    return false
end


# function init!(tracker::Tracker, t₀::Number; max_initial_step_size::Float64=Inf)
#   @unpack state, predictor, options = Pathtracker
#   state.code = TrackerCode.tracking
#   init!(state.segment_stepper, state.t, t₀)
#   Δs = initial_step_size(state, predictor, tracker.options)
#   Δs = min(Δs, max_initial_step_size)
#   propose_step!(state.segment_stepper, Δs)
#   state.Δs_prev = 0.0

#   tracker
# end 

###############
## Step Size ##
###############


# intial step size
function initial_step_size(
    state::TrackerState,
    predictor::Predictor,
    options::TrackerOptions,
)
    τ = trust_region(predictor)
    Δs₁ = 0.25τ
    min(Δs₁, options.max_step_size, options.max_initial_step_size)
end

function update_stepsize!(
    state::TrackerState,
    result::NewtonCorrectorResult,
    options::TrackerOptions,
    predictor::Predictor;
)
    @unpack ω, Δs_prev = state
    @unpack N = options.parameters

    p = predictor.order[]

    Δs = state.segment_stepper.Δs
    act_err = state.norm_Δx₀
    η = state.norm_Δx₀_prev / abs(state.Δs_prev)^p
    pred_err = η * abs(Δs)^p
    r = act_err / pred_err
    if isnan(r)
        r = 1.0
    end
    push!(state.act_err_pred_err_ratio, r)


    if is_converged(result)
        m = 2^(N - 1)

        α = max(max_value(state.act_err_pred_err_ratio), 1)

        ε = current_tol(state)
        δ1 = nthroot(ε * ω^(-m + 1), m)
        δ2 = 1 / (m * ω)
        Δs₁ = nthroot(min(δ1, δ2) / (α * result.norm_Δx₀), p) * abs(Δs)
        Δs₂ = 0.5trust_region(predictor)
        Δs′ = min(Δs₁, Δs₂)
        # increase step size by a factor of at most 8 in one step
        if state.last_step_failed
            Δs′ = min(Δs′, 1.25Δs)
        else
            Δs′ = min(Δs′, 8Δs)
        end
        step_success!(state.segment_stepper)
    else
        Δs′ = 0.5 * Δs
    end

    state.Δs_prev = Δs
    propose_step!(state.segment_stepper, Δs′)
    nothing
end


function check_terminated!(state::TrackerState, options::TrackerOptions)
    if state.extended_prec || !options.extended_precision
        @unpack N, target_accuracy = options.parameters
        @unpack ω = state
        m = 2^N

        μ_threshold = ((nthroot(target_accuracy, m) / ω^((m - 1) / m))^2 * ω)
        if state.μ > μ_threshold && state.statistics.last_steps_failed > 3
            state.code = TrackerCode.terminated_accuracy_limit
            return
        end
    end

    t′ = state.segment_stepper.t′
    t = state.segment_stepper.t
    if is_done(state.segment_stepper)
        state.code = TrackerCode.success
    elseif steps(state.statistics) ≥ options.max_steps
        state.code = TrackerCode.terminated_max_steps
    elseif state.segment_stepper.Δs < options.min_step_size
        state.code = TrackerCode.terminated_step_size_too_small
    elseif fast_abs(t′ - t) ≤ 2eps(fast_abs(t))
        state.code = TrackerCode.terminated_step_size_too_small
    elseif options.min_rel_step_size > 0 &&
           t′ != state.segment_stepper.target &&
           fast_abs(t′ - t) < fast_abs(t) * options.min_rel_step_size
        state.code = TrackerCode.terminated_step_size_too_small
    end
    nothing
end

function update_predictor!(tracker::Tracker)
    @unpack predictor, homotopy, state = tracker
    localize!(predictor, homotopy, state.x, state.t)
end


function update_precision!(tracker::Tracker, μ_low)
    @unpack homotopy, corrector, state, options = tracker
    @unpack μ, ω, x, t, norm = state
    @unpack N = options.parameters

    options.extended_precision || return false

    β = 0.0625
    μ_threshold = β * current_tol(state)
    if state.extended_prec && !state.keep_extended_prec && !isnan(μ_low) && μ_low > μ
        # check if we can go low again
        if μ_low < 0.125 * μ_threshold
            state.extended_prec = false
            state.μ = μ_low
        end
    elseif μ > μ_threshold
        use_extended_precision!(tracker)
    end

    state.extended_prec
end

function use_extended_precision!(tracker::Tracker)
    @unpack homotopy, corrector, state, options = tracker
    @unpack μ, ω, x, t, norm = state

    options.extended_precision || return state.μ
    !state.extended_prec || return state.μ

    state.extended_prec = true
    state.used_extended_prec = true

    newton_res = newton!(
        x,
        corrector,
        homotopy,
        x,
        t;
        extended_precision = true,
        tol = 1e-14,
        max_iters = 3,
    )
    state.μ = max(newton_res.accuracy, eps())
    state.μ
end

function refine_current_solution!(
    tracker,
    x = tracker.state.x,
    t = tracker.state.t;
    min_tol::Float64 = 4 * eps(),
)
    @unpack homotopy, corrector = tracker

    res = newton!(
        x,
        corrector,
        homotopy,
        x,
        t;
        extended_precision = true,
        tol = min_tol,
        max_iters = 3,
    )
    res.accuracy
end

"""
    step!(tracker::Tracker, debug::Bool = false)

Perform a single tracking step. Returns `true` if the step was accepted.
"""
function step!(tracker::Tracker)
    @unpack homotopy, corrector, predictor, state, options = tracker
    @unpack t, Δt, t′, x, x̄, norm = state

    # Use the current approximation of x(t) to obtain estimate
    # x̂ ≈ x(t + Δt) using the predictor
    predict!(predictor, Δt)
    x̂ = prediction(predictor)

    # Correct the predicted value x̂ to obtain x̄.
    # If the correction is successfull we have x̄ ≈ x(t+Δt).

    result = newton!(
        x̄,
        corrector,
        homotopy,
        x̂,
        t′; # = t + Δt;
        max_iters = options.parameters.N,
        tol = current_tol(state),
        extended_precision = state.extended_prec,
    )

    state.norm_Δx₀_prev = state.norm_Δx₀
    state.norm_Δx₀ = result.norm_Δx₀
    if is_converged(result)

        state.statistics.accepted_steps += 1
        state.statistics.ext_accepted_steps += state.extended_prec
        state.statistics.last_steps_failed = 0
        state.statistics.last_steps_accepted += 1
        state.accuracy = result.accuracy
        state.μ = max(result.accuracy, eps())
        ω = isnan(result.ω) ? state.ω : max(result.ω, 1)
        if isnan(ω)
            ω = 1
        end

        state.ω_prev = state.ω
        state.ω = ω
    else
        state.statistics.failed_attempts_to_go_to_zero += iszero(t′)
        # Step failed, so we have to try with a  (smaller) step size
        state.statistics.rejected_steps += 1
        state.statistics.ext_rejected_steps += state.extended_prec
        state.statistics.last_steps_failed += 1
        state.statistics.last_steps_accepted = 0
    end

    update_stepsize!(state, result, options, predictor)

    if is_converged(result)
        x .= x̄
        update_precision!(tracker, result.μ_low)
        if is_done(state.segment_stepper) &&
           options.extended_precision &&
           state.accuracy > 1e-12
            state.accuracy = refine_current_solution!(tracker; min_tol = 1e-14)
            state.refined_extended_prec = true
        end

        localize!(predictor, homotopy, state.x, state.t)
        state.τ = trust_region(predictor)

        update!(state.norm, x)
    end

    check_terminated!(state, options)

    !state.last_step_failed
end


function TrackerResult(H::AbstractHomotopy, state::TrackerState)
    TrackerResult(
        Symbol(state.code),
        get_solution(H, state.x, state.t),
        state.t,
        state.accuracy,
        state.statistics.accepted_steps,
        state.statistics.rejected_steps,
        state.extended_prec || state.refined_extended_prec,
        state.used_extended_prec || state.refined_extended_prec,
        state.ω,
        state.μ,
        state.τ,
    )
end

"""
     track(tracker::Tracker, x::AbstractVector, t₁ = 1.0, t₀ = 0.0; debug::Bool = false)

Track the given solution `x` at `t₁` using `tracker` to a solution at `t₀`.

    track(tracker::Tracker, r::TrackerResult, t₁ = 1.0, t₀ = 0.0; debug::Bool = false)

Track the solution of the result `r` from `t₁` to `t₀`.
"""
@inline function track(tracker::Tracker, x, t₁ = 1.0, t₀ = 0.0; kwargs...)
    track!(tracker, x, t₁, t₀; kwargs...)
    TrackerResult(tracker.homotopy, tracker.state)
end

"""
    start_parameters!(tracker::Tracker, p)

Set the start parameters of the homotopy of the tracker.
"""
start_parameters!(T::Tracker, p) = (start_parameters!(T.homotopy, p); T)

"""
    target_parameters!(tracker::Tracker, p)

Set the target parameters of the homotopy of the tracker.
"""
target_parameters!(T::Tracker, p) = (target_parameters!(T.homotopy, p); T)

parameters!(T::Tracker, p, q) = (parameters!(T.homotopy, p, q); T)

# PathIterator #
struct PathIterator{T<:Tracker}
    tracker::T
    t_real::Bool
end
Base.IteratorSize(::Type{<:PathIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:PathIterator}) = Base.HasEltype()

"""
    iterator(tracker::Tracker, x₁, t₁=1.0, t₀=0.0)

Prepare a tracker to make it usable as a (stateful) iterator. Use this if you want to inspect a specific
path. In each iteration the tuple `(x,t)` is returned.

## Example

Assume you have `Tracker` `tracker` and you wan to track `x₁` from 1.0 to 0.25:
```julia
for (x,t) in iterator(tracker, x₁, 1.0, 0.25)
    println("x at t=\$t:")
    println(x)
end
```

Note that this is a stateful iterator. You can still introspect the state of the tracker.
For example to check whether the tracker was successfull
(and did not terminate early due to some problem) you can do
```julia
println("Success: ", is_success(status(tracker)))
```
"""
function iterator(tracker::Tracker, x₁, t₁ = 1.0, t₀ = 0.0; kwargs...)
    init!(tracker, x₁, t₁, t₀; kwargs...)
    PathIterator(tracker, typeof(t₁ - t₀) <: Real)
end

function current_x_t(iter::PathIterator)
    @unpack x, t = iter.tracker.state
    (copy(x), iter.t_real ? real(t) : t)
end

function Base.iterate(iter::PathIterator, state = nothing)
    state === nothing && return current_x_t(iter), 1
    iter.tracker.state.code != TrackerCode.tracking && return nothing

    while is_tracking(iter.tracker.state.code)
        step_failed = !step!(iter.tracker)
        step_failed || break
    end
    current_x_t(iter), state + 1
end


include("path_info.jl")
