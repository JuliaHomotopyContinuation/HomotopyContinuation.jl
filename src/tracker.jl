export Tracker,
    TrackerResult,
    track,
    track!,
    is_success,
    solution,
    steps,
    accepted_steps,
    rejected_steps,
    iterator

"""
    TrackerOptions

The set of options set for a [`Tracker`](@ref). See the description of [`Tracker`](@ref)
for all possible options.
"""
Base.@kwdef mutable struct TrackerOptions
    max_steps::Int = 1_000
    max_step_size::Float64 = Inf
    a::Float64 = 0.125
    β_a::Float64 = 1.0
    β_ω::Float64 = 10.0
    β_τ::Float64 = 0.75
    extended_precision::Bool = true
    min_step_size::Float64 = exp2(-80)
end

Base.show(io::IO, opts::TrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::TrackerOptions) =
    opts

module TrackerReturnCode

@doc """
    TrackerReturnCode.codes

The possible states a `CoreTracker` can have are

* `TrackerReturnCode.tracking`: The tracking is still in progress
* `TrackerReturnCode.success`: Indicates a successfull completed tracking
* `TrackerReturnCode.terminated_max_iters`: Tracking terminated since maximal iterations reached.
* `TrackerReturnCode.terminated_ill_conditioned`: Tracking terminated since the path was too ill-conditioned.
* `TrackerReturnCode.terminated_invalid_startvalue`: Tracking terminated since the provided start value was invalid.
* `TrackerReturnCode.terminated_step_size_too_small`
"""
@enum codes begin
    tracking
    success
    terminated_max_iters
    terminated_accuracy_limit
    terminated_ill_conditioned
    terminated_invalid_startvalue
    terminated_step_size_too_small
    terminated_unknown
end

end

"""
    is_success(S::TrackerReturnCode.codes)

Returns `true` if `S` indicates a success in the path tracking.
"""
is_success(S) = S == TrackerReturnCode.success

"""
    is_terminated(S::TrackerReturnCode.codes)

Returns `true` if `S` indicates that the path tracking got terminated. This is not `true`
if `is_success(S)` is `true`.
"""
is_terminated(S::TrackerReturnCode.codes) =
    S ≠ TrackerReturnCode.tracking && S ≠ TrackerReturnCode.success

"""
    is_invalid_startvalue(S::TrackerReturnCode.codes)

Returns `true` if `S` indicates that the path tracking got terminated since the start
value was not a zero.
"""
is_invalid_startvalue(S::TrackerReturnCode.codes) =
    S == TrackerReturnCode.terminated_invalid_startvalue

"""
    is_tracking(S::TrackerReturnCode.codes)

Returns `true` if `S` indicates that the path tracking is not yet finished.
"""
is_tracking(S::TrackerReturnCode.codes) = S == TrackerReturnCode.tracking


# RESULT
"""
    TrackerResult

Containing the result of tracking a path with `Tracker`.
"""
struct TrackerResult{V<:AbstractVector{ComplexF64}}
    returncode::TrackerReturnCode.codes
    solution::V
    t::ComplexF64
    accuracy::Float64
    ω::Float64
    μ::Float64
    accepted_steps::Int
    rejected_steps::Int
    extended_precision_used::Bool
end

"""
    is_success(result::TrackerResult)

Returns `true` if the path tracking was successfull.
"""
is_success(result::TrackerResult) = is_success(result.returncode)

"""
    solution(result::TrackerResult)

Returns the solutions obtained by the `Tracker`.
"""
solution(result::TrackerResult) = result.solution

"""
    steps(result::TrackerResult)

Returns the number of steps done.
"""
steps(result::TrackerResult) = result.accepted_steps + result.rejected_steps

"""
    accepted_steps(result::TrackerResult)

Returns the number of accepted steps.
"""
accepted_steps(result::TrackerResult) = result.accepted_steps

"""
    rejected_steps(result::TrackerResult)

Returns the number of rejected_steps steps.
"""
rejected_steps(result::TrackerResult) = result.rejected_steps

Base.show(io::IO, result::TrackerResult) = print_fieldnames(io, result)
Base.show(io::IO, ::MIME"application/prs.juno.inline", result::TrackerResult) =
    result

mutable struct TrackerState{
    V<:AbstractVector{ComplexF64},
    M<:AbstractMatrix{ComplexF64},
}
    x::V # current x
    x̂::V # last prediction
    x̄::V # candidate for new x
    # internal step size
    segment::ComplexLineSegment
    s::DoubleF64 # current step length (0 ≤ s ≤ length(segment))
    s′::DoubleF64 # proposed s (0 ≤ s ≤ length(segment))
    Δs_prev::Float64 # previous step size
    # path tracking algorithm
    accuracy::Float64 # norm(x - x(t))
    ω::Float64 # liptschitz constant estimate, see arxiv:1902.02968
    μ::Float64 # limit accuracy
    τ::Float64 # trust region size
    norm_Δx₀::Float64 # debug info only
    use_extended_prec::Bool
    used_extended_prec::Bool
    norm::WeightedNorm{InfNorm}

    jacobian::Jacobian{M}
    cond_J_ẋ::Float64 # estimate of cond(H(x(t),t), ẋ(t))
    code::TrackerReturnCode.codes
    # first four derivatives of x(t)
    x¹::Vector{ComplexF64}
    x²::Vector{ComplexF64}
    x³::Vector{ComplexF64}
    x⁴::Vector{ComplexF64}
    # these are just tuples containing the derivative vectors
    # we have to preallocate them since otherwise this would result
    # in dynamic dispatch (at least Julia ≤ 1.3)
    dx¹::NTuple{1,Vector{ComplexF64}}
    dx²::NTuple{2,Vector{ComplexF64}}
    dx³::NTuple{3,Vector{ComplexF64}}
    # buffer for the computation of derivatives
    u::Vector{ComplexF64}

    # statistics
    accepted_steps::Int
    rejected_steps::Int
    last_step_failed::Bool
end

function TrackerState(
    H,
    x₁::AbstractVector,
    t₁,
    t₀,
    norm::WeightedNorm{InfNorm},
)
    x = isa(x₁, PVector) ? ComplexF64.(x₁) : Vector{ComplexF64}(x₁)
    x̂ = zero(x)
    x̄ = zero(x)
    x¹, x², x³, x⁴ = [zeros(ComplexF64, length(x)) for _ = 1:4]
    dx¹ = (x¹,)
    dx² = (x¹, x²)
    dx³ = (x¹, x², x³)
    segment = ComplexLineSegment(t₁, t₀)
    s = s′ = length(segment)
    Δs_prev = 0.0
    accuracy = 0.0
    μ = eps()
    ω = 1.0
    τ = Inf
    norm_Δx₀ = NaN
    used_extended_prec = use_extended_prec = false
    jacobian = Jacobian(zeros(ComplexF64, size(H)))
    cond_J_ẋ = NaN
    code = TrackerReturnCode.tracking
    u = zeros(ComplexF64, size(H, 1))
    accepted_steps = rejected_steps = 0
    last_step_failed = true

    TrackerState(
        x,
        x̂,
        x̄,
        segment,
        s,
        s′,
        Δs_prev,
        accuracy,
        ω,
        μ,
        τ,
        norm_Δx₀,
        use_extended_prec,
        used_extended_prec,
        norm,
        jacobian,
        cond_J_ẋ,
        code,
        x¹,
        x²,
        x³,
        x⁴,
        dx¹,
        dx²,
        dx³,
        u,
        accepted_steps,
        rejected_steps,
        last_step_failed,
    )
end

Base.show(io::IO, state::TrackerState) = print_fieldnames(io, state)
Base.show(io::IO, ::MIME"application/prs.juno.inline", state::TrackerState) =
    state
function Base.getproperty(state::TrackerState, sym::Symbol)
    if sym === :t
        return getfield(state, :segment)[getfield(state, :s)]
    elseif sym == :Δt
        segment = getfield(state, :segment)
        return step_size(segment, getfield(state, :s′) - getfield(state, :s))
    elseif sym == :t′
        return getfield(state, :segment)[getfield(state, :s′)]
    else # fallback to getfield
        return getfield(state, sym)
    end
end

steps(S::TrackerState) = S.accepted_steps + S.rejected_steps


function compute_derivatives!(
    state::TrackerState,
    H::AbstractHomotopy,
    x,
    t;
    log_scale::Bool = false,
)
    # unpack stuff to make the rest easier to read
    @unpack u, x¹, x², x³, x⁴, dx¹, dx², dx³, jacobian, norm = state

    # compute all taylor coefficients x¹, x², x³, x⁴
    diff_t!(u, H, x, t)
    u .= .-u
    LA.ldiv!(x¹, jacobian, u)

    # Check if we have to do iterative refinment for all the others as well
    δ = iterative_refinement!(x¹, jacobian, u, norm; fixed_precision = true)
    state.cond_J_ẋ = δ / eps()
    iterative_refinement = δ > sqrt(eps())
    if iterative_refinement
        if !derivative_refinement!(x¹, jacobian, u, norm)
            state.code = TrackerReturnCode.terminated_ill_conditioned
        end
    end

    diff_t!(u, H, x, t, dx¹)
    u .= .-u
    LA.ldiv!(x², jacobian, u)
    iterative_refinement && derivative_refinement!(x², jacobian, u, norm)

    diff_t!(u, H, x, t, dx²)
    u .= .-u
    LA.ldiv!(x³, jacobian, u)
    iterative_refinement && derivative_refinement!(x³, jacobian, u, norm)

    diff_t!(u, H, x, t, dx³)
    u .= .-u
    LA.ldiv!(x⁴, jacobian, u)
    iterative_refinement && derivative_refinement!(x⁴, jacobian, u, norm)

    state
end

function derivative_refinement!(x, jacobian, u, norm)
    δ̂ = iterative_refinement!(
        x,
        jacobian,
        u,
        norm;
        tol = sqrt(eps()),
        max_iters = 3,
    )
    δ̂ < sqrt(eps())
end

# TRACKER

struct Tracker{
    H<:AbstractHomotopy,
    # V and V̄ need to have the same container type
    V<:AbstractVector{ComplexF64},
    V̄<:AbstractVector{ComplexDF64},
    M<:AbstractMatrix{ComplexF64},
}
    homotopy::H
    predictor::Pade21
    corrector::NewtonCorrector{V̄}
    # these are mutable
    state::TrackerState{V,M}
    options::TrackerOptions
end

function Tracker(H::AbstractHomotopy; kwargs...)
    x = zeros(ComplexF64, size(H, 2))
    Tracker(H, x, 1.0, 0.0; kwargs...)
end
function Tracker(H::AffineChartHomotopy; kwargs...)
    d = dims(H.chart)
    x = PVector(randn(ComplexF64, sum(d) + length(d)), d)
    Tracker(H, x, 1.0, 0.0; kwargs...)
end

function Tracker(
    H::AbstractHomotopy,
    x₁::AbstractVector,
    t₁::Number,
    t₀::Number;
    norm::WeightedNorm{InfNorm} = WeightedNorm(InfNorm(), size(H, 2)),
    kwargs...,
)
    options = TrackerOptions(; kwargs...)
    state = TrackerState(H, x₁, t₁, t₀, norm)
    predictor = Pade21(size(H, 2))
    corrector = NewtonCorrector(options.a, state.x, size(H, 1))

    Tracker(H, predictor, corrector, state, options)
end

Base.show(io::IO, C::Tracker) = print(io, "Tracker")
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::Tracker) = x
Base.broadcastable(C::Tracker) = Ref(C)

"""
    state(tracker::Tracker)

Return the state of the tracker.
"""
state(tracker::Tracker) = tracker.state
status(tracker::Tracker) = tracker.state.code
LA.cond(tracker::Tracker) = tracker.state.cond_J_ẋ
# Step Size

_h(a) = 2a * (√(4 * a^2 + 1) - 2a)
# intial step size
function initial_step_size(
    state::TrackerState,
    predictor::Pade21,
    options::TrackerOptions,
)
    a = options.β_a * options.a
    ω = options.β_ω * state.ω
    e = state.norm(local_error(predictor))
    Δs₁ = nthroot((√(1 + 2 * _h(a)) - 1) / (ω * e), order(predictor))
    Δs₂ = options.β_τ * trust_region(predictor)
    min(Δs₁, Δs₂, options.max_step_size)
end

function update_stepsize!(
    state::TrackerState,
    result::NewtonCorrectorResult,
    options::TrackerOptions,
    predictor::Pade21,
)
    a = options.β_a * options.a
    ω = options.β_ω * state.ω
    p = order(predictor)
    if is_converged(result)
        e = state.norm(local_error(predictor))
        Δs₁ = nthroot((√(1 + 2 * _h(a)) - 1) / (ω * e), p)
        Δs₂ = options.β_τ * trust_region(predictor)
        Δs = min(Δs₁, Δs₂)
        # In the double double arithmetic implementation 1 + Inf = NaN, so we have
        # to be careful here to not add Inf
        s′ = DoubleF64(0.0)
        if !isinf(Δs)
            s′ = max(s′, state.s - Δs)
        end
        if !isinf(options.max_step_size)
            s′ = max(s′, state.s - options.max_step_size)
        end
        if state.last_step_failed
            s′ = max(state.s + state.Δs_prev, s′)
        end
        state.s′ = s′
    else
        j = result.iters - 2
        Θ_j = nthroot(result.θ, 1 << j)
        state.s′ =
            state.s +
            nthroot((√(1 + 2 * _h(0.5a)) - 1) / (√(1 + 2 * _h(Θ_j)) - 1), p) *
            (state.s′ - state.s)
    end
    nothing
end


function check_terminated!(state::TrackerState, options::TrackerOptions)
    if state.s == 0.0
        state.code = TrackerReturnCode.success
    elseif steps(state) ≥ options.max_steps
        state.code = TrackerReturnCode.terminated_max_iters
    elseif state.ω * state.μ > options.a * _h(options.a)
        state.code = TrackerReturnCode.terminated_accuracy_limit
    elseif state.s′ == state.s || fast_abs(state.Δt) < eps(fast_abs(state.t)) ||
           (state.s - state.s′) < options.min_step_size
        state.code = TrackerReturnCode.terminated_step_size_too_small
    elseif isnan(state.s′) # catch any NaNs produced somewhere
        state.code = TrackerReturnCode.terminated_unknown
    end
    nothing
end

function compute_derivatives_and_update_predictor!(tracker::Tracker)
    @unpack x, t, x¹, x², x³, x⁴ = tracker.state
    compute_derivatives!(tracker.state, tracker.homotopy, x, t)
    update!(tracker.predictor, x, x¹, x², x³, x⁴)
end

"""
    init!(tracker::NCT, x₁, t₁, t₀)

Setup `tracker` to track `x₁` from `t₁` to `t₀`.

    init!(tracker::NCT, t₀)

Setup `tracker` to continue tracking the current solution to `t₀`.
"""
init!(tracker::Tracker, r::TrackerResult, t₁, t₀; continuation::Bool = false) =
    init!(tracker, r.x, t₁, t₀, r.ω, r.μ; continuation = continuation)

function init!(
    tracker::Tracker,
    x₁::AbstractVector,
    t₁,
    t₀,
    ω::Float64 = NaN,
    μ::Float64 = NaN;
    keep_steps::Bool = false,
)
    @unpack state, predictor, corrector, homotopy, options = tracker
    @unpack x, x̄, norm, jacobian = state

    # intialize state
    x .= x₁
    state.segment = ComplexLineSegment(t₁, t₀)
    state.s = state.s′ = length(state.segment)
    state.Δs_prev = 0.0
    state.accuracy = eps()
    state.ω = 1.0
    state.used_extended_prec = state.use_extended_prec = false
    init!(norm, x)
    init!(jacobian)
    state.code = TrackerReturnCode.tracking
    if !keep_steps
        state.accepted_steps = state.rejected_steps = 0
    end
    state.last_step_failed = true

    # compute ω and limit accuracy μ for the start value
    t = state.t
    if isnan(ω) || isnan(μ)
        valid, ω, μ = init_newton!(
            x̄,
            corrector,
            homotopy,
            x,
            t,
            jacobian,
            norm;
            a = options.a,
        )
    else
        valid = true
    end
    if !isnan(ω)
        state.ω = ω
    end
    if valid
        state.accuracy = μ
        state.μ = max(μ, eps())
    end

    if !valid
        state.code = TrackerReturnCode.terminated_invalid_startvalue
    end

    # initialize the predictor
    evaluate_and_jacobian!(corrector.r, jacobian.J, homotopy, state.x, t)
    compute_derivatives_and_update_predictor!(tracker)
    state.τ = trust_region(predictor)
    # compute initial step size
    state.s′ =
        max(0.0, state.s - initial_step_size(state, predictor, tracker.options))
    tracker
end

function init!(tracker::Tracker, t₀)
    @unpack state, predictor, options = tracker
    state.segment = ComplexLineSegment(state.t, t₀)
    state.s = state.s′ = length(state.segment)
    state.Δs_prev = 0.0
    state.code = TrackerReturnCode.tracking
    state.s′ =
        max(0.0, state.s - initial_step_size(state, predictor, tracker.options))

    tracker
end

function update_precision!(tracker::Tracker, μ_low)
    @unpack homotopy, corrector, state, options = tracker
    @unpack μ, ω, x, t, jacobian, norm = state
    @unpack a = options

    options.extended_precision || return false

    if state.use_extended_prec && !isnan(μ_low)
        # check if we can go low again
        if μ_low * ω < a * (a^3)^2 * _h(a)
            state.use_extended_prec = false
            state.μ = μ_low
        end
    elseif μ * ω > a * (a^2)^2 * _h(a)
        state.use_extended_prec = true
        state.used_extended_prec = true
        # do one refinement step
        μ = extended_prec_refinement_step!(
            x,
            corrector,
            homotopy,
            x,
            t,
            jacobian,
            norm,
        )
        state.μ = max(μ, eps())
    end

    state.use_extended_prec
end

function step!(tracker::Tracker, debug::Bool = false)
    @unpack homotopy, corrector, predictor, state, options = tracker
    @unpack t, Δt, t′, x, x̂, x̄, jacobian, norm = state

    debug && printstyled(
        "\nt: ",
        round(t, sigdigits = 5),
        " Δt: ",
        round(Δt; sigdigits = 5),
        "\n";
        color = :yellow,
    )
    # Use the current approximation of x(t) to obtain estimate
    # x̂ ≈ x(t + Δt) using the predictor
    @unpack x¹, x², x³ = state
    predict!(x̂, predictor, homotopy, state.x, x¹, x², x³, Δt)

    update!(state.norm, x̂)

    # Correct the predicted value x̂ to obtain x̄.
    # If the correction is successfull we have x̄ ≈ x(t+Δt).
    result = newton!(
        x̄,
        corrector,
        homotopy,
        x̂,
        t′,
        jacobian,
        norm;
        ω = state.ω,
        μ = state.μ,
        extended_precision = state.use_extended_prec,
    )
    if debug
        printstyled(result, "\n"; color = is_converged(result) ? :green : :red)
    end
    if is_converged(result)
        # move forward
        x .= x̄
        state.Δs_prev = state.s′ - state.s
        state.s = state.s′
        state.accuracy = result.accuracy
        state.μ = max(result.accuracy, eps())
        state.ω = result.ω

        update_precision!(tracker, result.μ_low)
        compute_derivatives_and_update_predictor!(tracker)
        # Update other state
        state.τ = trust_region(predictor)
        state.accepted_steps += 1
    else
        # Step failed, so we have to try with a new (smaller) step size
        state.rejected_steps += 1
    end
    state.norm_Δx₀ = result.norm_Δx₀
    update_stepsize!(state, result, options, predictor)
    state.last_step_failed = !is_converged(result)

    check_terminated!(state, options)
    state.last_step_failed
end

function track!(
    tracker::Tracker,
    x::AbstractVector,
    t₁,
    t₀;
    ω::Float64 = NaN,
    μ::Float64 = NaN,
    keep_steps::Bool = false,
    debug::Bool = false,
)
    init!(tracker, x, t₁, t₀, ω, μ; keep_steps = keep_steps)

    while is_tracking(tracker.state.code)
        step!(tracker, debug)
    end

    (code = tracker.state.code, μ = tracker.state.μ, ω = tracker.state.ω)
end
function track!(tracker::Tracker, r::TrackerResult, t₁, t₀; debug::Bool = false)
    track!(tracker, solution(r), t₁, t₀; debug = debug, ω = r.ω, μ = r.μ)
end
function track!(tracker::Tracker, t₀; debug::Bool = false)
    init!(tracker, t₀)

    while is_tracking(tracker.state.code)
        step!(tracker, debug)
    end

    (code = tracker.state.code, μ = tracker.state.μ, ω = tracker.state.ω)
end

function TrackerResult(state::TrackerState)
    TrackerResult(
        state.code,
        copy(state.x),
        state.t,
        state.accuracy,
        state.ω,
        state.μ,
        state.accepted_steps,
        state.rejected_steps,
        state.used_extended_prec,
    )
end


@inline function track(tracker::Tracker, x, t₁, t₀; kwargs...)
    track!(tracker, x, t₁, t₀; kwargs...)
    TrackerResult(tracker.state)
end

# PathIterator #
struct PathIterator{H,P,M}
    tracker::Tracker{H,P,M}
    t_real::Bool
end
Base.IteratorSize(::Type{<:PathIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:PathIterator}) = Base.HasEltype()

"""
    iterator(tracker::CoreTracker, x₁, t₁=1.0, t₀=0.0)

Prepare a tracker to make it usable as a (stateful) iterator. Use this if you want to inspect a specific
path. In each iteration the tuple `(x,t)` is returned.

## Example

Assume you have `CoreTracker` `tracker` and you wan to track `x₁` from 1.0 to 0.25:
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
println("Success: ", status(tracker) == CoreTrackerStatus.success)
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
    iter.tracker.state.code != TrackerReturnCode.tracking && return nothing

    while is_tracking(iter.tracker.state.code)
        step_failed = step!(iter.tracker)
        step_failed || break
    end
    current_x_t(iter), state + 1
end
