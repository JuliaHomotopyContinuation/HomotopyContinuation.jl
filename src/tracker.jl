###########
# Options #
###########

"""
    TrackerOptions

The set of options set for a [`Tracker`](@ref). See the description of [`Tracker`](@ref)
for all possible options.
"""
Base.@kwdef struct TrackerOptions
    max_steps::Int = 1_000
    a::Float64 = 0.2
    β_a::Float64 = 0.1
    β_τ::Float64 = 0.7
end

Base.show(io::IO, opts::TrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::TrackerOptions) =
    opts


#########
# State #
#########

module TrackerCondition

@doc """
    TrackerCondition.conditions

The possible states a `CoreTracker` can have are

* `TrackerCondition.success`: Indicates a successfull completed tracking
* `TrackerCondition.tracking`: The tracking is still in progress
* `TrackerCondition.terminated_maximal_iterations`: Tracking terminated since maximal iterations reached.
* `TrackerCondition.terminated_ill_conditioned`: Tracking terminated since the path was too ill-conditioned.
* `TrackerCondition.terminated_invalid_startvalue`: Tracking terminated since the provided start value was invalid.
* `TrackerCondition.terminated_step_size_too_small`
"""
@enum conditions begin
    success
    tracking
    terminated_maximal_iterations
    terminated_accuracy_limit
    terminated_ill_conditioned
    terminated_invalid_startvalue
    terminated_step_size_too_small
end

end


"""
    is_success(S::TrackerCondition.conditions)

Returns `true` if `S` indicates a success in the path tracking.
"""
is_success(S) = S == TrackerCondition.success

"""
    is_terminated(S::TrackerCondition.conditions)

Returns `true` if `S` indicates that the path tracking got terminated. This is not `true`
if `is_success(S)` is `true`.
"""
is_terminated(S::TrackerCondition.conditions) =
    S ≠ TrackerCondition.tracking && S ≠ TrackerCondition.success

"""
    is_invalid_startvalue(S::TrackerCondition.conditions)

Returns `true` if `S` indicates that the path tracking got terminated since the start
value was not a zero.
"""
is_invalid_startvalue(S::TrackerCondition.conditions) =
    S == TrackerCondition.terminated_invalid_startvalue

"""
    is_tracking(S::TrackerCondition.conditions)

Returns `true` if `S` indicates that the path tracking is not yet finished.
"""
is_tracking(S::TrackerCondition.conditions) = S == TrackerCondition.tracking

mutable struct TrackerState{M<:AbstractMatrix{ComplexF64}}
    x::Vector{ComplexF64} # current x
    x̂::Vector{ComplexF64} # last prediction
    x̄::Vector{ComplexF64} # candidate for new x

    segment::ComplexLineSegment
    s::Double64 # current step length (0 ≤ s ≤ length(segment))
    s′::Double64 # proposed s (0 ≤ s ≤ length(segment))
    Δs_prev::Float64 # previous step size

    accuracy::Float64 # norm(x - x(t))
    ω::Float64 # liptschitz constant estimate, see arxiv:1902.02968
    μ::Float64 # limit accuracy
    τ::Float64 # trust region size
    norm_Δx₀::Float64 # debug info only

    norm::WeightedNorm{InfNorm}
    jacobian::Jacobian{Float64,M}

    condition::TrackerCondition.conditions

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
    x = Vector{ComplexF64}(x₁)
    x̂ = zero(x)
    x̄ = zero(x)

    segment = ComplexLineSegment(t₁, t₀)
    s = s′ = Double64(0.0)

    Δs_prev = 0.0

    accuracy = 0.0
    μ = eps()
    ω = 1.0
    τ = Inf
    norm_Δx₀ = NaN

    JM = Jacobian(zeros(ComplexF64, size(H)))

    condition = TrackerCondition.tracking

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
        norm,
        JM,
        condition,
        accepted_steps,
        rejected_steps,
        last_step_failed,
    )
end

Base.show(io::IO, state::TrackerState) = print_fieldnames(io, state)
Base.show(io::IO, ::MIME"application/prs.juno.inline", state::TrackerState) =
    state
steps(S::TrackerState) = S.accepted_steps + S.rejected_steps

struct TrackerResult{V<:AbstractVector}
    returncode::TrackerCondition.conditions
    x::V
    t::ComplexF64
    accuracy::Float64
    ω::Float64
    μ::Float64
    accepted_steps::Int
    rejected_steps::Int
end

function TrackerResult(state::TrackerState)
    TrackerResult(
        state.condition,
        copy(state.x),
        state.segment[state.s],
        state.accuracy,
        state.ω,
        state.μ,
        state.accepted_steps,
        state.rejected_steps,
    )
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
solution(result::TrackerResult) = result.x

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


struct Tracker{
    H<:AbstractHomotopy,
    P<:AbstractPredictorCache,
    M<:AbstractMatrix{ComplexF64},
}
    homotopy::H
    predictor::P
    corrector::NewtonCorrector
    # these are mutable
    state::TrackerState{M}
    options::TrackerOptions
end

function Tracker(H::AbstractHomotopy; kwargs...)
    x = zeros(ComplexF64, size(H, 2))
    Tracker(H, x, 1.0, 0.0; kwargs...)
end

function Tracker(
    H::AbstractHomotopy,
    x₁::AbstractVector,
    t₁::Number,
    t₀::Number;
    norm::WeightedNorm{InfNorm} = WeightedNorm(InfNorm(), size(H, 2)),
    predictor::AbstractPredictor = Pade21(),
    kwargs...,
)
    options = TrackerOptions(; kwargs...)
    state = TrackerState(H, Vector{ComplexF64}(x₁), t₁, t₀, norm)
    pred_cache = cache(predictor, size(H, 2))
    corrector = NewtonCorrector(options.a, size(H))

    Tracker(H, pred_cache, corrector, state, options)
end

Base.show(io::IO, C::Tracker) = print(io, "Tracker")
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::Tracker) = x
Base.broadcastable(C::Tracker) = Ref(C)

"""
    current_t(tracker::Tracker)

Current `t`.
"""
current_t(tracker::Tracker) = current_t(tracker.state)
current_t(state::TrackerState) = state.segment[state.s]

"""
    current_Δt(tracker::Tracker)

Current step_size `Δt`.
"""
current_Δt(tracker::Tracker) = current_Δt(tracker.state)
current_Δt(state::TrackerState) = state.segment[state.s′] - state.segment[state.s]

#
# function tracker_startsolutions(args...; kwargs...)
#     supported, rest = splitkwargs(
#         kwargs,
#         problem_startsolutions_supported_keywords,
#     )
#     prob, startsolutions = problem_startsolutions(args...; supported...)
#     tracker = Tracker(prob, start_solution_sample(startsolutions); rest...)
#     (tracker = tracker, startsolutions = startsolutions)
# end


# Step Size

_h(a) = 2a * (√(4 * a^2 + 1) - 2a)
## intial step size
function initial_step_size(
    state::TrackerState,
    predictor::AbstractPredictorCache,
    options::TrackerOptions,
)
    a = options.β_a * options.a
    e = state.norm(local_error(predictor))
    Δs₁ = nthroot((√(1 + 2 * _h(a)) - 1) / (state.ω * e), order(predictor))
    Δs₂ = options.β_τ * trust_region(predictor)
    min(Δs₁, Δs₂, length(state.segment))
end

function update_stepsize!(
    state::TrackerState,
    result::NewtonCorrectorResult,
    options::TrackerOptions,
    predictor::AbstractPredictorCache,
)

    a = options.β_a * options.a
    p = order(predictor)
    Δs = state.Δs_prev
    if is_converged(result)
        e = state.norm(local_error(predictor))
        e2 = state.norm(state.x̂, state.x) / Δs^p
        ē = max(e, e2)
        Δs₁ = nthroot((√(1 + 2 * _h(a)) - 1) / (state.ω * ē), p)
        Δs₂ = options.β_τ * trust_region(predictor)
        s′ = min(state.s + min(Δs₁, Δs₂), Double64(length(state.segment)))

        if state.last_step_failed
            state.s′ = min(state.s + Δs, s′)
        else
            state.s′ = s′
        end
    else
        j = result.iters - 2
        Θ_j = nthroot(result.θ, 1 << j)
        state.s′ = state.s + nthroot(
            (√(1 + 2 * _h(0.5a)) - 1) / (√(1 + 2 * _h(Θ_j)) - 1),
            p,
        ) * (state.s′ - state.s)
    end
    nothing
end


function check_terminated!(state::TrackerState, options::TrackerOptions)
    if state.s ≥ length(state.segment)
        state.condition = TrackerCondition.success
        state.s = length(state.segment)
    elseif steps(state) ≥ options.max_steps
        state.condition = TrackerCondition.terminated_maximal_iterations
    elseif state.ω * state.μ > options.a * _h(options.a)
        state.condition = TrackerCondition.terminated_accuracy_limit
    end
    nothing
end


"""
    init!(tracker::NCT, x₁, t₁, t₀)

Setup `tracker` to track `x₁` from `t₁` to `t₀`.
"""
function init!(tracker::Tracker, x₁::AbstractVector, t₁, t₀)
    @unpack state, predictor, corrector, homotopy, options = tracker
    @unpack x, x̄, norm, jacobian = state

    # intialize state
    x .= x₁
    state.segment = ComplexLineSegment(t₁, t₀)
    state.s = state.s′ = state.Δs_prev = 0.0
    state.accuracy = eps()
    state.ω = 1.0
    init!(norm, x)
    init!(jacobian)
    state.condition = TrackerCondition.tracking
    state.accepted_steps = state.rejected_steps = 0
    state.last_step_failed = true

    # compute ω and limit accuracy μ for the start value
    t = state.segment[state.s]
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
    if !isnan(ω)
        state.ω = ω
    end
    if valid
        state.accuracy = μ
        state.μ = max(μ, eps())
    end

    if !valid
        state.condition = TrackerCondition.terminated_invalid_startvalue
    end

    # initialize the predictor
    jacobian!(jacobian.J, homotopy, state.x, t)
    update!(predictor, homotopy, state.x, t, jacobian)
    state.τ = trust_region(predictor)
    # compute initial step size
    state.s′ = initial_step_size(state, predictor, tracker.options)

    tracker
end

function step!(tracker::Tracker, debug::Bool = false)
    @unpack homotopy, corrector, predictor, state, options = tracker
    @unpack x, x̂, x̄, jacobian, norm = state

    t, Δt = current_t(state), current_Δt(state)

    debug && printstyled(
        "\nt: ",
        round(t, sigdigits = 5),
        " Δt: ",
        round(Δt; sigdigits = 5),
        "\n";
        color = :yellow,
    )

    # Use the current approximation of x(t) to obtain estimate
    # x̂ ≈ x(t + Δt) using the provided predictor
    predict!(x̂, predictor, homotopy, x, t, Δt)
    # Correct the predicted value x̂ to obtain x̄.
    # If the correction is successfull we have x̄ ≈ x(t+Δt).
    result = newton!(
        x̄,
        corrector,
        homotopy,
        x̂,
        t + Δt,
        jacobian,
        norm;
        ω = state.ω,
        μ = state.μ,
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

        # If we use a weighted norm, now we update this
        update!(state.norm, x)
        # tell the predictors about the new derivative if they need to update something
        update!(predictor, homotopy, x, t + Δt, jacobian)

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

function track!(tracker::Tracker, x, t₁, t₀; debug::Bool = false)
    init!(tracker, x, t₁, t₀)

    while is_tracking(tracker.state.condition)
        step!(tracker, debug)
    end
end

function track(tracker::Tracker, x, t₁ = 0.0, t₀ = 1.0; debug::Bool = false)
    track!(tracker, x, t₁, t₀; debug = debug)
    TrackerResult(tracker.state)
end
