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

###
### Options and Parameters
###

"""
    TrackerParameters

Parameters that control the performance and robustness characteristics of the path tracking
algorithm. See [^Tim20] for an explanation and derivation of the parameters.
We provide three sets of parameters for common use cases:

* [`DEFAULT_TRACKER_PARAMETERS`](@ref)
* [`FAST_TRACKER_PARAMETERS`](@ref)
* [`CONSERVATIVE_TRACKER_PARAMETERS`](@ref)

[^Tim20]: Timme, S. "Mixed Precision Path Tracking for Polynomial Homotopy Continuation". arXiv:1902.02968 (2020)
"""
Base.@kwdef mutable struct TrackerParameters
    a::Float64 = 0.125
    β_a::Float64 = 1.0
    β_ω::Float64 = 40.0
    β_τ::Float64 = 0.6
    strict_β_τ::Float64 = min(β_τ, 0.45)
    min_newton_iters::Int = 2
end
Base.show(io::IO, TP::TrackerParameters) = print_fieldnames(io, TP)

"The default [`TrackerParameters`](@ref) which have a good balance between robustness and efficiency."
const DEFAULT_TRACKER_PARAMETERS = TrackerParameters()
"[`TrackerParameters`](@ref) which trade speed against a higher chance of path jumping."
const FAST_TRACKER_PARAMETERS = TrackerParameters(β_τ = 0.8, β_ω = 10.0)
"[`TrackerParameters`](@ref) which trade robustness against some speed."
const CONSERVATIVE_TRACKER_PARAMETERS =
    TrackerParameters(β_ω = 100.0, β_τ = 0.45, strict_β_τ = 0.35)


"""
    TrackerOptions(; options...)

The set of options for a [`Tracker`](@ref).

## Options

* `automatic_differentiation::Int = 3`: The value `automatic_differentiation` determines
  up to which order the derivative is computed using automatic differentiation.
  Otherwise numerical differentiation is used. The automatic differentiation results
  in additional compilation time, however for numerically challenging paths it is strongly
  recommended to use `automatic_differentiation = 3`.
* `max_steps::Int = 1_000`: The maximal number of steps a tracker attempts
* `max_step_size::Float64 = Inf`: The maximal size of a step
* `max_initial_step_size::Float64 = Inf`: The maximal size of the first step
* `min_step_size::Float64 = 1e-48`: The minimal step size. If a smaller step size would
  be necessary, then the tracking gets terminated.
* `extended_precision::Bool = true`: Whether to allow for the use of extended precision,
  if necessary, in some computations. This can greatly improve the ability to track
  numerically difficult paths.
* `terminate_cond::Float64 = 1e13`: If the relative component-wise condition number
  `cond(H_x, ẋ)` is larger than `terminate_cond` then the path is terminated as too
  ill-conditioned.
* `parameters::Union{Symbol,TrackerParameters} = :default` Set the
  [`TrackerParameters`](@ref) to control the performance of the path tracking algorithm.
  The values `:default`, `:conservative` and `:fast` are shorthands for using
  [`DEFAULT_TRACKER_PARAMETERS`](@ref), [`CONSERVATIVE_TRACKER_PARAMETERS`](@ref) resp.
  [`FAST_TRACKER_PARAMETERS`](@ref).
"""
mutable struct TrackerOptions
    automatic_differentiation::Int
    max_steps::Int
    max_step_size::Float64
    max_initial_step_size::Float64
    extended_precision::Bool
    min_step_size::Float64
    terminate_cond::Float64
    parameters::TrackerParameters
end
function TrackerOptions(;
    automatic_differentiation::Int = 3,
    max_steps::Int = 1_000,
    max_step_size::Float64 = Inf,
    max_initial_step_size::Float64 = Inf,
    extended_precision::Bool = true,
    min_step_size::Float64 = 1e-48,
    terminate_cond::Float64 = 1e13,
    parameters::Union{Symbol,TrackerParameters} = :default,
)

    if parameters isa Symbol
        if parameters == :default
            parameters = DEFAULT_TRACKER_PARAMETERS
        elseif parameters == :fast
            parameters = FAST_TRACKER_PARAMETERS
        elseif parameters == :conservative
            parameters = CONSERVATIVE_TRACKER_PARAMETERS
        else
            throw(ArgumentError("Unsupported `parameters` value $parameters"))
        end
    end

    TrackerOptions(
        automatic_differentiation,
        max_steps,
        max_step_size,
        max_initial_step_size,
        extended_precision,
        min_step_size,
        terminate_cond,
        parameters,
    )
end

Base.show(io::IO, opts::TrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::TrackerOptions) = opts


###
### TrackerCode
###

@doc """
    TrackerCode

The possible states a `CoreTracker` can have are of type `TrackerCode.codes` and can be

* `TrackerCode.success`: Indicates a successfull tracking.
* `TrackerCode.tracking`: The tracking is still in progress.
* `TrackerCode.terminated_accuracy_limit`: Tracking terminaed since the accuracy was insufficient.
* `TrackerCode.terminated_invalid_startvalue`: Tracking terminated since the provided start value was invalid.
* `TrackerCode.terminated_ill_conditioned`: Tracking terminated since the path was too ill-conditioned.
* `TrackerCode.terminated_max_steps`: Tracking terminated since maximal number of steps is reached.
* `TrackerCode.terminated_step_size_too_small`: Trackint terminated since the step size was too small.
* `TrackerCode.terminated_unknown`: An unintended error occured. Please consider reporting an issue.
"""
module TrackerCode

@enum codes begin
    tracking
    success
    terminated_max_steps
    terminated_accuracy_limit
    terminated_ill_conditioned
    terminated_invalid_startvalue
    terminated_step_size_too_small
    terminated_unknown
end

end

"""
    is_success(code::TrackerCode.codes)

Returns `true` if `code` indicates a success in the path tracking.
"""
is_success(S::TrackerCode.codes) = S == TrackerCode.success

"""
    is_terminated(code::TrackerCode.codes)

Returns `true` if `code` indicates that the path tracking got terminated.
"""
is_terminated(S::TrackerCode.codes) = S ≠ TrackerCode.tracking && S ≠ TrackerCode.success

"""
    is_invalid_startvalue(code::TrackerCode.codes)

Returns `true` if `code` indicates that the path tracking got terminated since the start
value was not a zero.
"""
is_invalid_startvalue(S::TrackerCode.codes) = S == TrackerCode.terminated_invalid_startvalue

"""
    is_tracking(code::TrackerCode.codes)

Returns `true` if `code` indicates that the path tracking is not yet finished.
"""
is_tracking(S::TrackerCode.codes) = S == TrackerCode.tracking

###
### RESULT
###

"""
    TrackerResult

Containing the result of tracking a path with a [`Tracker`](@ref).

## Fields

* `return_code::Symbol`: A code indicating whether the tracking was successfull (`:success`).
  See [`TrackerCode`](@ref) for all possible values.
* `solution::V`: The solution when the tracking stopped.
* `t::ComplexF64`: The value of `t` when the tracking stopped.
* `accuracy::Float64`: Estimate of the relative accuracy of the `solution`.
* `accepted_steps::Int`: Number of steps that got accepted.
* `rejected_steps::Int`: Number of steps that got rejected.
* `extended_precision::Bool`: Indicate whether extended precision is necessary to achieve
  the accuracy of the `solution`.
* `extended_precision_used::Bool`: This is `true` if during the tracking at any point
  extended precision was used.
"""
struct TrackerResult
    return_code::Symbol
    solution::Vector{ComplexF64}
    t::ComplexF64
    accuracy::Float64
    accepted_steps::Int
    rejected_steps::Int
    extended_precision::Bool
    extended_precision_used::Bool
    ω::Float64
    μ::Float64
    τ::Float64
end

Base.show(io::IO, result::TrackerResult) = print_fieldnames(io, result)
Base.show(io::IO, ::MIME"application/prs.juno.inline", result::TrackerResult) = result

"""
    is_success(result::TrackerResult)

Returns `true` if the path tracking was successfull.
"""
is_success(R::TrackerResult) = R.return_code == :success

"""
    is_invalid_startvalue(result::TrackerResult)

Returns `true` if the path tracking failed since the start value was invalid.
"""
is_invalid_startvalue(R::TrackerResult) = R.return_code == :invalid_startvalue

"""
    solution(result::TrackerResult)

Returns the solutions obtained by the `Tracker`.
"""
solution(result::TrackerResult) = result.solution

"""
    steps(result::TrackerResult)

Returns the number of steps done.
"""
steps(result::TrackerResult) = accepted_steps(result) + rejected_steps(result)

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

###
### STATE
###

mutable struct TrackerState{M<:AbstractMatrix{ComplexF64}}
    x::Vector{ComplexF64} # current x
    x̂::Vector{ComplexF64} # last prediction
    x̄::Vector{ComplexF64} # candidate for new x
    # internal step size
    segment_stepper::SegmentStepper
    Δs_prev::Float64 # previous step size
    # path tracking algorithm
    accuracy::Float64 # norm(x - x(t))
    ω::Float64 # liptschitz constant estimate, see arxiv:1902.02968
    μ::Float64 # limit accuracy
    τ::Float64 # trust region size
    norm_Δx₀::Float64 # debug info only
    extended_prec::Bool
    used_extended_prec::Bool
    norm::WeightedNorm{InfNorm}
    use_strict_β_τ::Bool

    jacobian::Jacobian{M}
    cond_J_ẋ::Float64 # estimate of cond(H(x(t),t), ẋ(t))
    code::TrackerCode.codes

    # buffer for the computation of derivatives
    tx⁰::TaylorVector{1,ComplexF64}
    tx¹::TaylorVector{2,ComplexF64}
    tx²::TaylorVector{3,ComplexF64}
    tx³::TaylorVector{4,ComplexF64}
    tx⁴::TaylorVector{5,ComplexF64}
    trust_tx::Vector{Bool}
    u::Vector{ComplexF64}

    # statistics
    accepted_steps::Int
    rejected_steps::Int
    last_steps_failed::Int
end

function TrackerState(H, x₁::AbstractVector, norm::WeightedNorm{InfNorm})
    x = convert(Vector{ComplexF64}, x₁)
    x̂ = zero(x)
    x̄ = zero(x)

    m, n = size(H)
    tx⁴ = TaylorVector{5}(ComplexF64, n)
    tx⁰ = TaylorVector{1}(tx⁴)
    tx¹ = TaylorVector{2}(tx⁴)
    tx² = TaylorVector{3}(tx⁴)
    tx³ = TaylorVector{4}(tx⁴)
    trust_tx = [true, true, true, true]
    u = zeros(ComplexF64, m)
    segment_stepper = SegmentStepper(1.0, 0.0)
    Δs_prev = 0.0
    accuracy = 0.0
    μ = eps()
    ω = 1.0
    τ = Inf
    norm_Δx₀ = NaN
    use_strict_β_τ = false
    used_extended_prec = extended_prec = false
    jacobian = Jacobian(zeros(ComplexF64, size(H)))
    cond_J_ẋ = NaN
    code = TrackerCode.tracking
    u = zeros(ComplexF64, size(H, 1))
    accepted_steps = rejected_steps = 0
    last_steps_failed = 0

    TrackerState(
        x,
        x̂,
        x̄,
        segment_stepper,
        Δs_prev,
        accuracy,
        ω,
        μ,
        τ,
        norm_Δx₀,
        extended_prec,
        used_extended_prec,
        norm,
        use_strict_β_τ,
        jacobian,
        cond_J_ẋ,
        code,
        tx⁰,
        tx¹,
        tx²,
        tx³,
        tx⁴,
        trust_tx,
        u,
        accepted_steps,
        rejected_steps,
        last_steps_failed,
    )
end

Base.show(io::IO, state::TrackerState) = print_fieldnames(io, state)
Base.show(io::IO, ::MIME"application/prs.juno.inline", state::TrackerState) = state
function Base.getproperty(state::TrackerState, sym::Symbol)
    if sym === :t
        return getfield(state, :segment_stepper).t
    elseif sym == :Δt
        return getfield(state, :segment_stepper).Δt
    elseif sym == :t′
        return getfield(state, :segment_stepper).t′
    elseif sym == :last_step_failed
        return getfield(state, :last_steps_failed) > 0
    else # fallback to getfield
        return getfield(state, sym)
    end
end

steps(S::TrackerState) = S.accepted_steps + S.rejected_steps

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
struct Tracker{H<:AbstractHomotopy,N,M<:AbstractMatrix{ComplexF64}}
    homotopy::H
    predictor::Pade21
    corrector::NewtonCorrector
    # these are mutable
    state::TrackerState{M}
    options::TrackerOptions
    AD::Val{N}
    ND::NumericalDifferentiation
end

Tracker(H::ModelKit.Homotopy; kwargs...) = Tracker(ModelKitHomotopy(H); kwargs...)
function Tracker(
    H::AbstractHomotopy,
    x::AbstractVector = zeros(size(H, 2));
    weighted_norm_options::WeightedNormOptions = WeightedNormOptions(),
    options::TrackerOptions = TrackerOptions(),
    AD::Int = options.automatic_differentiation,
)
    norm = WeightedNorm(ones(size(H, 2)), InfNorm(), weighted_norm_options)
    state = TrackerState(H, x, norm)
    predictor = Pade21(size(H, 2))
    corrector = NewtonCorrector(options.parameters.a, state.x, size(H, 1))

    if !(0 ≤ AD ≤ 4)
        throw(ArgumentError("`automatic_differentiation` has to be between 0 and 4."))
    end
    ND = NumericalDifferentiation(state.x, size(H, 1))
    Tracker(H, predictor, corrector, state, options, Val(AD), ND)
end

function Base.show(io::IO, C::Tracker{<:Any,N}) where {N}
    print(io, "Tracker with automatic_differentiation=$N")
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::Tracker) = x
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
status(tracker::Tracker) = tracker.state.code

LA.cond(tracker::Tracker) = tracker.state.cond_J_ẋ

function LA.cond(tracker::Tracker, x, t, d_l = nothing, d_r = nothing)
    J = tracker.state.jacobian
    evaluate_and_jacobian!(tracker.corrector.r, matrix(J), tracker.homotopy, x, t)
    updated!(J)
    LA.cond(J, d_l, d_r)
end
# Step Size

_h(a) = 2a * (√(4 * a^2 + 1) - 2a)

"Like `min(a,b)`` but ignoring any `NaN` values."
nanmin(a, b) = isnan(a) ? b : (isnan(b) ? a : min(a, b))

# intial step size
function initial_step_size(state::TrackerState, predictor::Pade21, options::TrackerOptions)
    a = options.parameters.β_a * options.parameters.a
    ω = options.parameters.β_ω * state.ω
    e = state.norm(local_error(predictor))
    τ = trust_region(predictor)
    Δs₁ = nthroot((√(1 + 2 * _h(a)) - 1) / (ω * e), order(predictor))
    Δs₂ = options.parameters.β_τ * τ
    Δs = nanmin(Δs₁, Δs₂)
    min(Δs, options.max_step_size, options.max_initial_step_size)
end

function update_stepsize!(
    state::TrackerState,
    result::NewtonCorrectorResult,
    options::TrackerOptions,
    predictor::Pade21;
    ad_for_error_estimate::Bool = true,
)
    a = options.parameters.β_a * options.parameters.a
    ω = options.parameters.β_ω * state.ω
    p = order(predictor)
    τ = state.τ #trust_region(predictor)
    if is_converged(result)
        # e₁ = state.norm(local_error(predictor))
        # If we don't use automatic_differentiation for the 4th derivative
        # this can be completely of. So check against the actual error of the
        # last step
        if ad_for_error_estimate
            e = state.norm(local_error(predictor))
            Δs₁ = nthroot((√(1 + 2 * _h(a)) - 1) / (ω * e), p)
        else
            e₁ = nthroot(state.norm(local_error(predictor)), p)
            e₂ = nthroot(state.norm(state.x, state.x̂), p) / state.Δs_prev
            e = nanmin(e₁, e₂)
            Δs₁ = nthroot((√(1 + 2 * _h(a)) - 1) / ω, p) / e
        end
        Δs₂ = options.parameters.β_τ * τ
        if state.use_strict_β_τ || dist_to_target(state.segment_stepper) < Δs₂
            Δs₂ = options.parameters.strict_β_τ * τ
        end
        Δs = min(nanmin(Δs₁, Δs₂), options.max_step_size)
        if state.last_step_failed
            Δs = min(Δs, state.Δs_prev)
        end
    else
        j = result.iters - 2
        Θ_j = nthroot(result.θ, 1 << j)
        if isnan(Θ_j) || result.return_code == NEWT_SINGULARITY || isnan(result.accuracy)
            Δs = 0.25 * state.segment_stepper.Δs
        else
            Δs =
                nthroot((√(1 + 2 * _h(0.5a)) - 1) / (√(1 + 2 * _h(Θ_j)) - 1), p) *
                state.segment_stepper.Δs
        end
    end
    propose_step!(state.segment_stepper, Δs)
    nothing
end


function check_terminated!(state::TrackerState, options::TrackerOptions)
    if state.extended_prec || !options.extended_precision
        @unpack a, min_newton_iters = options.parameters
        tol_acc = a^(2^min_newton_iters - 1) * _h(a)
    else
        tol_acc = Inf
    end

    if is_done(state.segment_stepper)
        state.code = TrackerCode.success
    elseif steps(state) ≥ options.max_steps
        state.code = TrackerCode.terminated_max_steps
    elseif state.ω * state.μ > tol_acc
        state.code = TrackerCode.terminated_accuracy_limit
    elseif state.last_steps_failed ≥ 3 && state.cond_J_ẋ > options.terminate_cond
        state.code = TrackerCode.terminated_ill_conditioned
    elseif state.segment_stepper.Δs < options.min_step_size
        state.code = TrackerCode.terminated_step_size_too_small
    end
    nothing
end

function compute_derivatives_and_update_predictor!(tracker::Tracker)
    compute_derivatives!(tracker)
    update!(tracker.predictor, tracker.state.tx⁴, tracker.state.trust_tx)
end

function compute_derivatives!(
    tracker::Tracker;
    min_acc::Float64 = sqrt(eps()),
    max_iters::Int = 5,
)
    @unpack homotopy, state, options, AD, ND = tracker
    # unpack stuff to make the rest easier to read
    @unpack u, x, t, tx⁰, tx¹, tx², tx³, tx⁴, jacobian, norm, τ = state
    x⁰, x¹, x², x³, x⁴ = vectors(tx⁴)
    τ̂ = 1.0 #0.5τ

    dist_target = dist_to_target(state.segment_stepper)

    x⁰ .= x
    # compute all taylor coefficients x¹, x², x³, x⁴
    state.trust_tx[1] = taylor!(
        u,
        Val(1),
        homotopy,
        tx⁰,
        t,
        AD,
        ND;
        cond = state.cond_J_ẋ,
        dist_to_target = dist_target,
    )
    u .= .-u
    LA.ldiv!(x¹, jacobian, u)

    # Check if we have to do iterative refinment for all the others as well
    δ = fixed_precision_iterative_refinement!(x¹, workspace(jacobian), u, norm)
    state.cond_J_ẋ = δ / eps()
    iterative_refinement = δ > min_acc
    if iterative_refinement
        (acc, diverged) = iterative_refinement!(
            x¹,
            jacobian,
            u,
            norm;
            tol = min_acc,
            max_iters = max_iters,
        )
        if diverged
            state.code = TrackerCode.terminated_ill_conditioned
        end
    end

    state.trust_tx[2] = taylor!(
        u,
        Val(2),
        homotopy,
        tx¹,
        t,
        AD,
        ND;
        cond = state.cond_J_ẋ,
        dist_to_target = dist_target,
        incremental = true,
    )
    u .= .-u
    LA.ldiv!(x², jacobian, u)
    iterative_refinement &&
    iterative_refinement!(x², jacobian, u, norm; tol = min_acc, max_iters = max_iters)

    state.trust_tx[3] = taylor!(
        u,
        Val(3),
        homotopy,
        tx²,
        t,
        AD,
        ND;
        cond = state.cond_J_ẋ,
        dist_to_target = dist_target,
        incremental = true,
    )
    u .= .-u
    LA.ldiv!(x³, jacobian, u)
    iterative_refinement &&
    iterative_refinement!(x³, jacobian, u, norm; tol = min_acc, max_iters = max_iters)

    state.trust_tx[4] = taylor!(
        u,
        Val(4),
        homotopy,
        tx³,
        t,
        AD,
        ND;
        dist_to_target = dist_target,
        cond = state.cond_J_ẋ,
        incremental = true,
    )
    u .= .-u
    LA.ldiv!(x⁴, jacobian, u)

    state
end

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
    ω::Float64 = NaN,
    μ::Float64 = NaN,
    τ::Float64 = Inf,
    max_initial_step_size::Float64 = Inf,
    keep_steps::Bool = false,
)
    @unpack state, predictor, corrector, homotopy, options = tracker
    @unpack x, x̄, norm, jacobian = state

    # intialize state
    set_solution!(x, homotopy, x₁, t₁)
    init!(state.segment_stepper, t₁, t₀)
    state.Δs_prev = 0.0
    state.accuracy = eps()
    state.ω = 1.0
    state.used_extended_prec = state.extended_prec = false
    state.use_strict_β_τ = false
    init!(norm, x)
    init!(jacobian)
    state.code = TrackerCode.tracking
    if !keep_steps
        state.accepted_steps = state.rejected_steps = 0
    end
    state.last_steps_failed = 0

    # compute ω and limit accuracy μ for the start value
    t = state.t
    if isnan(ω) || isnan(μ)
        a = options.parameters.a
        valid, ω, μ = init_newton!(x̄, corrector, homotopy, x, t, jacobian, norm; a = a)
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
        state.code = TrackerCode.terminated_invalid_startvalue
        return false
    end

    # initialize the predictor
    state.τ = τ
    evaluate_and_jacobian!(corrector.r, workspace(jacobian), homotopy, state.x, t)
    updated!(jacobian)
    compute_derivatives_and_update_predictor!(tracker)
    state.τ = trust_region(predictor)
    # compute initial step size
    Δs = initial_step_size(state, predictor, tracker.options)
    Δs = min(Δs, max_initial_step_size)
    propose_step!(state.segment_stepper, Δs)

    is_tracking(state.code)
end

function init!(tracker::Tracker, t₀::Number)
    @unpack state, predictor, options = tracker
    state.code = TrackerCode.tracking
    init!(state.segment_stepper, state.t, t₀)
    Δs = initial_step_size(state, predictor, tracker.options)
    propose_step!(state.segment_stepper, Δs)
    state.Δs_prev = 0.0

    tracker
end

function update_precision!(tracker::Tracker, μ_low)
    @unpack homotopy, corrector, state, options = tracker
    @unpack μ, ω, x, t, jacobian, norm = state
    @unpack a = options.parameters

    options.extended_precision || return false

    if state.extended_prec && !isnan(μ_low)
        # check if we can go low again
        if μ_low * ω < a^7 * _h(a)
            state.extended_prec = false
            state.μ = μ_low
        end
    elseif μ * ω > a^5 * _h(a)
        state.extended_prec = true
        state.used_extended_prec = true
        # do two refinement steps
        for i = 1:2
            μ = extended_prec_refinement_step!(x, corrector, homotopy, x, t, jacobian, norm)
        end
        state.μ = max(μ, eps())
    end

    state.extended_prec
end

function refine_current_solution!(tracker)
    @unpack homotopy, corrector, state, options = tracker
    @unpack x, t, jacobian, norm = state

    extended_prec_refinement_step!(x, corrector, homotopy, x, t, jacobian, norm)
    extended_prec_refinement_step!(x, corrector, homotopy, x, t, jacobian, norm)
end

"""
    step!(tracker::Tracker, debug::Bool = false)

Perform a single tracking step. Returns `true` if the step was accepted.
"""
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
    predict!(x̂, predictor, homotopy, state.tx³, Δt)
    # if t′ == 0
    #     e = state.norm(local_error(predictor)) * state.segment_stepper.Δs^4
    # end
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
        extended_precision = state.extended_prec,
    )
    if debug
        printstyled(result, "\n"; color = is_converged(result) ? :green : :red)
    end
    if is_converged(result)
        # move forward
        x .= x̄
        state.Δs_prev = state.segment_stepper.Δs
        step_success!(state.segment_stepper)
        state.accuracy = result.accuracy
        state.μ = max(result.accuracy, eps())

        no_ω_update = result.ω == state.ω

        state.ω = max(result.ω, 0.5 * state.ω)
        update_precision!(tracker, result.μ_low)
        compute_derivatives_and_update_predictor!(tracker)

        # # Update other state
        # if no_ω_update && options.automatic_differentiation < 3
        #     state.τ = max(trust_region(predictor), 5 * state.τ)
        # else
        state.τ = trust_region(predictor)
        # end

        state.accepted_steps += 1
        state.last_steps_failed = 0
    else
        # Step failed, so we have to try with a new (smaller) step size
        state.rejected_steps += 1
        state.last_steps_failed += 1
    end
    state.norm_Δx₀ = result.norm_Δx₀
    update_stepsize!(
        state,
        result,
        options,
        predictor;
        ad_for_error_estimate = tracker.AD isa Val{4},
    )

    check_terminated!(state, options)
    !state.last_step_failed
end

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
    ω::Float64 = NaN,
    μ::Float64 = NaN,
    τ::Float64 = Inf,
    keep_steps::Bool = false,
    max_initial_step_size::Float64 = Inf,
    debug::Bool = false,
)
    init!(
        tracker,
        x,
        t₁,
        t₀;
        ω = ω,
        μ = μ,
        τ = τ,
        keep_steps = keep_steps,
        max_initial_step_size = max_initial_step_size,
    )

    while is_tracking(tracker.state.code)
        step!(tracker, debug)
    end

    tracker.state.code
end

function track!(tracker::Tracker, r::TrackerResult, t₁ = 1.0, t₀ = 0.0; debug::Bool = false)
    track!(tracker, solution(r), t₁, t₀; debug = debug, ω = r.ω, μ = r.μ, τ = r.τ)
end
function track!(tracker::Tracker, t₀; debug::Bool = false)
    init!(tracker, t₀)

    while is_tracking(tracker.state.code)
        step!(tracker, debug)
    end

    tracker.state.code
end

function TrackerResult(H::AbstractHomotopy, state::TrackerState)
    TrackerResult(
        Symbol(state.code),
        get_solution(H, state.x, state.t),
        state.t,
        state.accuracy,
        state.accepted_steps,
        state.rejected_steps,
        state.extended_prec,
        state.used_extended_prec,
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
