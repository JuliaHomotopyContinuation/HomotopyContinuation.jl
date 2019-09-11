export CoreTracker,
       CoreTrackerResult,
       CoreTrackerStatus,
       is_success,
       is_terminated,
       is_tracking,
       CoreTrackerOptions,
       coretracker,
       coretracker_startsolutions,
       affine_tracking,
       track,
       track!,
       setup!, # deprecated
       init!,
       iterator,
       current_x,
       current_t,
       current_Δt,
       steps,
       iters,
       status,
       accuracy,
       max_corrector_iters,
       max_step_size,
       refinement_accuracy,
       max_refinement_iters,
       set_accuracy!,
       set_max_corrector_iters!,
       set_refinement_accuracy!,
       set_max_refinement_iters!,
       set_max_step_size!,
       digits_lost,
       options


const coretracker_supported_keywords = [
    :corrector,
    :predictor,
    :patch,
    :initial_step_size,
    :min_step_size,
    :max_step_size,
    :accuracy,
    :refinement_accuracy,
    :max_corrector_iters,
    :max_refinement_iters,
    :max_steps,
    :simple_step_size_alg,
    :auto_scaling,
    :terminate_ill_conditioned,
    :log_transform,
    :precision,
    :steps_jacobian_info_update,
    :max_lost_digits
]


####################
# CoreTrackerState #
####################

@doc """
    enum CoreTrackerStatus

The possible states the `CoreTracker` can achieve are

* `CT_SUCCESS`: Indicates a successfull completed tracking
* `CT_TRACKING`: The tracking is still in progress
* `CT_TERMINATED_MAX_ITERS`: Tracking terminated since maximal iterations reached.
* `CT_TERMINATED_INVALID_STARTVALUE`: Tracking terminated since the provided start value was invalid.
* `CT_TERMINATED_STEP_SIZE_TOO_SMALL`
* `CT_TERMINATED_SINGULARITY`
* `CT_TERMINATED_ILL_CONDITIONED`
""" @enum CoreTrackerStatus begin
    CT_SUCCESS
    CT_TRACKING
    CT_TERMINATED_MAX_ITERS
    CT_TERMINATED_INVALID_STARTVALUE
    CT_TERMINATED_STEP_SIZE_TOO_SMALL
    CT_TERMINATED_SINGULARITY
    CT_TERMINATED_ILL_CONDITIONED
end

"""
    is_success(S::CoreTrackerStatus)

Returns `true` if `S` indicates a success in the path tracking.
"""
is_success(S::CoreTrackerStatus) = S == CT_SUCCESS

"""
    is_terminated(S::CoreTrackerStatus)

Returns `true` if `S` indicates that the path tracking got terminated. This is not `true`
if `is_success(S)` is `true`.
"""
is_terminated(S::CoreTrackerStatus) = S ≠ CT_TRACKING && S ≠ CT_SUCCESS

"""
    is_tracking(S::CoreTrackerStatus)

Returns `true` if `S` indicates that the path tracking is not yet finished.
"""
is_tracking(S::CoreTrackerStatus) = S == CT_TRACKING

####################
# CoreTrackerResult #
####################

"""
     CoreTrackerResult{V<:AbstractVector}

Containing the result of a tracked path. The fields are
* `returncode::CTStatus` If the tracking was successfull then it is `CT_SUCCESS`.
* `x::V` The result.
* `t::ComplexF64` The `t` when the path tracker stopped.
* `accuracy::Float64`: The estimated accuracy of `x`.
"""
struct CoreTrackerResult{V<:AbstractVector}
    returncode::CoreTrackerStatus
    x::V
    t::ComplexF64
    accuracy::Float64
    accepted_steps::Int
    rejected_steps::Int
end

function CoreTrackerResult(tracker)
    state = tracker.state
    CoreTrackerResult(
        state.status,
        copy(state.x),
        state.segment[state.s],
        state.accuracy,
        state.accepted_steps,
        state.rejected_steps
    )
end

"""
    is_success(R::CoreTrackerResult)

Returns `true` if the path tracking was successfull.
"""
is_success(R::CoreTrackerResult) = is_success(R.returncode)

Base.show(io::IO, result::CoreTrackerResult) = print_fieldnames(io, result)
Base.show(io::IO, ::MIME"application/prs.juno.inline", result::CoreTrackerResult) = result


######################
# CoreTrackerOptions #
######################


"""
    CoreTrackerOptions

The set of options set for a [`CoreTracker`](@ref). See the description of [`CoreTracker`](@ref)
for all possible options.
"""
mutable struct CoreTrackerOptions
    accuracy::Float64
    max_corrector_iters::Int
    refinement_accuracy::Float64
    max_refinement_iters::Int
    max_steps::Int
    initial_step_size::Float64
    min_step_size::Float64
    max_step_size::Float64
    simple_step_size_alg::Bool
    update_patch::Bool
    max_lost_digits::Float64
    terminate_ill_conditioned::Bool
    steps_jacobian_info_update::Int
    # Not changeable options
    precision::PrecisionOption
    logarithmic_time_scale::Bool
end

function CoreTrackerOptions(
    ;
    accuracy = 1e-7,
    initial_step_size = 0.1,
    max_corrector_iters::Int = 2,
    max_refinement_iters = 5,
    max_step_size = Inf,
    min_step_size = 1e-14,
    parameter_homotopy = false,
    max_steps = parameter_homotopy ? 10_000 : 1_000,
    precision::PrecisionOption = PRECISION_FIXED_64,
    max_lost_digits = default_max_lost_digits(precision, accuracy),
    refinement_accuracy = 1e-8,
    simple_step_size_alg = false,
    steps_jacobian_info_update::Int = 3,
    terminate_ill_conditioned::Bool = true,
    update_patch = true,
    logarithmic_time_scale = false
)

    CoreTrackerOptions(
        accuracy,
        max_corrector_iters,
        refinement_accuracy,
        max_refinement_iters,
        max_steps,
        initial_step_size,
        min_step_size,
        max_step_size,
        simple_step_size_alg,
        update_patch,
        max_lost_digits,
        terminate_ill_conditioned,
        steps_jacobian_info_update,
        precision,
        logarithmic_time_scale
    )
end

Base.show(io::IO, opts::CoreTrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::CoreTrackerOptions) = opts

const eps_double_64 = Float64(eps(Double64))
function default_max_lost_digits(prec::PrecisionOption, accuracy::Float64)
    if prec == PRECISION_FIXED_64
        # maximal_digits_available - digits necessary - buffer
        -log10(eps()) + log10(accuracy) + 1
    else
        # if higher precision is available we will more like be constrained
        # by the fact that the jacobian cannot be too ill-conditioned
        min(-log10(eps()) - 3, -log10(eps_double_64) + log10(accuracy) - 3)
    end
end

####################
# CoreTrackerState #
####################

mutable struct CoreTrackerState{
    T,
    AV<:AbstractVector{Complex{T}},
    AN<:AbstractNorm,
    MaybePatchState<:Union{AbstractAffinePatchState,Nothing}
}
    x::AV # current x
    x̂::AV # last prediction
    x̄::AV # candidate for new x
    ẋ::Vector{Complex{T}} # derivative at current x

    segment::ComplexSegment
    s::Float64 # current step length (0 ≤ s ≤ length(segment))
    Δs::Float64 # current step size
    Δs_prev::Float64 # previous step size

    accuracy::Float64 # norm(x - x(t))
    ω::Float64 # liptschitz constant estimate, see arxiv:1902.02968
    limit_accuracy::Float64

    norm::AN
    jacobian::JacobianMonitor{T}
    residual::Vector{Float64}

    steps_jacobian_info_update::Int
    status::CoreTrackerStatus
    patch::MaybePatchState

    accepted_steps::Int
    rejected_steps::Int
    last_step_failed::Bool
    consecutive_successfull_steps::Int
end

function CoreTrackerState(
    H,
    x₁::AbstractVector,
    t₁,
    t₀,
    options::CoreTrackerOptions,
    patch::Union{Nothing,AbstractAffinePatchState},
    norm::AbstractNorm
)
    x = isa(x₁, SVector) ? Vector(x₁) : copy(x₁)
    x̂, x̄ = copy(x), copy(x)
    ẋ = Vector(x)
    segment = ComplexSegment(promote(complex(t₁), complex(t₀))...)
    s = 0.0
    Δs = convert(
        Float64,
        min(options.initial_step_size, length(segment), options.max_step_size)
    )
    Δs_prev = 0.0
    accuracy = limit_accuracy = 0.0
    ω = 1.0
    JM = JacobianMonitor(jacobian(H, x, t₁))
    residual = zeros(first(size(H)))
    steps_jacobian_info_update = 0
    accepted_steps = rejected_steps = 0
    status = CT_TRACKING
    last_step_failed = false
    consecutive_successfull_steps = 0
    CoreTrackerState(
        x,
        x̂,
        x̄,
        ẋ,
        segment,
        s,
        Δs,
        Δs_prev,
        accuracy,
        ω,
        limit_accuracy,
        norm,
        JM,
        residual,
        steps_jacobian_info_update,
        status,
        patch,
        accepted_steps,
        rejected_steps,
        last_step_failed,
        consecutive_successfull_steps
    )
end

Base.show(io::IO, state::CoreTrackerState) = print_fieldnames(io, state)
Base.show(io::IO, ::MIME"application/prs.juno.inline", state::CoreTrackerState) = state


##################
## CoreTracker ##
##################
"""
    CoreTracker(problem::AbstractProblem, x₁; kwargs...)

Construct a `CoreTracker` from the given `problem` to track elements of type `x₁`.
The path is tracked using a predictor-corrector scheme. The recommended methods to construct
a `CoreTracker` are [`CoreTracker`](@ref) and [`coretracker_startsolutions`](@ref).
Note that a `CoreTracker` is also a (mutable) iterator.

## Keyword arguments
* `accuracy=1e-7`: The accuracy required during the path tracking.
* `corrector::AbstractCorrector`: The corrector used during in the predictor-corrector scheme. The default is [`NewtonCorrector`](@ref).
* `initial_step_size=0.1`: The size of the first step.
* `max_corrector_iters=2`: The maximal number of correction steps in a single step.
* `max_lost_digits::Real`: The tracking is terminated if we estimate that we loose more than `max_lost_digits` in the linear algebra steps. This threshold depends on the `precision` argument.
* `max_refinement_iters=5`: The maximal number of correction steps used to refine the final value.
* `max_steps=1_000`: The maximal number of iterations the path tracker has available. Note that this changes to `10_000` for parameter homotopies.
* `max_step_size=Inf`: The maximal step size.
* `min_step_size=1e-14`: The minimal step size.
* `precision::PrecisionOption=PRECISION_FIXED_64`: The precision used for evaluating the residual in Newton's method.
* `predictor::AbstractPredictor`: The predictor used during in the predictor-corrector scheme. The default is [`Heun`](@ref)()`.
* `refinement_accuracy=1e-8`: The precision required for the final value.
* `simple_step_size_alg=false`: Use a more simple step size algorithm.
* `steps_jacobian_info_update::Int=1`: Every n-th step a linear system will be solved using a QR factorization to obtain an estimate for the condition number of the Jacobian.
* `terminate_ill_conditioned::Bool=true`: Indicates whether the path tracking should be terminated for ill-conditioned paths. A path is considerd ill-conditioned if the condition number of the Jacobian is larger than ≈1e14 or if it is larger than 1e`max_lost_digits`.
"""
struct CoreTracker{
    T,
    AV<:AbstractVector{Complex{T}},
    AN<:AbstractNorm,
    P<:AbstractPredictorCache,
    AP<:Union{Nothing,AbstractAffinePatchState},
    H<:HomotopyWithCache
}
    homotopy::H
    predictor::P
    corrector::NewtonCorrectorCache{T}
    # these are mutable
    state::CoreTrackerState{T,AV,AN,AP}
    options::CoreTrackerOptions
    # buffer
    dt::Vector{Complex{T}}
end

function CoreTracker(
    prob::AbstractProblem,
    x₁;
    log_homotopy::Bool = false,
    logarithmic_time_scale = nothing, kwargs...
)
    if log_homotopy
        H = LogHomotopy(prob.homotopy)
        CoreTracker(
            H,
            embed(prob, x₁),
            complex(0.0),
            complex(36.0),
            prob;
            logarithmic_time_scale = unpack(logarithmic_time_scale, true), kwargs...
        )
    else
        CoreTracker(
            prob.homotopy,
            embed(prob, x₁),
            complex(1.0),
            complex(0.0),
            prob;
            logarithmic_time_scale = unpack(logarithmic_time_scale, false),
            kwargs...
        )
    end
end

# Tracking in affine space
function CoreTracker(
    homotopy::AbstractHomotopy,
    x₁::AbstractVector,
    t₁::Number,
    t₀::Number,
    prob::AbstractProblem;
    patch = nothing,
    corrector::NewtonCorrector = NewtonCorrector(),
    predictor::AbstractPredictor = default_predictor(x₁),
    auto_scaling::Bool = true,
    norm::AbstractNorm = InfNorm(),
    kwargs...
)

    patch === nothing || throw(ArgumentError("You can only pass `patch=$(patch)` if `affine_tracking=false`."))

    options = CoreTrackerOptions(
        ;
        parameter_homotopy = isa(homotopy, ParameterHomotopy), kwargs...
    )

    # We close over the homotopy and its cache to be able to pass things around more easily
    H = HomotopyWithCache(homotopy, Vector(x₁), t₁)

    used_norm = auto_scaling ? WeightedNorm(norm, x₁) : norm
    # We have to make sure that the element type of x is invariant under evaluation
    state = CoreTrackerState(
        H,
        indempotent_x(H, x₁, t₁),
        t₁,
        t₀,
        options,
        nothing,
        used_norm
    )

    pred_cache = cache(predictor, H, state.x, state.ẋ, t₁)
    corr_cache = cache(corrector, H, state.x, t₁)
    dt = Vector{eltype(state.ẋ)}(undef, first(size(H)))

    CoreTracker(H, pred_cache, corr_cache, state, options, dt)
end


# Tracking in Projective Space
function CoreTracker(
    homotopy::AbstractHomotopy,
    x₁::ProjectiveVectors.PVector,
    t₁::Number,
    t₀::Number,
    prob::AbstractProblem;
    patch = has_dedicated_homvars(prob.vargroups) ? EmbeddingPatch() : OrthogonalPatch(),
    corrector::NewtonCorrector = NewtonCorrector(),
    predictor::AbstractPredictor = default_predictor(x₁),
    auto_scaling::Bool = true,
    norm::AbstractNorm = InfNorm(),
    simple_step_size_alg = !isa(patch, EmbeddingPatch),
    kwargs...
)

    options = CoreTrackerOptions(
        ;
        parameter_homotopy = isa(homotopy, ParameterHomotopy), kwargs...
    )

    if homotopy isa PatchedHomotopy
        error("You cannot pass a `PatchedHomotopy` to CoreTracker. Instead pass the homotopy and patch separate.")
    end

    patch_state = state(patch, x₁)

    # We close over the patch state, the homotopy and its cache
    # to be able to pass things around more easily
    H = HomotopyWithCache(PatchedHomotopy(homotopy, patch_state), x₁, t₁)

    if is_global_patch(patch) && auto_scaling
        used_norm = WeightedNorm(norm, x₁)
    else
        used_norm = norm
    end
    # We have to make sure that the element type of x is invariant under evaluation
    ct_state = CoreTrackerState(
        H,
        indempotent_x(H, x₁, t₁),
        t₁,
        t₀,
        options,
        patch_state,
        used_norm
    )

    pred_cache = cache(predictor, H, ct_state.x, ct_state.ẋ, t₁)
    corr_cache = cache(corrector, H, ct_state.x, t₁)
    dt = Vector{eltype(ct_state.ẋ)}(undef, first(size(H)))

    CoreTracker(H, pred_cache, corr_cache, ct_state, options, dt)
end


default_predictor(x::AbstractVector) = Pade21()
# Do not really understand this but Heun doesn't work that great for multi-homogeneous tracking
default_predictor(x::ProjectiveVectors.PVector{T,1}) where {T} = Pade21()
default_predictor(x::ProjectiveVectors.PVector{T,N}) where {T,N} = Euler()

"""
    indempotent_x(H, x₁, t₁)

This returns a vector similar to `x₁` but with an element type which is invariant under evaluation.
"""
function indempotent_x(H::AbstractHomotopy, x₁, t₁)
    u = Vector{Any}(undef, size(H)[1])
    evaluate!(u, H, x₁, t₁)

    if isa(x₁, SVector)
        indem_x = Vector{promote_type(typeof(u[1]), ComplexF64)}(undef, length(x₁))
    else
        indem_x = similar(x₁, promote_type(typeof(u[1]), ComplexF64))
    end
    indem_x .= x₁
    indem_x
end

Base.show(io::IO, C::CoreTracker) =
    print(io, "CoreTracker tracking a path of type $(typeof(C.state.x))")
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::CoreTracker) = x


###########
## Alias ##
###########
const CT = CoreTracker
const CTS = CoreTrackerState
const CTO = CoreTrackerOptions
const CTR = CoreTrackerResult

#######################
## Correct & Tangent ##
#######################

function compute_ẋ!(tracker::CT; update_jacobian::Bool = true)
    @unpack homotopy, state, options, dt = tracker

    if update_jacobian
        jacobian_and_dt!(jacobian(state.jacobian), dt, homotopy, state.x, current_t(state))
        updated!(state.jacobian)
    else
        dt!(dt, homotopy, state.x, current_t(state))
    end
    # flip sign of entries in dt
    @inbounds for i in eachindex(dt)
        dt[i] = -dt[i]
    end
    # solve the linear system to compute ẋ
    LA.ldiv!(state.ẋ, state.jacobian, dt)
    state.ẋ
end


@inline function correct!(
    x̄::AbstractVector,
    tracker::CT,
    x::AbstractVector = tracker.state.x,
    t::Number = current_t(tracker.state)
)

    newton!(
        x̄,
        tracker.homotopy,
        x,
        t,
        tracker.state.jacobian,
        tracker.state.norm,
        tracker.corrector;
        tol = tracker.options.accuracy, max_iters = tracker.options.max_corrector_iters + 1
    )
end

@inline function limit_accuracy!(
    tracker::CT,
    x::AbstractVector = tracker.state.x,
    t::Number = current_t(tracker.state),
    norm::AbstractNorm = tracker.state.norm
)

    limit_acc = limit_accuracy!(
        tracker.state.residual,
        tracker.homotopy,
        x,
        t,
        tracker.state.jacobian,
        norm,
        tracker.corrector;
        accuracy = tracker.options.accuracy
    )
    tracker.state.limit_accuracy = limit_acc
end

####################
## INITIALIZATION ##
####################

embed!(x::ProjectiveVectors.PVector, y) = ProjectiveVectors.embed!(x, y)
embed!(x::AbstractVector, y) = x .= y

function check_start_value(tracker::CT, x, t)
    embed!(tracker.state.x̄, x)
    init!(tracker.state.norm, tracker.state.x̄)
    result = correct!(tracker.state.x̄, tracker, x, t)
    is_converged(result)
end

function check_start_value!(tracker::CT)
    result = correct!(tracker.state.x̄, tracker)
    if is_converged(result)
        tracker.state.x .= tracker.state.x̄
    else
        tracker.state.status = CT_TERMINATED_INVALID_STARTVALUE
    end
    nothing
end


"""
    init!(tracker::CT, x₁, t₁, t₀;
      setup_patch::Bool = tracker.options.update_patch,
      loop::Bool = false,
      check_start_value::Bool = !loop)

Setup `tracker` to track `x₁` from `t₁` to `t₀`. Use this if you want to use the
`tracker` as an iterator.
"""
function init!(
    tracker::CT,
    x₁::AbstractVector,
    t₁,
    t₀;
    setup_patch::Bool = tracker.options.update_patch,
    loop::Bool = false,
    check_start_value::Bool = !loop
)
    @unpack state, predictor, homotopy = tracker

    init!(state, x₁, t₁, t₀, tracker.options, setup_patch, loop)
    check_start_value && check_start_value!(tracker)
    if !loop
        compute_ẋ!(tracker; update_jacobian = true)
        init!(predictor, homotopy, state.x, state.ẋ, current_t(state), state.jacobian)
    end
    tracker
end

function init!(
    state::CTS,
    x₁::AbstractVector,
    t₁,
    t₀,
    options::CTO,
    setup_patch::Bool,
    loop::Bool
)
    state.segment = ComplexSegment(t₁, t₀)
    state.s = 0.0
    state.Δs = min(options.initial_step_size, length(state.segment), options.max_step_size)
    state.Δs_prev = 0.0
    state.status = CT_TRACKING
    embed!(state.x, x₁)
    setup_patch && state.patch !== nothing && init!(state.patch, state.x)
    state.last_step_failed = false
    state.consecutive_successfull_steps = 0
    if !loop
        state.accuracy = 0.0
        state.ω = 1.0
        state.accepted_steps = state.rejected_steps = 0
        init!(state.norm, state.x)
        init!(state.jacobian)
        state.residual .= 0.0
        state.steps_jacobian_info_update = 0
    end

    state
end

@deprecate setup!(tracker::CoreTracker, x₁, t₁, t₀; kwargs...) init!(
    tracker,
    x₁,
    t₁,
    t₀;
    kwargs...
)

##############
## TRACKING ##
##############

"""
    track(tracker, x₁, t₁=1.0, t₀=0.0; options...)::CTR

Track a value `x₁` from `t₁` to `t₀` using the given `CoreTracker` `tracker`.
This returns a `CoreTrackerResult`. This modifies `tracker`.
See [`track!`](@ref) for the possible options.
"""
function track(tracker::CT, x₁::AbstractVector, t₁ = 1.0, t₀ = 0.0; kwargs...)
    track!(tracker, x₁, t₁, t₀; kwargs...)
    CoreTrackerResult(tracker)
end

"""
    track!(tracker, x₁, t₁=1.0, t₀=0.0; setup_patch=true, check_start_value=true, loop::Bool=false)::CTR

Track a value `x₁` from `t₁` to `t₀` using the given `CoreTracker` `tracker`.
Returns one of the enum values of `CoreTrackerStatus` indicating the status.
Check [`is_success`](@ref) and [`is_terminated`](@ref) to test for the status.

If `setup_patch` is `true` then [`init!`](@ref) is called at the beginning
of the tracking.

    track!(x₀, tracker, x₁, t₁=1.0, t₀=0.0; options...)::CTStatus

Additionally also stores the result in `x₀` if the tracking was successfull.
"""
function track!(
    x₀::AbstractVector,
    tracker::CT,
    x₁::AbstractVector,
    t₁ = 1.0,
    t₀ = 0.0;
    setup_patch::Bool = tracker.options.update_patch,
    loop::Bool = false,
    check_start_value::Bool = !loop,
    debug::Bool = false
)

    _track!(tracker, x₁, t₁, t₀, setup_patch, loop, check_start_value, debug)
    retcode = status(tracker)
    if is_success(status(tracker))
        x₀ .= current_x(tracker)
    end
    retcode
end
function track!(
    tracker::CT,
    x₁,
    t₁ = 1.0,
    t₀ = 0.0;
    setup_patch::Bool = tracker.options.update_patch,
    loop::Bool = false,
    check_start_value::Bool = !loop,
    debug::Bool = false
)

    _track!(tracker, x₁, t₁, t₀, setup_patch, loop, check_start_value, debug)
end

@inline function _track!(
    tracker::CT,
    x₁::AbstractVector,
    t₁::Number,
    t₀::Number,
    setup_patch::Bool,
    loop::Bool,
    check_start_value::Bool,
    debug::Bool
)

    @unpack state = tracker

    t₁ == t₀ && return (state.status = CT_SUCCESS)

    init!(
        tracker,
        x₁,
        t₁,
        t₀;
        setup_patch = setup_patch, check_start_value = check_start_value, loop = loop
    )

    while is_tracking(state.status)
        step!(tracker, debug)
        check_terminated!(tracker)
    end

    is_success(state.status) && refine!(tracker)

    state.status
end


"""
    step!(tracker::CT)

Perform one step of the path tracker. This tries to move forward by `Δt`, i.e.,
obtaining an approximation of `x(t+Δt)`.
Returns `true` if the step was successfull.
"""
function step!(tracker::CT, debug::Bool = false)
    @unpack homotopy, predictor, state, options = tracker
    @unpack x, x̂, x̄, ẋ = state

    # @unpack currently doesn't work with getproperty
    t, Δt = current_t(state), current_Δt(state)

    debug && printstyled(
        "\nt: ",
        round(t, sigdigits = 5),
        " Δt: ",
        round(Δt; sigdigits = 5),
        "\n";
        color = :yellow
    )

    # Use the current approximation of x(t) to obtain estimate
    # x̂ ≈ x(t + Δt) using the provided predictor
    predict!(x̂, predictor, homotopy, x, t, Δt, ẋ, state.jacobian)

    # Correct the predicted value x̂ to obtain x̄.
    # If the correction is successfull we have x̄ ≈ x(t+Δt).
    result = correct!(x̄, tracker, x̂, t + Δt)

    debug && show(result)

    if is_converged(result)
        # move forward
        x .= x̄
        state.s += state.Δs
        state.accuracy = result.accuracy

        # update step size
        state.Δs_prev = state.Δs
        update_stepsize!(tracker, result)

        # If we use dynamic affine charts, now we update to a new one
        if state.patch !== nothing && options.update_patch
            changepatch!(state.patch, x)
        end

        # If we use a weighted norm, now we update this
        update!(state.norm, x)

        # TODO: THIS SHOULDN'T BE DONE IN EVERY STEP!
        # Update informations regarding limit accuracy
        limit_accuracy!(tracker)

        if debug
            println("limit_acc: ", round(state.limit_accuracy; sigdigits = 5))
        end

        # Compute the new tangent ẋ(t)
        compute_ẋ!(tracker)

        # tell the predictors about the new derivative if they need to update something
        update!(predictor, homotopy, x, ẋ, t + Δt, state.jacobian)

        # Update other state
        state.accepted_steps += 1
        state.steps_jacobian_info_update += 1
    else
        # Step failed, so we have to try with a new (smaller) step size
        update_stepsize!(tracker, result)

        # Update other state
        state.rejected_steps += 1
    end
    state.last_step_failed = is_converged(result)
end

## Step size update

function update_stepsize!(tracker::CT, result::NewtonCorrectorResult)
    @unpack state, options = tracker

    if !isnan(result.ω₀)
        # Use result.ω₀ in the converged case to avoid wrong values of ω due to
        # hitting the accuracy limit.
        # In the divergend case this should not happen
        state.ω = is_converged(result) ? result.ω₀ : result.ω
    end

    Δs = let
        if options.simple_step_size_alg
            simple_step_size_alg!(state, options, result)
        else
            adaptive_step_size_alg!(state, options, result, order(tracker.predictor))
        end
    end

    # Make sure to not overshoot.
    if !is_min_step_size(Δs, state.s, options)
        state.status = CT_TERMINATED_STEP_SIZE_TOO_SMALL
    else
        state.Δs = min(length(state.segment) - state.s, Δs, options.max_step_size)
    end

    nothing
end

function is_min_step_size(Δs, s, options::CTO)
    if options.logarithmic_time_scale
        exp(-s) - exp(-(s + Δs)) ≥ options.min_step_size
    else
        Δs ≥ options.min_step_size
    end
end


## Simple step size algorithm

function simple_step_size_alg!(state::CTS, options::CTO, result::NewtonCorrectorResult)

    if !is_converged(result)
        state.consecutive_successfull_steps = 0
        return state.Δs / 2
    end

    state.consecutive_successfull_steps = (state.consecutive_successfull_steps + 1) % 5
    state.consecutive_successfull_steps == 0 ? 2state.Δs : state.Δs
end

## Adaptive step size algorithm

"Use the adaptive step size control described in https://arxiv.org/abs/1902.02968"
function adaptive_step_size_alg!(
    state::CTS,
    options::CTO,
    result::NewtonCorrectorResult,
    ord::Int
)
    @unpack s, ω = state

    d_x̂_x = result.iters > 1 ? state.norm(state.x̂, state.x) : result.norm_Δx₀
    Δx₀ = result.norm_Δx₀

    if is_converged(result)
        # This is to handle the edge case that g(δ_N_ω) >> (ω * d_x̂_x̄) > 0 but
        # at the same time δ_N_ω < eps(). Since then g(δ_N_ω) = 0
        δ_N_ω = max(δ(options, ω, 0.5), 1e-15)
        Δs = nthroot(g(δ_N_ω) / (ω * d_x̂_x), ord) * state.Δs
    else
        δ_N_ω = δ(options, ω, 0.5)
        ω_η = 0.5ω * Δx₀
        if δ_N_ω < ω_η
            Δs = min(nthroot(g(δ_N_ω) / g(ω_η), ord), 0.9) * state.Δs
        else
            Δs = 0.5state.Δs
        end
    end

    Δs
end

g(Θ::Real) = sqrt(1. + 4Θ) - 1.0
δ(opts::CTO, ω::Real, μ::Real) = min(√(0.5ω) * τ(opts, μ), 0.25)
function τ(opts::CTO, μ::Real)
    # most common case: 2(2 - 0.5) == 3
    if opts.max_corrector_iters == 2 && μ == 0.5
        cbrt(opts.accuracy)
    else
        opts.accuracy^(1 / (2(opts.max_corrector_iters - μ)))
    end
end


function check_terminated!(tracker::CT)
    Δs = length(tracker.state.segment) - tracker.state.s
    if Δs < eps(length(tracker.state.segment))
        tracker.state.status = CT_SUCCESS
        tracker.state.s = length(tracker.state.segment)
    elseif steps(tracker) ≥ tracker.options.max_steps
        tracker.state.status = CT_TERMINATED_MAX_ITERS
    end
    nothing
end

function refine!(tracker::CT, accuracy = tracker.options.refinement_accuracy)
    # if tracker.state.accuracy < accuracy
    #     return
    # end
    # result = correct!(tracker.state.x̄, tracker;
    #         accuracy = accuracy,
    #         max_iters = tracker.options.max_refinement_iters)
    # if is_converged(result)
    #     tracker.state.x .= tracker.state.x̄
    #     tracker.state.accuracy = result.accuracy
    # end
    nothing
end

function Base.iterate(tracker::CoreTracker, state::Union{Nothing,Int}=nothing)
    state === nothing && return tracker, 1

    if is_tracking(tracker.state.status)
        step!(tracker)
        check_terminated!(tracker)

        if is_success(tracker.state.status)
            refine!(tracker)
        end

        tracker, state + 1
    else
        nothing
    end
end


"""
    affine_tracking(tracker)


Returns `true` if the path tracker tracks in affine space. If `false` then then the path tracker
tracks in some affine chart of the projective space.
"""
affine_tracking(tracker::CT) = !isa(tracker.state.x, ProjectiveVectors.PVector)

#################
## Query State ##
#################
"""
    current_t(tracker::CT)

Current `t`.
"""
current_t(tracker::CT) = current_t(tracker.state)
current_t(state::CTS) = state.segment[state.s]

"""
 current_Δt(tracker::CT)

Current step_size `Δt`.
"""
current_Δt(tracker::CT) = current_Δt(tracker.state)
current_Δt(state::CTS) = state.segment[state.Δs] - state.segment.start

"""
    steps(tracker::CT)

Current number of steps.
"""
steps(tracker::CT) = steps(tracker.state)
steps(state::CTS) = state.accepted_steps + state.rejected_steps

@deprecate iters(tracker::CT) steps(tracker)


"""
    status(tracker::CT)

Current status.
"""
status(tracker::CT) = status(tracker.state)
status(state::CTS) = state.status

"""
    current_x(tracker::CT)

Return the current value of `x`.
"""
current_x(tracker::CT) = current_x(tracker.state)
current_x(state::CTS) = state.x


"""
    LinearAlgebra.cond(tracker::CT)

Returns the currently computed approximation of the condition number of the
Jacobian.
"""
cond(tracker::CT) = LinearAlgebra.cond(tracker.state)
cond(state::CTS) = state.jacobian.cond


"""
    digits_lost(tracker::CT)

Returns the currently computed approximation of the number of digits lost
during the linear system solving in Newton's method.
"""
digits_lost(tracker::CT) = digits_lost(tracker.state)
digits_lost(state::CTS) = unpack(state.jacobian.digits_lost, 0.0)


"""
    norm(tracker::CT)

Returns the norm used to compute distances during the path tracking.
"""
norm(tracker::CT) = norm(tracker.state)
norm(state::CTS) = state.norm

"""
    options(tracker::CT)

Returns the options used in the tracker.
"""
options(tracker::CT) = tracker.options

##################
# Modify options #
##################
"""
    accuracy(tracker::CT)

Current accuracy.
"""
accuracy(tracker::CT) = tracker.options.accuracy

"""
    set_accuracy!(tracker::CT, accuracy, update_max_lost_digits=true)

Set the current accuracy to `accuracy`. If `update_max_lost_digits` is `true` then
the setting `max_lost_digits` will be updated to the default setting.
"""
function set_accuracy!(tracker::CT, accuracy, update_max_lost_digits::Bool = true)
    @unpack options = tracker
    options.accuracy = accuracy
    if update_max_lost_digits
        options.max_lost_digits = default_max_lost_digits(options.precision, accuracy)
    end
    tracker
end

"""
    refinement_accuracy(tracker::CT)

Current refinement accuracy.
"""
refinement_accuracy(tracker::CT) = tracker.options.refinement_accuracy

"""
    set_max_refinement_iters!(tracker::CT, accuracy)

Set the current refinement accuracy to `accuracy`.
"""
set_refinement_accuracy!(T::CT, accuracy) = T.options.refinement_accuracy = accuracy

"""
    max_refinement_iters(tracker::CT)

Current refinement max_steps.
"""
max_refinement_iters(T::CT) = T.options.max_refinement_iters

"""
    set_max_refinement_iters!(tracker::CT, n)

Set the current refinement max_steps to `n`.
"""
set_max_refinement_iters!(T::CT, n) = T.options.max_refinement_iters = n

"""
    max_corrector_iters(tracker::CT)

Current correction max_steps.
"""
max_corrector_iters(T::CT) = T.options.max_corrector_iters

"""
 set_max_corrector_iters!(tracker::CT, n)

Set the correction max_steps to `n`.
"""
set_max_corrector_iters!(T::CT, n) = T.options.max_corrector_iters = n

"""
    max_step_size (tracker::CT)

Current maximal step size.
"""
max_step_size(T::CT) = T.options.max_step_size

"""
    set_max_corrector_iters!(tracker::CT, Δs)

Set the maximal step size to `Δs`.
"""
set_max_step_size!(T::CT, Δs) = T.options.max_step_size = Δs

################
# PathIterator #
################
struct PathIterator{Tracker<:CoreTracker}
    tracker::Tracker
    t_real::Bool
end
Base.IteratorSize(::Type{<:PathIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:PathIterator}) = Base.HasEltype()

"""
    iterator(tracker::CT, x₁, t₁=1.0, t₀=0.0; affine=true)

Prepare a tracker to make it usable as a (stateful) iterator. Use this if you want to inspect a specific
path. In each iteration the tuple `(x,t)` is returned.
If `affine == true` then `x` is the affine solution (internally we compute in projective space).

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
println("Success: ", status(tracker) == CT_SUCCESS)
```
"""
function iterator(tracker::CT, x₁, t₁ = 1.0, t₀ = 0.0; kwargs...)
    init!(tracker, x₁, t₁, t₀; kwargs...)
    PathIterator(tracker, typeof(t₁ - t₀) <: Real)
end

function current_x_t(iter::PathIterator)
    x = current_x(iter.tracker)
    t = current_t(iter.tracker)
    (copy(x), iter.t_real ? real(t) : t)
end

function Base.iterate(iter::PathIterator, state = nothing)
    state === nothing && return current_x_t(iter), 1
    iter.tracker.state.status != CT_TRACKING && return nothing

    t = current_t(iter.tracker)
    step_done = false
    while !step_done && is_tracking(iter.tracker.state.status)

        step!(iter.tracker)
        check_terminated!(iter.tracker)

        if is_success(iter.tracker.state.status)
            refine!(iter.tracker)
        end

        step_done = current_t(iter.tracker) != t
    end
    current_x_t(iter), state + 1
end


#############################
# Convencience constructors #
#############################
"""
    coretracker_startsolutions(args...; kwargs...)

Construct a [`CoreTracker`](@ref) and `startsolutions` in the same way `solve`
does it. This also takes the same input arguments as `solve`. This is convenient if you want
to investigate single paths.
"""
function coretracker_startsolutions(args...; kwargs...)
    supported, rest = splitkwargs(kwargs, problem_startsolutions_supported_keywords)
    prob, startsolutions = problem_startsolutions(args...; supported...)
    tracker = CoreTracker(prob, start_solution_sample(startsolutions); rest...)
    (tracker = tracker, startsolutions = startsolutions)
end

"""
    coretracker(args...; kwargs...)

Construct a [`CoreTracker`](@ref) in the same way `solve`
does it. This also takes the same input arguments as `solve` with the exception that you do not need to specify startsolutions.
This is convenient if you want to investigate single paths.

## Examples

### Obtain single solution
We want to construct a path tracker to track a parameterized system `f` with parameters `p`
from the parameters `a` to `b`.
```julia
tracker = CoreTracker(f, parameters=p, p₁=a, p₀=b)
```
You then can obtain a single solution at `b` by using
```julia
x_b = track(tracker, x_a).x
```

### Trace a path
To trace a path you can use the [`iterator`](@ref) method.

```julia
tracker = CoreTracker(f, parameters=p, p₁=a, p₀=b, max_step_size =0.01)
for (x, t) in iterator(tracker, x₁)
@show (x,t)
end
```

If we want to guarantee smooth traces we can limit the maximal step size.
```julia
tracker = CoreTracker(f, parameters=p, p₁=a, p₀=b, max_step_size =0.01)
for (x, t) in iterator(tracker, x₁)
@show (x,t)
end
```
"""
function coretracker(args...; kwargs...)
    tracker, _ = coretracker_startsolutions(args...; kwargs...)
    tracker
end
