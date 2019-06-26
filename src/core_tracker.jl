export CoreTracker, CoreTrackerResult, CoreTrackerStatus, CoreTrackerOptions,
        coretracker, coretracker_startsolutions, affine_tracking,
        track, track!, setup!, iterator,
        current_x, current_t, current_Δt, iters, status, accuracy, max_corrector_iters,
        max_step_size , refinement_accuracy, max_refinement_iters,
        set_accuracy!, set_max_corrector_iters!, set_refinement_accuracy!,
        set_max_refinement_iters!, set_max_step_size!, digits_lost, options

const coretracker_supported_keywords = [:corrector, :predictor, :patch,
    :initial_step_size, :min_step_size , :max_step_size,
    :accuracy, :refinement_accuracy, :max_corrector_iters, :max_refinement_iters,
    :max_steps, :simple_step_size_alg, :auto_scaling, :terminate_ill_conditioned,
    :log_transform, :precision, :steps_jacobian_info_update,
    :max_lost_digits]


####################
# CoreTrackerState #
####################

module CoreTrackerStatus

    @doc """
        CoreTrackerStatus.states

    The possible states the coretracker can achieve are

    * `CoreTrackerStatus.success`
    * `CoreTrackerStatus.tracking`
    * `CoreTrackerStatus.terminated_maximal_iterations`
    * `CoreTrackerStatus.terminated_invalid_startvalue`
    * `CoreTrackerStatus.terminated_step_size_too_small`
    * `CoreTrackerStatus.terminated_singularity`
    * `CoreTrackerStatus.terminated_ill_conditioned`
    """
    @enum states begin
        success
        tracking
        terminated_maximal_iterations
        terminated_invalid_startvalue
        terminated_step_size_too_small
        terminated_singularity
        terminated_ill_conditioned
    end
end


####################
# CoreTrackerResult #
####################

"""
     CoreTrackerResult{V<:AbstractVector}

Containing the result of a tracked path. The fields are
* `returncode::CoreTrackerStatus.states` If the tracking was successfull then it is `CoreTrackerStatus.success`.
* `x::V` The result.
* `t::ComplexF64` The `t` when the path tracker stopped.
* `accuracy::Float64`: The estimated accuracy of `x`.
"""
struct CoreTrackerResult{V <: AbstractVector}
     returncode::CoreTrackerStatus.states
     x::V
     t::ComplexF64
     accuracy::Float64
     accepted_steps::Int
     rejected_steps::Int
end

function CoreTrackerResult(tracker)
    state = tracker.state
     CoreTrackerResult(state.status,
          copy(state.x), state.segment[state.s],
          state.accuracy,
          state.accepted_steps,
          state.rejected_steps)
end

Base.show(io::IO, result::CoreTrackerResult) = print_fieldnames(io, result)
Base.show(io::IO, ::MIME"application/prs.juno.inline", result::CoreTrackerResult) = result


##################
## PATH TRACKER ##
##################

"""
    AutoScalingOptions(;scale_min=0.01,
                     scale_abs_min=min(scale_min^2, 200 * sqrt(eps()),
                     scale_max=1.0 / eps() / sqrt(2)

Parameters for the auto scaling of the variables.
"""
struct AutoScalingOptions
    scale_min::Float64
    scale_abs_min::Float64
    scale_max::Float64
end
function AutoScalingOptions(;scale_min=0.01,
                             scale_abs_min=min(scale_min^2, 200 * sqrt(eps())),
                             scale_max=1.0 / eps() / sqrt(2))
    AutoScalingOptions(scale_min, scale_abs_min, scale_max)
end
Base.show(io::IO, opts::AutoScalingOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::AutoScalingOptions) = opts

mutable struct StepSizeModel
    ω::Float64
    expected_Δx₀::Float64
end
StepSizeModel() = StepSizeModel(NaN, NaN)

function reset!(model::StepSizeModel)
    model.ω = model.expected_Δx₀ = NaN
    model
end
Base.show(io::IO, opts::StepSizeModel) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::StepSizeModel) = opts

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
    min_step_size ::Float64
    max_step_size ::Float64
    simple_step_size_alg::Bool
    update_patch::Bool
    max_lost_digits::Float64
    auto_scaling::Bool
    auto_scaling_options::AutoScalingOptions
    terminate_ill_conditioned::Bool
    precision::PrecisionOption
    steps_jacobian_info_update::Int
end

function CoreTrackerOptions(;
    accuracy=1e-7,
    auto_scaling=true,
    auto_scaling_options=AutoScalingOptions(),
    initial_step_size=0.1,
    max_corrector_iters::Int=2,
    max_refinement_iters=5,
    max_step_size=Inf,
    min_step_size=1e-14,
    parameter_homotopy=false,
    max_steps=parameter_homotopy ? 10_000 : 1_000,
    precision::PrecisionOption=PRECISION_FIXED_64,
    max_lost_digits=default_max_lost_digits(precision, accuracy),
    refinement_accuracy=1e-8,
    simple_step_size_alg=false,
    steps_jacobian_info_update::Int=1,
    terminate_ill_conditioned::Bool=true,
    update_patch=true)

    CoreTrackerOptions(accuracy, max_corrector_iters, refinement_accuracy,
            max_refinement_iters, max_steps, initial_step_size, min_step_size,
            max_step_size, simple_step_size_alg, update_patch, max_lost_digits,
            auto_scaling, auto_scaling_options, terminate_ill_conditioned,
            precision, steps_jacobian_info_update)
end

Base.show(io::IO, opts::CoreTrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::CoreTrackerOptions) = opts

function default_max_lost_digits(prec::PrecisionOption, accuracy::Float64)
    if prec == PRECISION_FIXED_64
        # maximal_digits_available - digits necessary - buffer
        -log10(eps()) + log10(accuracy) + 1
    else
        # if higher precision is available we will more like be constrained
        # by the fact that the jacobian cannot be too ill-conditioned
        min(-log10(eps()) - 3, -log10(eps(Double64)) + log10(accuracy) - 3)
    end
end

####################
# CoreTrackerState #
####################

mutable struct CoreTrackerState{T, AV<:AbstractVector{T}, MaybePatchState <: Union{AbstractAffinePatchState, Nothing}, IP}
    x::AV # current x
    x̂::AV # last prediction
    x̄::AV # canidate for new x
    ẋ::Vector{T} # derivative at current x
    inner_product::IP
    step_size::StepSizeModel
    segment::ComplexSegment
    s::Float64 # current step length (0 ≤ s ≤ length(segment))
    Δs::Float64 # current step size
    Δs_prev::Float64 # previous step size
    accuracy::Float64
    jacobian::Jacobian{T}
    steps_jacobian_info_update::Int
    status::CoreTrackerStatus.states
    patch::MaybePatchState
    accepted_steps::Int
    rejected_steps::Int
    last_step_failed::Bool
    consecutive_successfull_steps::Int
end

function CoreTrackerState(H, x₁::AbstractVector, t₁, t₀, options::CoreTrackerOptions,
                           patch::Union{Nothing, AbstractAffinePatchState}=nothing,
                           inner_product::AbstractInnerProduct=EuclideanIP())
    if isa(x₁, SVector)
        x = Vector(x₁)
    else
        x = copy(x₁)
    end
    x̂, x̄ = copy(x), copy(x)
    ẋ = Vector(x)
    step_size_model = StepSizeModel()
    segment = ComplexSegment(promote(complex(t₁), complex(t₀))...)
    s = 0.0
    Δs = convert(Float64, min(options.initial_step_size, length(segment), options.max_step_size ))
    Δs_prev = 0.0
    accuracy = 0.0
    Jac = Jacobian(jacobian(H, x, t₁))
    steps_jacobian_info_update = 0
    accepted_steps = rejected_steps = 0
    status = CoreTrackerStatus.tracking
    last_step_failed = false
    consecutive_successfull_steps = 0
    CoreTrackerState(x, x̂, x̄, ẋ, inner_product, step_size_model, segment, s, Δs, Δs_prev, accuracy,
        Jac, steps_jacobian_info_update, status, patch,
        accepted_steps, rejected_steps, last_step_failed, consecutive_successfull_steps)
end

Base.show(io::IO, state::CoreTrackerState) = print_fieldnames(io, state)
Base.show(io::IO, ::MIME"application/prs.juno.inline", state::CoreTrackerState) = state

function reset!(state::CoreTrackerState, x₁::AbstractVector, t₁, t₀, options::CoreTrackerOptions, setup_patch::Bool, loop::Bool)
    state.segment = ComplexSegment(promote(t₁, t₀)...)
    state.s = 0.0
    state.Δs = min(options.initial_step_size, length(state.segment), options.max_step_size )
    state.Δs_prev = 0.0
    state.status = CoreTrackerStatus.tracking
    embed!(state.x, x₁)
    setup_patch && state.patch !== nothing && setup!(state.patch, state.x)
    state.last_step_failed = false
    state.consecutive_successfull_steps = 0
    if !loop
        state.accuracy = 0.0
        state.accepted_steps = state.rejected_steps = 0
        reset!(state.step_size)

        if options.auto_scaling
            init_auto_scaling!(state.inner_product, state.x, options.auto_scaling_options)
        end
        reset!(state.jacobian)
        state.steps_jacobian_info_update = 0
    end

    state
end

embed!(x::ProjectiveVectors.PVector, y) = ProjectiveVectors.embed!(x, y)
embed!(x::AbstractVector, y) = x .= y

function init_auto_scaling!(ip::WeightedIP, x::AbstractVector, opts::AutoScalingOptions)
    point_norm = euclidean_norm(x)
    w = ip.weight
    for i in 1:length(w)
        wᵢ = abs(x[i])
        if wᵢ < opts.scale_min * point_norm
            wᵢ = opts.scale_min * point_norm
        elseif wᵢ > opts.scale_max * point_norm
            wᵢ = opts.scale_max * point_norm
        end
        w[i] = max(wᵢ, opts.scale_abs_min)
    end
    nothing
end
init_auto_scaling!(ip::EuclideanIP, x::AbstractVector, opts::AutoScalingOptions) = nothing

####################
# CoreTrackerCache #
####################
mutable struct CoreTrackerCache{H<:HomotopyWithCache, P<:AbstractPredictorCache,
             C<:AbstractCorrectorCache, T}
    homotopy::H
    predictor::P
    corrector::C
    out::Vector{T}
    r::Vector{T}
end
function CoreTrackerCache(H::HomotopyWithCache, predictor, corrector, state::CoreTrackerState)
    t = state.segment[state.s]
    pcache = cache(predictor, H, state.x, state.ẋ, t)
    ccache = cache(corrector, H, state.x, t)
    out = H(state.x, t)
    r = copy(out)
    CoreTrackerCache(H, pcache, ccache, out, r)
end


###############
# CoreTracker #
###############
"""
    CoreTracker(problem::AbstractProblem, x₁; kwargs...)

Construct a `CoreTracker` from the given `problem` to track elements of type `x₁`.
The path is tracked using a predictor-corrector scheme. The recommended methods to construct
a `CoreTracker` are [`coretracker`](@ref) and [`coretracker_startsolutions`](@ref).
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
struct CoreTracker{H<:AbstractHomotopy,
    Predictor<:AbstractPredictor,
    Corrector<:AbstractCorrector,
    MaybePatch<:Union{Nothing, AbstractAffinePatch},
    S<:CoreTrackerState,
    C<:CoreTrackerCache}
    # these are fixed
    homotopy::H
    predictor::Predictor
    corrector::Corrector
    affine_patch::MaybePatch
    # these are mutable
    state::S
    options::CoreTrackerOptions
    cache::C
end

function CoreTracker(prob::AbstractProblem, x₁; kwargs...)
    CoreTracker(prob.homotopy, embed(prob, x₁), complex(1.0), complex(0.0), prob; kwargs...)
end

# Tracking in Projective Space
function CoreTracker(homotopy::AbstractHomotopy, x₁::ProjectiveVectors.PVector, t₁, t₀, prob::AbstractProblem;
    patch=has_dedicated_homvars(prob.vargroups) ? EmbeddingPatch() : OrthogonalPatch(),
    corrector::AbstractCorrector=NewtonCorrector(),
    predictor::AbstractPredictor=default_predictor(x₁),
    log_transform=false,
    simple_step_size_alg=!isa(patch, EmbeddingPatch), kwargs...)

    options = CoreTrackerOptions(;
                    parameter_homotopy=isa(homotopy, ParameterHomotopy),
                    kwargs...)

    if homotopy isa PatchedHomotopy
        error("You cannot pass a `PatchedHomotopy` to CoreTracker. Instead pass the homotopy and patch separate.")
    end

    H = log_transform ? LogHomotopy(homotopy) : homotopy
    patch_state = state(patch, x₁)
    # We close over the patch state, the homotopy and its cache
    # to be able to pass things around more easily
    HC = HomotopyWithCache(PatchedHomotopy(H, patch_state), x₁, t₁)
    if is_global_patch(patch)
        inner_product = WeightedIP(x₁)
    else
        inner_product = EuclideanIP()
    end
    # We have to make sure that the element type of x is invariant under evaluation
    tracker_state = CoreTrackerState(HC, indempotent_x(HC, x₁, t₁), t₁, t₀, options, patch_state, inner_product)
    cache = CoreTrackerCache(HC, predictor, corrector, tracker_state)

    CoreTracker(H, predictor, corrector, patch, tracker_state, options, cache)
end

# Tracking in affine space
function CoreTracker(homotopy::AbstractHomotopy, x₁::AbstractVector, t₁, t₀, prob::AbstractProblem;
    patch=nothing,
    corrector::AbstractCorrector=NewtonCorrector(),
    predictor::AbstractPredictor=default_predictor(x₁),
    log_transform=false, kwargs...)

    if patch !== nothing
        throw(ArgumentError("You can only pass `patch=$(patch)` if `affine_tracking=false`."))
    end

    options = CoreTrackerOptions(;parameter_homotopy=isa(homotopy, ParameterHomotopy), kwargs...)

    H = log_transform ? LogHomotopy(homotopy) : homotopy
    # We close over the patch state, the homotopy and its cache
    # to be able to pass things around more easily
    HC = HomotopyWithCache(H, Vector(x₁), t₁)

    inner_product = options.auto_scaling ? WeightedIP(x₁) : EuclideanIP()
    # We have to make sure that the element type of x is invariant under evaluation
    tracker_state = CoreTrackerState(HC, indempotent_x(HC, x₁, t₁), t₁, t₀, options, nothing, inner_product)
    cache = CoreTrackerCache(HC, predictor, corrector, tracker_state)

    CoreTracker(H, predictor, corrector, nothing, tracker_state, options, cache)
end

default_predictor(x::AbstractVector) = Heun()
# Do not really understand this but Heun doesn't work that great for multi-homogeneous tracking
default_predictor(x::ProjectiveVectors.PVector{T,1}) where {T} = Heun()
default_predictor(x::ProjectiveVectors.PVector{T,N}) where {T, N} = Euler()

"""
    indempotent_x(H, x₁, t₁)

This returns a vector similar to `x₁` but with an element type which is invariant under evaluation.
"""
function indempotent_x(H, x₁, t₁)
    u = Vector{Any}(undef, size(H)[1])
    evaluate!(u, H, x₁, t₁)

    if isa(x₁, SVector)
        indem_x = Vector{promote_type(typeof(u[1]), ComplexF64)}(undef, length(x₁))
    else
        indem_x = similar(x₁, promote_type(typeof(u[1]), ComplexF64))
    end
    indem_x .= x₁
end

Base.show(io::IO, C::CoreTracker) = print(io, "CoreTracker tracking a path of type $(typeof(C.state.x))")
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::CoreTracker) = x

##############
## TRACKING ##
##############

"""
    track(tracker, x₁, t₁=1.0, t₀=0.0; options...)::CoreTrackerResult

Track a value `x₁` from `t₁` to `t₀` using the given `CoreTracker` `tracker`.
This returns a `CoreTrackerResult`. This modifies `tracker`.
See [`track!`](@ref) for the possible options.
"""
function track(tracker::CoreTracker, x₁::AbstractVector, t₁=1.0, t₀=0.0; kwargs...)
     track!(tracker, x₁, t₁, t₀; kwargs...)
     CoreTrackerResult(tracker)
end

"""
     track!(tracker, x₁, t₁=1.0, t₀=0.0; setup_patch=true, checkstartvalue=true, loop::Bool=false)::CoreTrackerResult

Track a value `x₁` from `t₁` to `t₀` using the given `CoreTracker` `tracker`.
Returns one of the enum values of `CoreTrackerStatus.states` indicating the status.
If the tracking was successfull it is `CoreTrackerStatus.success`.
If `setup_patch` is `true` then [`setup!`](@ref) is called at the beginning
of the tracking.

    track!(x₀, tracker, x₁, t₁=1.0, t₀=0.0; options...)::CoreTrackerStatus.states

Additionally also stores the result in `x₀` if the tracking was successfull.
"""
function track!(x₀, tracker::CoreTracker, x₁, t₁=1.0, t₀=0.0; setup_patch::Bool=tracker.options.update_patch,
            loop::Bool=false, checkstartvalue::Bool=!loop)
     _track!(tracker, x₁, t₁, t₀, setup_patch, checkstartvalue, loop)
     retcode = status(tracker)
     if retcode == CoreTrackerStatus.success
         x₀ .= current_x(tracker)
     end
     retcode
end
@inline function track!(tracker::CoreTracker, x₁, t₁=1.0, t₀=0.0;
        setup_patch::Bool=tracker.options.update_patch,
        checkstartvalue::Bool=true,
        loop::Bool=false)
    _track!(tracker, x₁, t₁, t₀, setup_patch, checkstartvalue, loop)
end

function _track!(tracker::CoreTracker, x₁, t₁, t₀,
                    setup_patch::Bool, checkstartvalue::Bool, loop::Bool)
    if t₁ == t₀
        tracker.state.status = CoreTrackerStatus.success
        return tracker.state.status
    end

    setup!(tracker, x₁, t₁, t₀, setup_patch, checkstartvalue, loop)

    while tracker.state.status == CoreTrackerStatus.tracking
        step!(tracker)
        check_terminated!(tracker)
    end
    if tracker.state.status == CoreTrackerStatus.success
        refine!(tracker)
    end

    tracker.state.status
end

"""
    setup!(coretracker, x₁, t₁=1.0, t₀=0.0, setup_patch=coretracker.options.update_patch, checkstartvalue=true, compute_ẋ=true)

Setup `coretracker` to track `x₁` from `t₁` to `t₀`. Use this if you want to use the
coretracker as an iterator.
"""
function setup!(tracker::CoreTracker, x₁::AbstractVector, t₁=1.0, t₀=0.0, setup_patch=tracker.options.update_patch, checkstartvalue=true, loop::Bool=false)
    @unpack state, cache = tracker

    try
        reset!(state, x₁, t₁, t₀, tracker.options, setup_patch, loop)
        reset!(cache.predictor, state.x, t₁)
        checkstartvalue && checkstartvalue!(tracker)
        if !loop
            compute_ẋ!(state, cache, tracker.options)
            setup!(cache.predictor, cache.homotopy, state.x, state.ẋ, current_t(state), state.jacobian)
        end
    catch err
        if !(err isa LinearAlgebra.SingularException)
            rethrow(err)
        end
        tracker.state.status = CoreTrackerStatus.terminated_singularity
    end
    tracker
end

function checkstartvalue!(tracker::CoreTracker)
    result = correct!(tracker.state.x̄, tracker; update_jacobian_infos=true)
    if isconverged(result)
        tracker.state.x .= tracker.state.x̄
    else
        tracker.state.status = CoreTrackerStatus.terminated_invalid_startvalue
    end
    nothing
end

function compute_ẋ!(state, cache, options::CoreTrackerOptions)
    @inbounds jacobian_and_dt!(state.jacobian.J, cache.out, cache.homotopy, state.x, current_t(state))
    # apply row scaling to J and compute factorization
    updated_jacobian!(state.jacobian)

    @inbounds for i in eachindex(cache.out)
        cache.out[i] = -cache.out[i]
    end
    solve!(state.ẋ, state.jacobian, cache.out)
    nothing
end


@inline function correct!(x̄, tracker::CoreTracker,
        x=tracker.state.x,
        t=tracker.state.segment[tracker.state.s],
        norm=tracker.state.inner_product;
        accuracy::Float64=tracker.options.accuracy,
        max_iters::Int=tracker.options.max_corrector_iters,
        precision::PrecisionOption=tracker.options.precision,
        update_jacobian_infos::Bool=false,
        use_qr::Bool=false)

    correct!(x̄, tracker.corrector, tracker.cache.corrector,
             tracker.cache.homotopy, x, t, norm, tracker.state.jacobian,
             accuracy, max_iters, tracker.state.step_size;
             precision=precision,
             update_jacobian_infos=update_jacobian_infos, use_qr=use_qr)
end

function step!(tracker::CoreTracker)
    @unpack state, cache, options = tracker
    @unpack x, x̂, x̄, ẋ = state
    H = cache.homotopy

    try
        t, Δt = current_t(state), current_Δt(state)
        predict!(x̂, tracker.predictor, cache.predictor, H, x, t, Δt, ẋ, tracker.state.jacobian)
        # check if we need to update the jacobian_infos
        update_jacobian_infos =
                state.last_step_failed ||
                state.steps_jacobian_info_update ≥ options.steps_jacobian_info_update
        # reset counter
        update_jacobian_infos && (state.steps_jacobian_info_update = 0)
        result = correct!(x̄, tracker, x̂, t + Δt; update_jacobian_infos=update_jacobian_infos)
        if isconverged(result)
            # Step is accepted, assign values
            state.accepted_steps += 1
            x .= x̄
            state.s += state.Δs
            state.accuracy = result.accuracy
            # Step size change
            state.Δs_prev = state.Δs
            update_stepsize!(tracker, result)
            if state.patch !== nothing && options.update_patch
                changepatch!(state.patch, x)
            end
            options.auto_scaling && auto_scaling!(state, options)
            # update derivative
            compute_ẋ!(state, cache, options)
            # tell the predictors about the new derivative if they need to update something
            update!(cache.predictor, H, x, ẋ, t + Δt, state.jacobian)
            state.steps_jacobian_info_update += 1

            # Check termination criterion: we became too ill-conditioned
            if options.terminate_ill_conditioned && is_ill_conditioned(state, options) && state.last_step_failed
                state.status = CoreTrackerStatus.terminated_ill_conditioned
                return nothing
            end

            state.last_step_failed = false
        else
            # We have to reset the patch
            state.rejected_steps += 1
            update_rank!(state.jacobian)
            # Step failed, so we have to try with a new (smaller) step size
            update_stepsize!(tracker, result)
            Δt = current_Δt(state)
            state.last_step_failed = true
        end
    catch err
        if !(err isa LinearAlgebra.SingularException)
            rethrow(err)
        end
        tracker.state.status = CoreTrackerStatus.terminated_singularity
    end
    nothing
end

function is_ill_conditioned(state::CoreTrackerState, options::CoreTrackerOptions)
    d = max(log₁₀(cond(state)), digits_lost(state))
    d > options.max_lost_digits || state.jacobian.corank > 0
end

function is_overdetermined_tracking(tracker)
    m, n = size(tracker.cache.homotopy)
    m > n
end

function check_min_step_size!(tracker::CoreTracker{<:LogHomotopy})
    @unpack state, options = tracker
    if abs((exp(-state.Δs) - 1) * exp(-state.s)) < options.min_step_size
        state.status = CoreTrackerStatus.terminated_step_size_too_small
    end
end
function check_min_step_size!(tracker::CoreTracker)
    @unpack state, options = tracker
    if state.Δs < options.min_step_size
        state.status = CoreTrackerStatus.terminated_step_size_too_small
    end
end

g(Θ) = sqrt(1+4Θ) - 1
# Choose 0.25 instead of 1.0 due to Newton-Kantorovich theorem
δ(opts::CoreTrackerOptions, ω, μ) = @fastmath min(√(0.5ω) * τ(opts, μ), 0.25)
τ(opts::CoreTrackerOptions, μ) = @fastmath opts.accuracy^(1/(2opts.max_corrector_iters - μ))

function update_stepsize!(tracker::CoreTracker, result::CorrectorResult)
    @unpack state, options = tracker
    model = state.step_size
    ord = order(tracker.predictor)

    if options.simple_step_size_alg
        simple_step_size_alg!(tracker, result)
        return nothing
    end

    # The step size control is described in https://arxiv.org/abs/1902.02968

    # we have to handle the special case that there is only 1 iteration
    # in this case we cannot estimate ω and therefore just assume ω = 2
    # Also note ||x̂-x|| = ||Δx₀||
    if result.iters == 1
        ω = isnan(model.ω) ? 2.0 : 0.5model.ω
        d_x̂_x̄ = result.norm_Δx₀
    else
        ω = result.ω₀
        d_x̂_x̄ = distance(state.x̂, state.x, state.inner_product)
    end
    Δx₀ = result.norm_Δx₀
    if isconverged(result)
        δ_N_ω = δ(options, ω, 0.25)
        # This is to handle the edge case that g(δ_N_ω) >> (ω * d_x̂_x̄) > 0 but
        # at the same time δ_N_ω < eps(). Since then g(δ_N_ω) = 0
        δ_N_ω = max(δ_N_ω, 1e-15) #
        λ = g(δ_N_ω) / (ω * d_x̂_x̄)
        Δs = nthroot(λ, ord) * state.Δs
        if state.last_step_failed
            Δs = min(Δs, state.Δs)
        end
    else
        if isfinite(model.ω)
            ω = max(ω, model.ω)
        end
        δ_N_ω = δ(options, ω, 0.25)
        ω_η = 0.5ω * Δx₀
        if δ_N_ω < ω_η
            λ = g(δ_N_ω) / g(ω_η)
            Δs = min(nthroot(λ, ord), 0.9) * state.Δs
        else
            Δs = 0.5 * state.Δs
        end
        # If we fail consecutively reduce step size more aggressively.
        if state.last_step_failed
            Δs = min(Δs, 0.25 * state.Δs)
        end
    end
    model.expected_Δx₀ = 2δ_N_ω / ω
    model.ω = ω
    Δs = min(Δs, options.max_step_size)
    if (length(state.segment) - state.s) < Δs
        state.Δs = length(state.segment) - state.s
    else
        state.Δs = Δs
        check_min_step_size!(tracker)
    end

    nothing
end

function simple_step_size_alg!(tracker::CoreTracker, result::CorrectorResult)
    @unpack state, options = tracker
    if isconverged(result)
        state.consecutive_successfull_steps += 1
        if state.consecutive_successfull_steps == 5
            Δs′ = 2 * state.Δs
            state.consecutive_successfull_steps = 0
        else
            Δs′ = state.Δs
        end
    else
        state.consecutive_successfull_steps = 0
        Δs′ = 0.5 * state.Δs
    end

    state.Δs = min(Δs′, length(state.segment) - state.s)

    if !isconverged(result)
        check_min_step_size!(tracker)
    end
end

function auto_scaling!(state::CoreTrackerState, options::CoreTrackerOptions)
    auto_scaling!(state.inner_product, state.x, options.auto_scaling_options)
end
function auto_scaling!(ip::WeightedIP, x::AbstractVector, opts::AutoScalingOptions)
    norm_x = ip(x)
    for i in 1:length(x)
        wᵢ = (abs(x[i]) + ip.weight[i]) / 2
        if wᵢ < opts.scale_min * norm_x
            wᵢ = opts.scale_min * norm_x
        elseif wᵢ > opts.scale_max * norm_x
            wᵢ = opts.scale_max * norm_x
        end
        ip.weight[i] = max(wᵢ, opts.scale_abs_min)
    end
    nothing
end
auto_scaling!(ip::EuclideanIP, x::AbstractVector, opts::AutoScalingOptions) = nothing

function check_terminated!(tracker::CoreTracker)
    if abs(tracker.state.s - length(tracker.state.segment)) < 2eps(length(tracker.state.segment))
        tracker.state.status = CoreTrackerStatus.success
    elseif iters(tracker) ≥ tracker.options.max_steps
        tracker.state.status = CoreTrackerStatus.terminated_maximal_iterations
    end
    nothing
end

function refine!(tracker::CoreTracker, accuracy=tracker.options.refinement_accuracy)
    if tracker.state.accuracy < accuracy
        return
    end
    result = correct!(tracker.state.x̄, tracker;
        accuracy=accuracy,
        max_iters=tracker.options.max_refinement_iters)
    if isconverged(result)
        tracker.state.x .= tracker.state.x̄
        tracker.state.accuracy = result.accuracy
    end
    nothing
end

function Base.iterate(tracker::CoreTracker, state=nothing)
    state === nothing && return tracker, 1

    if tracker.state.status == CoreTrackerStatus.tracking
        step!(tracker)
        check_terminated!(tracker)

        if tracker.state.status == CoreTrackerStatus.success
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
affine_tracking(tracker::CoreTracker) = !isa(tracker.state.x, ProjectiveVectors.PVector)

#################
## Query State ##
#################
"""
     current_t(tracker::CoreTracker)

Current `t`.
"""
current_t(tracker::CoreTracker) = current_t(tracker.state)
current_t(state::CoreTrackerState) = state.segment[state.s]

"""
     current_Δt(tracker::CoreTracker)

Current step_size `Δt`.
"""
current_Δt(tracker::CoreTracker) = current_Δt(tracker.state)
current_Δt(state::CoreTrackerState) = state.segment[state.Δs] - state.segment.start

"""
     iters(tracker::CoreTracker)

Current number of iterations.
"""
iters(tracker::CoreTracker) = iters(tracker.state)
iters(state::CoreTrackerState) = state.accepted_steps + state.rejected_steps

"""
     status(tracker::CoreTracker)

Current status.
"""
status(tracker::CoreTracker) = status(tracker.state)
status(state::CoreTrackerState) = state.status

"""
    current_x(tracker::CoreTracker)

Return the current value of `x`.
"""
current_x(tracker::CoreTracker) = current_x(tracker.state)
current_x(state::CoreTrackerState) = state.x


"""
    LinearAlgebra.cond(tracker::CoreTracker)

Returns the currently computed approximation of the condition number of the
Jacobian.
"""
cond(tracker::CoreTracker) = LinearAlgebra.cond(tracker.state)
cond(state::CoreTrackerState) = state.jacobian.cond


"""
    digits_lost(tracker::CoreTracker)

Returns the currently computed approximation of the number of digits lost
during the linear system solving in Newton's method.
"""
digits_lost(tracker::CoreTracker) = digits_lost(tracker.state)
digits_lost(state::CoreTrackerState) = unpack(state.jacobian.digits_lost, 0.0)


"""
    inner(tracker::CoreTracker)

Returns the inner product used to compute distance during the path tracking.
"""
inner(tracker::CoreTracker) = inner(tracker.state)
inner(state::CoreTrackerState) = state.inner_product

"""
    options(tracker::CoreTracker)

Returns the options used in the tracker.
"""
options(tracker::CoreTracker) = tracker.options

##################
# Modify options #
##################
"""
     accuracy(tracker::CoreTracker)

Current accuracy.
"""
accuracy(tracker::CoreTracker) = tracker.options.accuracy

"""
     set_accuracy!(tracker::CoreTracker, accuracy; update_max_lost_digits=true)

Set the current accuracy to `accuracy`. If `update_max_lost_digits` is `true` then
the setting `max_lost_digits` will be updated to the default setting.
"""
function set_accuracy!(tracker::CoreTracker, accuracy; update_max_lost_digits::Bool=true)
    @unpack options = tracker
    options.accuracy = accuracy
    if update_max_lost_digits
        options.max_lost_digits = default_max_lost_digits(options.precision, accuracy)
    end
    tracker
end

"""
     refinement_accuracy(tracker::CoreTracker)

Current refinement accuracy.
"""
refinement_accuracy(tracker::CoreTracker) = tracker.options.refinement_accuracy

"""
     set_max_refinement_iters!(tracker::CoreTracker, accuracy)

Set the current refinement accuracy to `accuracy`.
"""
set_refinement_accuracy!(T::CoreTracker, accuracy) = T.options.refinement_accuracy = accuracy

"""
     max_refinement_iters(tracker::CoreTracker)

Current refinement max_steps.
"""
max_refinement_iters(T::CoreTracker) = T.options.max_refinement_iters

"""
     set_max_refinement_iters!(tracker::CoreTracker, n)

Set the current refinement max_steps to `n`.
"""
set_max_refinement_iters!(T::CoreTracker, n) = T.options.max_refinement_iters = n

"""
     max_corrector_iters(tracker::CoreTracker)

Current correction max_steps.
"""
max_corrector_iters(T::CoreTracker) = T.options.max_corrector_iters

"""
     set_max_corrector_iters!(tracker::CoreTracker, n)

Set the correction max_steps to `n`.
"""
set_max_corrector_iters!(T::CoreTracker, n) = T.options.max_corrector_iters = n

"""
     max_step_size (tracker::CoreTracker)

Current maximal step size.
"""
max_step_size(T::CoreTracker) = T.options.max_step_size

"""
     set_max_corrector_iters!(tracker::CoreTracker, Δs)

Set the maximal step size to `Δs`.
"""
set_max_step_size!(T::CoreTracker, Δs) = T.options.max_step_size  = Δs

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
    iterator(tracker::CoreTracker, x₁, t₁=1.0, t₀=0.0; affine=true)

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
println("Success: ", status(tracker) == CoreTrackerStatus.success)
```
"""
function iterator(tracker::CoreTracker, x₁, t₁=1.0, t₀=0.0; kwargs...)
    setup!(tracker, x₁, t₁, t₀; kwargs...)
    PathIterator(tracker, typeof(t₁ - t₀) <: Real)
end

function current_x_t(iter::PathIterator)
    x = current_x(iter.tracker)
    t = current_t(iter.tracker)
    (x, iter.t_real ? real(t) : t)
end

function Base.iterate(iter::PathIterator, state=nothing)
    state === nothing && return current_x_t(iter), 1
    iter.tracker.state.status != CoreTrackerStatus.tracking && return nothing

    step_done = false
    while !step_done && (iter.tracker.state.status == CoreTrackerStatus.tracking)

        step!(iter.tracker)
        check_terminated!(iter.tracker)

        if iter.tracker.state.status == CoreTrackerStatus.success
            refine!(iter.tracker)
        end

        step_done = !iter.tracker.state.last_step_failed
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
    supported, rest = splitkwargs(kwargs,problem_startsolutions_supported_keywords)
    prob, startsolutions = problem_startsolutions(args...; supported...)
    tracker = CoreTracker(prob, start_solution_sample(startsolutions); rest...)
    (tracker=tracker, startsolutions=startsolutions)
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
tracker = coretracker(f, parameters=p, p₁=a, p₀=b)
```
You then can obtain a single solution at `b` by using
```julia
x_b = track(tracker, x_a).x
```

### Trace a path
To trace a path you can use the [`iterator`](@ref) method.

```julia
tracker = coretracker(f, parameters=p, p₁=a, p₀=b, max_step_size =0.01)
for (x, t) in iterator(tracker, x₁)
    @show (x,t)
end
```

If we want to guarantee smooth traces we can limit the maximal step size.
```julia
tracker = coretracker(f, parameters=p, p₁=a, p₀=b, max_step_size =0.01)
for (x, t) in iterator(tracker, x₁)
    @show (x,t)
end
```
"""
function coretracker(args...; kwargs...)
    tracker, _ = coretracker_startsolutions(args...; kwargs...)
    tracker
end
