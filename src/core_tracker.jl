export CoreTrackerStatus,
       is_success,
       is_terminated,
       is_tracking,
       is_invalid_startvalue,
       CoreTrackerResult,
       solution,
       CoreTrackerOptions,
       CoreTracker,
       coretracker,
       coretracker_startsolutions,
       affine_tracking,
       track,
       track!,
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
       set_accuracy!,
       set_max_corrector_iters!,
       set_max_step_size!,
       options,
       # deprecated
       digits_lost,
       setup!,
       refinement_accuracy,
       max_refinement_iters,
       set_max_refinement_iters!,
       set_refinement_accuracy!,
       start_parameters!,
       target_parameters!



const coretracker_supported_keywords = [
    :accuracy,
    :auto_scaling,
    :initial_step_size,
    :min_step_size,
    :max_cond,
    :max_corrector_iters,
    :max_step_size,
    :max_steps,
    :norm,
    :patch,
    :precision,
    :predictor,
    :simple_step_size_alg,
    :terminate_ill_conditioned,
    :log_homotopy,
    :log_transform,
    :logarithmic_time_scale,
    :from_infinity,
    # deprecated,
    :max_lost_digits,
    :max_refinement_iters,
    :refinement_accuracy,
]


####################
# CoreTrackerState #
####################

module CoreTrackerStatus

@doc """
    CoreTrackerStatus.states

The possible states a `CoreTracker` can have are

* `CoreTrackerStatus.success`: Indicates a successfull completed tracking
* `CoreTrackerStatus.tracking`: The tracking is still in progress
* `CoreTrackerStatus.terminated_maximal_iterations`: Tracking terminated since maximal iterations reached.
* `CoreTrackerStatus.terminated_ill_conditioned`: Tracking terminated since the path was too ill-conditioned.
* `CoreTrackerStatus.terminated_invalid_startvalue`: Tracking terminated since the provided start value was invalid.
* `CoreTrackerStatus.terminated_step_size_too_small`
"""
@enum states begin
    success
    tracking
    terminated_maximal_iterations
    terminated_accuracy_limit
    terminated_ill_conditioned
    terminated_invalid_startvalue
    terminated_step_size_too_small
        # no more used
    terminated_singularity
end
end

"""
    is_success(S::CoreTrackerStatus.states)

Returns `true` if `S` indicates a success in the path tracking.
"""
is_success(S) = S == CoreTrackerStatus.success

"""
    is_terminated(S::CoreTrackerStatus.states)

Returns `true` if `S` indicates that the path tracking got terminated. This is not `true`
if `is_success(S)` is `true`.
"""
is_terminated(S::CoreTrackerStatus.states) =
    S ≠ CoreTrackerStatus.tracking && S ≠ CoreTrackerStatus.success

"""
    is_invalid_startvalue(S::CoreTrackerStatus.states)

Returns `true` if `S` indicates that the path tracking got terminated since the start
value was not a zero.
"""
is_invalid_startvalue(S::CoreTrackerStatus.states) =
    S == CoreTrackerStatus.terminated_invalid_startvalue

"""
    is_tracking(S::CoreTrackerStatus.states)

Returns `true` if `S` indicates that the path tracking is not yet finished.
"""
is_tracking(S::CoreTrackerStatus.states) = S == CoreTrackerStatus.tracking

####################
# CoreTrackerResult #
####################

"""
     CoreTrackerResult{V<:AbstractVector}

The result of [`track(::CoreTracker, x₁, t₁, t₀)](@ref). You can use
[`is_success`](@ref) to check whether the tracking was successfull and [`solution`](@ref)
to obtain the solution.

# Fields

* `returncode` The [`CoreTrackerStatus.states`](@ref) enum.
* `x::V` The solution at `t`.
* `t::ComplexF64` The time `t` when the path tracker stopped.
* `accuracy::Float64`: The estimated accuracy of `x`.
* `accepted_steps::Int`: The number of accepted steps during the tracking.
* `rejected_steps::Int`: The number of rejected_steps steps during the tracking.
"""
struct CoreTrackerResult{V<:AbstractVector}
    returncode::CoreTrackerStatus.states
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
        state.rejected_steps,
    )
end

"""
    is_success(result::CoreTrackerResult)

Returns `true` if the path tracking was successfull.
"""
is_success(result::CoreTrackerResult) = is_success(result.returncode)

"""
    solution(result::CoreTrackerResult)

Returns the solutions obtained by the `CoreTracker`.
"""
solution(result::CoreTrackerResult) = result.x


Base.show(io::IO, result::CoreTrackerResult) = print_fieldnames(io, result)
Base.show(io::IO, ::MIME"application/prs.juno.inline", result::CoreTrackerResult) = result


######################
# CoreTrackerOptions #
######################
"""
    PrecisionOption

Controls the used precision for computing the residual in Newton's method.
See [[Tisseur01]](https://epubs.siam.org/doi/abs/10.1137/S0895479899359837) for the analysis behind this approach.

## Values
* `PRECISION_FIXED_64`: Only use default 64 bit machine precision
* `PRECISION_FIXED_128`: Always use emulated 128 bit floating point numbers. These are provided by the [DoubleFloats.jl](https://github.com/JuliaMath/DoubleFloats.jl) package.
* `PRECISION_ADAPTIVE`: Adaptively switch between 64 and 128 bit floating point numbers.
"""
@enum PrecisionOption begin
    PRECISION_FIXED_64
    PRECISION_FIXED_128
    PRECISION_ADAPTIVE
end

"""
    CoreTrackerOptions

The set of options set for a [`CoreTracker`](@ref). See the description of [`CoreTracker`](@ref)
for all possible options.
"""
mutable struct CoreTrackerOptions
    accuracy::Float64
    initial_step_size::Union{Nothing,Float64}
    max_cond::Float64
    max_corrector_iters::Int
    max_step_size::Float64
    max_steps::Int
    min_step_size::Float64
    precision::PrecisionOption
    simple_step_size_alg::Bool
    steps_jacobian_info_update::Int
    terminate_ill_conditioned::Bool
    track_cond::Bool
    update_patch::Bool
    # Not changeable options
    logarithmic_time_scale::Bool
    from_infinity::Bool
end

function CoreTrackerOptions(
    ;
    # internal
    parameter_homotopy = false,
    update_patch = true,
    logarithmic_time_scale = false,
    from_infinity = false,
    # public
    accuracy = 1e-7,
    initial_step_size = nothing,
    max_cond::Float64 = 1e12,
    max_corrector_iters::Int = 2,
    max_step_size = Inf,
    max_steps = parameter_homotopy ? 10_000 : 1_000,
    min_step_size = 1e-14,
    precision::Union{PrecisionOption,Symbol} = :double,
    simple_step_size_alg = false,
    steps_jacobian_info_update::Int = 5,
    terminate_ill_conditioned::Bool = true,
    track_cond::Bool = true,
     # deprecated
    max_lost_digits = nothing,
    max_refinement_iters = nothing,
    refinement_accuracy = nothing,
)

    max_lost_digits !== nothing && @warn(
        "Passing `max_lost_digits` to `CoreTracker` is deprecated.",
    )
    refinement_accuracy !== nothing && @warn(
        "Passing `refinement_accuracy` to `CoreTracker` is deprecated. " *
        "Solutions are now automatically refined to the maximal achievable accuracy.",
    )
    max_refinement_iters !== nothing && @warn(
        "Passing `max_refinement_iters` to `CoreTracker` is deprecated. " *
        "Solutions are now automatically refined to the maximal achievable accuracy.",
    )

    CoreTrackerOptions(
        accuracy,
        initial_step_size,
        max_cond,
        max_corrector_iters,
        max_step_size,
        max_steps,
        min_step_size,
        make_precision(precision),
        simple_step_size_alg,
        steps_jacobian_info_update,
        terminate_ill_conditioned,
        track_cond,
        update_patch,
        logarithmic_time_scale,
        from_infinity,
    )
end

make_precision(p::PrecisionOption) = p
function make_precision(p::Symbol)
    if p == :double
        PRECISION_FIXED_64
    elseif p == :double_double
        PRECISION_FIXED_128
    elseif p == :adaptive
        PRECISION_ADAPTIVE
    else
        throw(ArgumentError("Unsupported argument `precision=$p`. " *
                            "Possible values are `:double`, `:double_double` and `:adaptive`."))
    end
end

Base.show(io::IO, opts::CoreTrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::CoreTrackerOptions) = opts

Base.copy(opts::CoreTrackerOptions) = copy!(CoreTrackerOptions(), opts)
function Base.copy!(opts::CoreTrackerOptions, opts2::CoreTrackerOptions)
    @unpack accuracy,
        initial_step_size,
        min_step_size,
        max_cond,
        max_corrector_iters,
        max_step_size,
        max_steps,
        precision,
        simple_step_size_alg,
        steps_jacobian_info_update,
        terminate_ill_conditioned,
        track_cond,
        update_patch,
        logarithmic_time_scale,
        from_infinity = opts2
    @pack! opts = accuracy,
        initial_step_size,
        min_step_size,
        max_cond,
        max_corrector_iters,
        max_step_size,
        max_steps,
        precision,
        simple_step_size_alg,
        steps_jacobian_info_update,
        terminate_ill_conditioned,
        track_cond,
        update_patch,
        logarithmic_time_scale,
        from_infinity
    opts
end

####################
# CoreTrackerState #
####################

mutable struct CoreTrackerState{
    T,
    AV<:AbstractVector{Complex{T}},
    AN<:AbstractNorm,
    MaybePatchState<:Union{AbstractAffinePatchState,Nothing},
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
    norm_Δx₀::Float64 # norm of the first newton update
    ω::Float64 # liptschitz constant estimate, see arxiv:1902.02968
    limit_accuracy::Float64
    eval_err::Float64

    norm::AN
    jacobian::JacobianMonitor{T}
    residual::Vector{Float64}

    steps_jacobian_info_update::Int
    status::CoreTrackerStatus.states
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
    norm::AbstractNorm,
)
    x = isa(x₁, SVector) ? Vector(x₁) : copy(x₁)
    x̂, x̄ = copy(x), copy(x)
    ẋ = Vector(x)
    segment = ComplexSegment(promote(complex(t₁), complex(t₀))...)
    s = 0.0
    Δs = convert(
        Float64,
        min(unpack(options.initial_step_size, Inf), length(segment), options.max_step_size),
    )
    Δs_prev = 0.0
    norm_Δx₀ = accuracy = 0.0
    limit_accuracy = eval_err = eps()
    ω = 1.0
    JM = JacobianMonitor(jacobian(H, x, t₁))
    residual = zeros(first(size(H)))
    steps_jacobian_info_update = 0
    accepted_steps = rejected_steps = 0
    status = CoreTrackerStatus.tracking
    last_step_failed = true
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
        norm_Δx₀,
        ω,
        limit_accuracy,
        eval_err,
        norm,
        JM,
        residual,
        steps_jacobian_info_update,
        status,
        patch,
        accepted_steps,
        rejected_steps,
        last_step_failed,
        consecutive_successfull_steps,
    )
end

Base.show(io::IO, state::CoreTrackerState) = print_fieldnames(io, state)
Base.show(io::IO, ::MIME"application/prs.juno.inline", state::CoreTrackerState) = state


##################
## CoreTracker ##
##################
"""
    CoreTracker

A `CoreTracker` is the low-level path tracker. Its job is to track an initial solution
`x₁` from `t₁` to `t₀` by following a solution path ``x(t)``. For this a *predictor-corrector* scheme is used.
See our [introduction guide](https://www.juliahomotopycontinuation.org/guides/introduction/#tracking-solution-paths)
for a high-level explanation of the predictor corrector scheme.
The `CoreTracker` accepts for a value `t` a solution ``x^*`` if ``||x(t) - x^*|| < τ`` where ``τ``
is controlled by the `accuracy` option.

To interact with a `CoreTracker` take a look at [`track`](@ref) resp. [`track!`](@ref).
You can also use `CoreTracker` as an iterator. For this you need to call
[`init!(::CoreTracker, x₁, t₁, t₀)`](@ref) first. There is also [`iterator`](@ref) if you
are only interested in pairs `(x,t)`.

The `CoreTracker` **cannot** handle singular solutions or divergent paths. For this you have
to use a [`PathTracker`](@ref) (which uses internally again a `CoreTracker`).

The `CoreTracker` can track solutions in projective space if the underlying homotopy is homogeneous
and the provided start solution is a projective vector (with type [`PVector`](https://www.juliahomotopycontinuation.org/ProjectiveVectors.jl/stable/#ProjectiveVectors.PVector))
from the [`ProjectiveVectors.jl`](https://github.com/JuliaHomotopyContinuation/ProjectiveVectors.jl) package.

The `CoreTracker` can also perform a reparamterization of the parameter ``t`` from ``[0,1]`` towards
``[0,∞]`` by ``t = e^{-s}``. If the homotopy assumes ``t ∈ [0,1]`` this can be done by
passing the `log_transform = true` option. If you already have a homotopy which assumes
``s ∈ [0,∞]`` you can use the `logarithmic_time_scale = true` to tell the core tracker this.
In this case the `min_step_size` criterion is still applied to ``t`` and not ``s``.


## Options

`CoreTracker` accepts many different options. Here are the options in alphabetic order:

* `accuracy` (default `1e-7`): The maximal acceptable distance to the true solution at a
  time step `t`.
* `auto_scaling` (default `true`): Enables an automatic scaling of the variables by
  introducing a *weighted norm*. See also [`WeightedNorm`](@ref).
* `initial_step_size`: The size of the first step. By default an automatic estimate is
  computed but this can be overwritten with this option.
* `max_cond` (default `1e10`): The maximal condition number before a path is terminated due
  to ill-conditioning. This is only applicable if `terminate_ill_conditioned` is `true`.
* `max_corrector_iters (default `2`): The maximal number of Newton iteration steps until
  which the desired accuracy has to be achieved.
* `max_step_size` (default `Inf`): Limit the path tracker to a maximal step size.
* `max_steps` (default `1000`): The maximal number of steps the path tracker makes.
  Note that this changes to `10_000` for parameter homotopies.
* `min_step_size` (default `1e-14`): The path tracking is terminated if the step size is
  smaller than this value.
* `norm` (default [`InfNorm()`](@ref)): The norm used in the computations. Note that if the
  `auto_scaling` is true, this will automatically be converted to a [`WeightedNorm`](@ref).
* `patch`: It is possible to perform tracking in projective space. For the actual
  compuations we still need to choose an affine chart. There are multiple charts possible,
  see [`AbstractAffinePatch`](@ref). The default is [`RandomPatch`](@ref).
* `precision` (default `:double`): This controls the precision used in the computation
  of ``H(x,t)``. If it is set to `:double` only double precision arithmetic is used
  (`Complex{Float64}`).  If it set to `:double_double` double double arithmetic
  as implemented by the [`DoubleFloats`](https://github.com/JuliaMath/DoubleFloats.jl)
  package is **always** used. The option `:adaptive` changes adaptively the precision.
* `predictor` (default [`Pade21`](@ref)): The predictor used.
* `simple_step_size_alg` (default `false`): By default the step size algorithm presented in
  [arXiv:1902.02968](https://arxiv.org/abs/1902.02968) is used. If this option is true
  a more simple step size algorithm is used which doubles the step size after 5 consecutive
  successive steps and halfes the step size after a failed step.
* `terminate_ill_conditioned` (default `true`): Indicates whether the path tracking should be
  terminated for ill-conditioned paths. A path is considered ill-conditioned if the desirable
  accuracy is no more achievable or if the condition number of the Jacobian is larger than
  `1e14`.

"""
struct CoreTracker{
    T,
    AV<:AbstractVector{Complex{T}},
    AN<:AbstractNorm,
    P<:AbstractPredictorCache,
    AP<:Union{Nothing,AbstractAffinePatchState},
    H<:HomotopyWithCache,
    NCC<:NewtonCorrectorCache,
}
    homotopy::H
    predictor::P
    corrector::NCC
    # these are mutable
    state::CoreTrackerState{T,AV,AN,AP}
    options::CoreTrackerOptions
    # buffer
    dt::Vector{Complex{T}}
end

function CoreTracker(
    prob::AbstractProblem,
    x₁;
    log_transform::Bool = false,
    log_homotopy::Bool = log_transform,
    logarithmic_time_scale = log_transform || log_homotopy,
    kwargs...,
)
    if log_homotopy
        H = LogHomotopy(prob.homotopy)
        CoreTracker(
            H,
            embed(prob, x₁),
            complex(0.0),
            complex(36.0),
            prob;
            logarithmic_time_scale = logarithmic_time_scale,
            kwargs...,
        )
    else
        CoreTracker(
            prob.homotopy,
            embed(prob, x₁),
            complex(1.0),
            complex(0.0),
            prob;
            logarithmic_time_scale = logarithmic_time_scale,
            kwargs...,
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
    predictor::AbstractPredictor = Pade21(),
    auto_scaling::Bool = true,
    norm::AbstractNorm = InfNorm(),
    kwargs...,
)

    patch === nothing || throw(ArgumentError("You can only pass `patch=$(patch)` if `affine_tracking=false`."))

    options = CoreTrackerOptions(
        ;
        parameter_homotopy = isa(homotopy, ParameterHomotopy),
        kwargs...,
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
        used_norm,
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
    patch = RandomPatch(),
    corrector::NewtonCorrector = NewtonCorrector(),
    predictor::AbstractPredictor = Pade21(),
    auto_scaling::Bool = true,
    norm::AbstractNorm = InfNorm(),
    simple_step_size_alg = isa(patch, OrthogonalPatch),
    kwargs...,
)

    options = CoreTrackerOptions(
        ;
        parameter_homotopy = isa(homotopy, ParameterHomotopy),
        kwargs...,
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
        used_norm,
    )

    pred_cache = cache(predictor, H, ct_state.x, ct_state.ẋ, t₁)
    corr_cache = cache(corrector, H, ct_state.x, t₁)
    dt = Vector{eltype(ct_state.ẋ)}(undef, first(size(H)))

    CoreTracker(H, pred_cache, corr_cache, ct_state, options, dt)
end

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

Base.broadcastable(C::CoreTracker) = Ref(C)

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
    t::Number = current_t(tracker.state);
    simplified_last_step::Bool = true,
    min_iters::Int = 1,
    max_iters::Int = tracker.options.max_corrector_iters + 1,
    full_steps::Int = max_iters - simplified_last_step,
    tol::Float64 = tracker.options.accuracy,

)
    if tracker.options.precision == PRECISION_ADAPTIVE
        double_64_evaluation = at_limit_accuracy(
            tracker.state,
            tracker.options;
            safety_factor = 100.0,
        )
    else
        double_64_evaluation = tracker.options.precision == PRECISION_FIXED_128
    end

    newton!(
        x̄,
        tracker.homotopy,
        x,
        t,
        tracker.state.jacobian,
        tracker.state.norm,
        tracker.corrector;
        tol = tol,
        max_iters = max_iters,
        min_iters = min_iters,
        double_64_evaluation = double_64_evaluation,
        full_steps = full_steps,
    )
end

function limit_accuracy!(
    tracker::CT,
    x::AbstractVector = tracker.state.x,
    t::Number = current_t(tracker.state),
    norm::AbstractNorm = tracker.state.norm;
    update_cond::Bool = tracker.options.track_cond,
    update_eval_err::Bool = true,
)
    @unpack state, options = tracker

    for i = 1:2
        # We do a second iteration if we would terminate the path tracking to make sure
        # that this is not just an artifact
        i == 2 && !terminate_limit_accuracy(state, options) && continue
        limit_acc, eval_err = limit_accuracy!(
            tracker.homotopy,
            x,
            t,
            state.jacobian,
            norm,
            tracker.corrector;
            x_accuracy = state.accuracy,
            update_jacobian = false,
            compute_cond = update_cond,
            compute_eval_err = update_eval_err,
        )
        state.accuracy = min(state.accuracy, limit_acc)
        state.limit_accuracy = max(state.accuracy, limit_acc)
        if update_eval_err
            state.eval_err = exp((log(state.eval_err) + log(eval_err)) / 2)
        end
    end

    nothing
end

function terminate_limit_accuracy(state, options)
    (options.precision == PRECISION_FIXED_64) && at_limit_accuracy(
        state,
        options;
        safety_factor = 10.0,
    )
end
function at_limit_accuracy(state::CTS, opts::CTO; safety_factor::Float64 = 10.0)
    opts.accuracy < safety_factor * state.limit_accuracy
end
####################
## INITIALIZATION ##
####################

embed!(x::ProjectiveVectors.PVector, y) = ProjectiveVectors.embed!(x, y)
embed!(x::AbstractVector, y) = x .= y

function check_start_value!(tracker::CT, x = tracker.state.x, t = current_t(tracker.state))
    # We perturb the provided soution by sligthly more than the desired accuracy
    # to check that the Newton iteration converges properly and to obtain an estimate
    # of ω
    perturb!(tracker.state.x̄, x, 10 * tracker.options.accuracy, tracker.state.norm)

    # Solutions could need higher order precision already at the beginning
    # Therefore, if the standard check fails, try again with higher precision (if applicable)
    double_64_evaluation = false
    for i = 1:2
        result = newton!(
            tracker.state.x̄,
            tracker.homotopy,
            tracker.state.x̄,
            t,
            tracker.state.jacobian,
            tracker.state.norm,
            tracker.corrector;
            tol = tracker.options.accuracy,
            max_iters = 2,
            simplified_last_step = false,
            double_64_evaluation = double_64_evaluation,
        )
        if is_converged(result)
            tracker.state.ω = result.ω₀
            x .= tracker.state.x̄
            tracker.state.accuracy = result.accuracy
            limit_accuracy!(tracker)
            break
        elseif i == 1 && tracker.options.precision != PRECISION_FIXED_64
            double_64_evaluation = true
        else
            tracker.state.status = CoreTrackerStatus.terminated_invalid_startvalue
            break
        end
    end
    nothing
end

perturb!(y, x, Δ::Real, w::WeightedNorm) = y .= x .+ Δ .* weights(w)
perturb!(y, x, Δ::Real, norm::AbstractNorm) = y .= x .+ Δ

"""
    init!(tracker::CT, x₁, t₁, t₀;
      setup_patch::Bool = tracker.options.update_patch,
      loop::Bool = false,
      check_start_value::Bool = !loop)

Setup `tracker` to track `x₁` from `t₁` to `t₀`. Use this if you want to use the
`tracker` as an iterator.
"""
function init!(
    tracker::CoreTracker,
    x₁::AbstractVector,
    t₁,
    t₀;
    setup_patch::Bool = tracker.options.update_patch,
    loop::Bool = false,
    keep_steps::Bool = false,
    check_start_value::Bool = !loop,
)
    @unpack state, predictor, homotopy = tracker

    init!(state, x₁, t₁, t₀, tracker.options, setup_patch, loop, keep_steps)
    if !loop
        check_start_value!(tracker)
        compute_ẋ!(tracker; update_jacobian = false)
        init!(predictor, homotopy, state.x, state.ẋ, current_t(state), state.jacobian)
    end
    initial_step_size!(state, predictor, tracker.options)
    tracker
end

function init!(
    state::CTS,
    x₁::AbstractVector,
    t₁,
    t₀,
    options::CTO,
    setup_patch::Bool,
    loop::Bool,
    keep_steps::Bool,
)
    state.segment = ComplexSegment(t₁, t₀)
    state.s = 0.0
    state.Δs_prev = 0.0
    state.status = CoreTrackerStatus.tracking
    embed!(state.x, x₁)
    setup_patch && state.patch !== nothing && init!(state.patch, state.x)
    if !loop
        state.consecutive_successfull_steps = 0
        state.last_step_failed = true
        state.accuracy = 0.0
        state.ω = 1.0
        state.eval_err = state.limit_accuracy = eps()
        if !keep_steps
            state.accepted_steps = state.rejected_steps = 0
        end
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
    kwargs...,
)

##############
## TRACKING ##
##############

"""
    track(tracker::CoreTracker, x₁, t₁=1.0, t₀=0.0)::CoreTrackerResult

Track a value `x₁` from `t₁` to `t₀` using the given `CoreTracker` `tracker`.
This returns a `CoreTrackerResult`. This modifies `tracker`.
See [`track!`](@ref) for the possible options.
To investigate the behaviour of a particular take a look at [`path_info`](@ref).
"""
function track(
    tracker::CT,
    x₁::AbstractVector,
    t₁::Number = 1.0,
    t₀::Number = 0.0;
    kwargs...,
)
    track!(tracker, x₁, t₁, t₀; kwargs...)
    CoreTrackerResult(tracker)
end

"""
    track!(tracker::CoreTracker, x₁, t₁=1.0, t₀=0.0)::CoreTrackerStatus.states

Track a value `x₁` from `t₁` to `t₀` using the given `CoreTracker` `tracker`.
Returns one of the enum values of `CoreTrackerStatus` indicating the status.
Check [`is_success`](@ref) and [`is_terminated`](@ref) to test for the status.

If `setup_patch` is `true` then [`init!`](@ref) is called at the beginning
of the tracking.

    track!(x₀, tracker::CoreTracker, x₁, t₁=1.0, t₀=0.0)::CoreTrackerStatus.states

Additionally also stores the result in `x₀` if the tracking was successfull.
"""
function track!(
    tracker::CT,
    x₁::AbstractVector,
    t₁::Number = 1.0,
    t₀::Number = 0.0;
    setup_patch::Bool = tracker.options.update_patch,
    loop::Bool = false,
    check_start_value::Bool = !loop,
    debug::Bool = false,
)
    _track!(tracker, x₁, t₁, t₀, setup_patch, loop, check_start_value, debug)
end

function track!(
    x₀::AbstractVector,
    tracker::CT,
    x₁::AbstractVector,
    t₁::Number = 1.0,
    t₀::Number = 0.0;
    setup_patch::Bool = tracker.options.update_patch,
    loop::Bool = false,
    check_start_value::Bool = !loop,
    debug::Bool = false,
)

    _track!(tracker, x₁, t₁, t₀, setup_patch, loop, check_start_value, debug)
    retcode = status(tracker)
    if is_success(status(tracker))
        x₀ .= current_x(tracker)
    end
    retcode
end
@inline function _track!(
    tracker::CT,
    x₁::AbstractVector,
    t₁::Number,
    t₀::Number,
    setup_patch::Bool,
    loop::Bool,
    check_start_value::Bool,
    debug::Bool,
)
    t₁ == t₀ && return (tracker.state.status = CoreTrackerStatus.success)

    init!(
        tracker,
        x₁,
        t₁,
        t₀;
        setup_patch = setup_patch,
        check_start_value = check_start_value,
        loop = loop,
    )

    while is_tracking(tracker.state.status)
        step!(tracker, debug)
    end

    tracker.state.status
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
        color = :yellow,
    )

    # Use the current approximation of x(t) to obtain estimate
    # x̂ ≈ x(t + Δt) using the provided predictor
    predict!(x̂, predictor, homotopy, x, t, Δt, ẋ, state.jacobian)
    # Correct the predicted value x̂ to obtain x̄.
    # If the correction is successfull we have x̄ ≈ x(t+Δt).
    result = correct!(x̄, tracker, x̂, t + Δt; simplified_last_step = true)

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

        # Compute the new tangent ẋ(t).
        compute_ẋ!(tracker; update_jacobian = true)

        # Make an additional newton step to check the limiting accuracy
        if (state.steps_jacobian_info_update += 1) ≥ options.steps_jacobian_info_update ||
           state.last_step_failed
            limit_accuracy!(tracker)
            state.steps_jacobian_info_update = 0
        end
        # tell the predictors about the new derivative if they need to update something
        update!(predictor, homotopy, x, ẋ, t + Δt, state.jacobian, state.eval_err)

        # Update other state
        state.accepted_steps += 1
        state.last_step_failed = false
    else
        # Step failed, so we have to try with a new (smaller) step size
        update_stepsize!(tracker, result)

        # Update other state
        state.rejected_steps += 1
        state.last_step_failed = true
    end
    state.norm_Δx₀ = result.norm_Δx₀

    check_terminated!(state, options)

    # If we terminate, always update the limit accuracy (and by this refine the solution
    # to maximal accuracy).
    if is_success(state.status) && state.steps_jacobian_info_update ≠ 0
        limit_accuracy!(tracker)
    end

    state.last_step_failed
end

###############
## STEP SIZE ##
###############

## intial step size
function initial_step_size!(state::CTS, predictor::AbstractPredictorCache, options::CTO)
    if options.initial_step_size !== nothing
        return state.Δs = max(
            min(options.initial_step_size, length(state.segment), options.max_step_size),
            sqrt(options.min_step_size),
        )
    end
    ω = max(state.ω, 1e-8)
    p = order(predictor)
    deriv, k = unpack(highest_derivative(predictor), (state.ẋ, 1))
    η_p = nthroot(state.norm(deriv) / factorial(k), k)
    δ_N_ω = δ(options, ω, options.max_corrector_iters / 2)
    Δs = nthroot(g(δ_N_ω) / ω, p) / η_p
    if isnan(Δs)
        Δs = 0.05 * length(state.segment)
    end

    state.Δs = max(
        min(Δs, length(state.segment), options.max_step_size),
        options.min_step_size,
    )
end

## Step size update

function update_stepsize!(tracker::CT, result::NewtonCorrectorResult)
    @unpack state, options = tracker

    if !isnan(result.ω₀)
        # Use result.ω₀ in the converged case to avoid wrong values of ω due to
        # hitting the accuracy limit.
        if is_converged(result) || isnan(state.ω)
            state.ω = result.ω₀
        else
            state.ω = max(result.ω₀, state.ω)
        end
    end

    # Use the simple step size algorithm in the case that we get very close to
    # the accuracy limit. The motivation is that we can run into danger that the computed
    # estimate of ω is too high due to hitting limit accuracy.
    near_accuracy_limit = at_limit_accuracy(state, options; safety_factor = 100.0)
    if options.simple_step_size_alg ||
       (near_accuracy_limit && options.precision == PRECISION_FIXED_64)
        Δs = simple_step_size_alg!(state, options, result)
    else
        Δs = adaptive_step_size_alg!(state, options, result, order(tracker.predictor))
    end

    # Make sure to not overshoot.
    if !is_min_step_size(Δs, state.s, length(state.segment), options) && steps(state) > 5
        state.status = CoreTrackerStatus.terminated_step_size_too_small
    else
        state.Δs = min(length(state.segment) - state.s, Δs, options.max_step_size)
    end

    nothing
end

function is_min_step_size(Δs, s, length_segment, options::CTO)
    if options.logarithmic_time_scale && options.from_infinity
        ŝ = length_segment - s
        t = exp(-ŝ)
        Δt = exp(-ŝ + Δs) - t
        Δt ≥ min(options.min_step_size, t / 100)
    elseif options.logarithmic_time_scale
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
    state.consecutive_successfull_steps == 0 ? 2 * state.Δs : state.Δs
end

## Adaptive step size algorithm

"Use the adaptive step size control described in https://arxiv.org/abs/1902.02968"
function adaptive_step_size_alg!(
    state::CTS,
    options::CTO,
    result::NewtonCorrectorResult,
    ord::Int,
)
    @unpack s, ω = state

    d_x̂_x = result.iters > 1 ? state.norm(state.x̂, state.x) : result.norm_Δx₀
    Δx₀ = result.norm_Δx₀

    if is_converged(result)
        # This is to handle the edge case that g(δ_N_ω) >> (ω * d_x̂_x̄) > 0 but
        # at the same time δ_N_ω < eps(). Since then g(δ_N_ω) = 0
        δ_N_ω = max(δ(options, ω, options.max_corrector_iters / 4), 1e-15)
        Δs = nthroot(g(δ_N_ω) / (ω * d_x̂_x), ord) * state.Δs
        if state.last_step_failed
            Δs = min(Δs, state.Δs)
        end
    else
        δ_N_ω = δ(options, ω, options.max_corrector_iters / 4)
        ω_η = 0.5ω * Δx₀
        if δ_N_ω < ω_η
            Δs = min(nthroot(g(δ_N_ω) / g(ω_η), ord), 0.9) * state.Δs
        else
            Δs = 0.5 * state.Δs
        end
    end

    Δs
end

g(Θ::Real) = sqrt(1.0 + 4Θ) - 1.0
δ(opts::CTO, ω::Real, μ::Real) = min(√(0.5ω) * τ(opts, μ), 0.25)
function τ(opts::CTO, μ::Real)
    # most common case: 2(2 - 0.5) == 3
    if opts.max_corrector_iters == 2 && μ == 0.5
        cbrt(opts.accuracy)
    else
        opts.accuracy^(1 / (2 * (opts.max_corrector_iters - μ)))
    end
end


function check_terminated!(state::CTS, options::CTO)
    Δs = length(state.segment) - state.s
    if Δs < options.min_step_size
        state.status = CoreTrackerStatus.success
        state.s = length(state.segment)
    elseif steps(state) ≥ options.max_steps
        state.status = CoreTrackerStatus.terminated_maximal_iterations
    elseif options.terminate_ill_conditioned && terminate_limit_accuracy(state, options)
        state.status = CoreTrackerStatus.terminated_accuracy_limit
    elseif options.terminate_ill_conditioned && cond(state) > options.max_cond
        state.status = CoreTrackerStatus.terminated_ill_conditioned
    end
    nothing
end


function Base.iterate(tracker::CoreTracker, state::Int = 0)
    if is_tracking(tracker.state.status)
        step!(tracker)
        tracker, state + 1
    else
        nothing
    end
end

function is_valid_start_value(tracker::CoreTracker, x::AbstractVector, t::Number)
    @unpack state = tracker
    embed!(state.x̄, x)
    state.patch !== nothing && init!(state.patch, state.x̄)
    check_start_value!(tracker, state.x̄, t)
    !is_invalid_startvalue(status(tracker))
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
    current_t(tracker::CoreTracker)

Current `t`.
"""
current_t(tracker::CoreTracker) = current_t(tracker.state)
current_t(state::CTS) = state.segment[state.s]

"""
    current_Δt(tracker::CoreTracker)

Current step_size `Δt`.
"""
current_Δt(tracker::CoreTracker) = current_Δt(tracker.state)
current_Δt(state::CTS) = state.segment[state.Δs] - state.segment.start

"""
    steps(tracker::CoreTracker)

Current number of steps.
"""
steps(tracker::CoreTracker) = steps(tracker.state)
steps(state::CTS) = state.accepted_steps + state.rejected_steps

@deprecate iters(tracker::CoreTracker) steps(tracker)


"""
    status(tracker::CoreTracker)

Current status.
"""
status(tracker::CoreTracker) = status(tracker.state)
status(state::CTS) = state.status

"""
    current_x(tracker::CoreTracker)

Return the current value of `x`.
"""
current_x(tracker::CoreTracker) = current_x(tracker.state)
current_x(state::CTS) = state.x


"""
    LinearAlgebra.cond(tracker::CoreTracker)

Returns the currently computed approximation of the condition number of the
Jacobian.
"""
LA.cond(tracker::CoreTracker) = cond(tracker.state)
LA.cond(state::CTS) = cond(state.jacobian)


"""
    norm(tracker::CoreTracker)

Returns the norm used to compute distances during the path tracking.
"""
LA.norm(tracker::CoreTracker) = LA.norm(tracker.state)
LA.norm(state::CTS) = state.norm

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
    set_accuracy!(tracker::CoreTracker, accuracy)

Set the current accuracy to `accuracy`.
"""
function set_accuracy!(tracker::CoreTracker, accuracy)
    @unpack options = tracker
    options.accuracy = accuracy
    tracker
end


"""
    max_corrector_iters(tracker::CoreTracker)

Current correction max_steps.
"""
max_corrector_iters(T::CT) = T.options.max_corrector_iters

"""
 set_max_corrector_iters!(tracker::CoreTracker, n)

Set the correction max_steps to `n`.
"""
set_max_corrector_iters!(T::CT, n) = T.options.max_corrector_iters = n

"""
    max_step_size (tracker::CoreTracker)

Current maximal step size.
"""
max_step_size(T::CT) = T.options.max_step_size

"""
    set_max_corrector_iters!(tracker::CoreTracker, Δs)

Set the maximal step size to `Δs`.
"""
set_max_step_size!(T::CT, Δs) = T.options.max_step_size = Δs

## Deprecated ##
"""
    refinement_accuracy(tracker::CoreTracker)

Current refinement accuracy.
"""
function refinement_accuracy(tracker::CoreTracker)
    @warn("Solutions are now automatically refined to maximal achievable accuracy.")
    return tracker.options.accuracy
end

"""
    set_max_refinement_iters!(tracker::CoreTracker, accuracy)

Set the current refinement accuracy to `accuracy`.
"""
function set_refinement_accuracy!(T::CT, accuracy)
    @warn(
        "`set_refinement_accuracy!` is deprecated. Solutions are now automatically refined to maximal achievable accuracy.",
    )
end

"""
    max_refinement_iters(tracker::CoreTracker)

Current refinement max_steps.
"""
function max_refinement_iters(T::CT)
    @warn(
        "`max_refinement_iters` is deprecated. Solutions are now automatically refined to the maximal achievable accuracy.",
    )
    T.options.max_corrector_iters
end

"""
    set_max_refinement_iters!(tracker::CoreTracker, n)

Set the current refinement max_steps to `n`.
"""
function set_max_refinement_iters!(T::CT, n)
    @warn(
        "`set_max_refinement_iters` is deprecated. Solutions are now automatically refined to the maximal achievable accuracy.",
    )
end


"""
    start_parameters!(tracker::CoreTracker, p)

Set the start parameters of the homotopy in in `tracker` to `p`.
"""
start_parameters!(tracker::CoreTracker, p::AbstractVector) =
    (set_start_parameters!(basehomotopy(tracker.homotopy), p); tracker)

"""
    target_parameters!(tracker::CoreTracker, p)

Set the target parameters of the homotopy in in `tracker` to `p`.
"""
target_parameters!(tracker::CoreTracker, p::AbstractVector) =
    (set_target_parameters!(basehomotopy(tracker.homotopy), p); tracker)


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
    iter.tracker.state.status != CoreTrackerStatus.tracking && return nothing

    t = current_t(iter.tracker)
    step_done = false
    while !step_done && is_tracking(iter.tracker.state.status)
        step!(iter.tracker)
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
We want to construct a path tracker to track a parametrized system `f` with parameters `p`
from the parameters `a` to `b`.
```julia
tracker = coretracker(f; parameters = p, start_parameters = a, target_parameters = b)
```
You then can obtain a single solution at `b` by using
```julia
x_b = track(tracker, x_a).x
```

### Trace a path
To trace a path you can use the [`iterator`](@ref) method.

```julia
tracker = coretracker(f; parameters = p, start_parameters = a, target_parameters = b)
for (x, t) in iterator(tracker, x₁)
    @show (x,t)
end
```

If we want to guarantee smooth traces we can limit the maximal step size.
```julia
tracker = coretracker(f; parameters = p, max_step_size = 0.01,
                         start_parameters = a, target_parameters = b)
for (x, t) in iterator(tracker, x₁)
    @show (x,t)
end
```
"""
coretracker(args...; kwargs...) = first(coretracker_startsolutions(args...; kwargs...))
