export PathTracker, PathTrackerResult, PathTrackerStatus,
        pathtracker, pathtracker_startsolutions, affine_tracking,
        track, track!, setup!, iterator,
        currx, currt, currΔt, curriters, currstatus, accuracy, max_corrector_iters,
        max_step_size , refinement_accuracy, max_refinement_iters,
        set_accuracy!, set_max_corrector_iters!, set_refinement_accuracy!,
        set_max_refinement_iters!, set_max_step_size!

const pathtracker_allowed_keywords = [:corrector, :predictor,
    :initial_step_size, :min_step_size , :max_step_size,
    :accuracy, :refinement_accuracy, :max_corrector_iters, :max_refinement_iters,
    :max_steps, :simple_step_size_alg,
    # deprecated
    :initial_steplength, :minimal_steplength, :maximal_steplength,
    :initial_step_size, :minimal_step_size, :maximal_step_size,
    :tol, :refinement_tol, :corrector_maxiters, :refinement_maxiters,
    :maxiters, :simple_step_size]


####################
# PathTrackerState #
####################

module PathTrackerStatus

    @doc """
        PathTrackerStatus.states

    The possible states the pathtracker can achieve are

    * `PathTrackerStatus.success`
    * `PathTrackerStatus.tracking`
    * `PathTrackerStatus.terminated_maximal_iterations`
    * `PathTrackerStatus.terminated_invalid_startvalue`
    * `PathTrackerStatus.terminated_step_size_too_small`
    * `PathTrackerStatus.terminated_singularity`
    * `PathTrackerStatus.terminated_ill_conditioned`
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
# PathTrackerResult #
####################

"""
     PathTrackerResult{V<:AbstractVector}

Containing the result of a tracked path. The fields are
* `returncode::PathTrackerStatus.states` If the tracking was successfull then it is `PathTrackerStatus.success`.
* `x::V` The result.
* `t::ComplexF64` The `t` when the path tracker stopped.
* `accuracy::Float64`: The estimated accuracy of `x`.
"""
struct PathTrackerResult{V <: AbstractVector}
     returncode::PathTrackerStatus.states
     x::V
     t::ComplexF64
     accuracy::Float64
     accepted_steps::Int
     rejected_steps::Int
end

function PathTrackerResult(tracker)
    state = tracker.state
     PathTrackerResult(state.status,
          copy(state.x), state.segment[state.s],
          state.accuracy,
          state.accepted_steps,
          state.rejected_steps)
end

Base.show(io::IO, result::PathTrackerResult) = print_fieldnames(io, result)
Base.show(io::IO, ::MIME"application/prs.juno.inline", result::PathTrackerResult) = result


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


######################
# PathTrackerOptions #
######################
mutable struct PathTrackerOptions
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
    maximal_lost_digits::Float64
    auto_scaling::Bool
    auto_scaling_options::AutoScalingOptions
end



function PathTrackerOptions(::Type{Precision}; accuracy=1e-7,
    refinement_accuracy=1e-8,
    max_corrector_iters::Int=2,
    max_refinement_iters=max_corrector_iters,
    max_steps=1_000,
    initial_step_size=0.1,
    min_step_size=1e-14,
    max_step_size=Inf,
    simple_step_size_alg=false,
    update_patch=true,
    maximal_lost_digits=default_maximal_lost_digits(Precision),
    auto_scaling=true,
    auto_scaling_options=AutoScalingOptions(),
    # deprecated in 0.6
    tol=nothing,
    refinement_tol=nothing,
    corrector_maxiters=nothing,
    refinement_maxiters=nothing,
    maxiters=nothing,
    initial_steplength=nothing,
    minimal_steplength=nothing,
    maximal_steplength=nothing,
    simple_step_size=nothing
    ) where {Precision<:Real}

    @deprecatekwarg tol accuracy
    @deprecatekwarg refinement_tol refinement_accuracy
    @deprecatekwarg corrector_maxiters max_corrector_iters
    @deprecatekwarg refinement_maxiters max_refinement_iters
    @deprecatekwarg maxiters max_steps
    @deprecatekwarg initial_steplength initial_step_size
    @deprecatekwarg minimal_steplength min_step_size
    @deprecatekwarg maximal_steplength max_step_size
    @deprecatekwarg simple_step_size simple_step_size_alg

    PathTrackerOptions(accuracy, max_corrector_iters, refinement_accuracy,
            max_refinement_iters, max_steps, initial_step_size, min_step_size,
            max_step_size, simple_step_size_alg, update_patch, float(maximal_lost_digits),
            auto_scaling, auto_scaling_options)
end

default_maximal_lost_digits(::Type{T}) where T = -log10(eps(T)) - 3

Base.show(io::IO, opts::PathTrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::PathTrackerOptions) = opts


####################
# PathTrackerState #
####################

mutable struct PathTrackerState{T, AV<:AbstractVector{T}, MaybePatchState <: Union{AbstractAffinePatchState, Nothing}, IP}
    x::AV # current x
    x̂::AV # last prediction
    x̄::AV # canidate for new x
    ẋ::Vector{T} # derivative at current x
    inner_product::IP
    η::Float64
    ω::Float64
    segment::ComplexSegment
    s::Float64 # current step length (0 ≤ s ≤ length(segment))
    Δs::Float64 # current step size
    Δs_prev::Float64 # previous step size
    accuracy::Float64
    # The relative number of digits lost during the solution of the linear systems
    # in Newton's method. See `solve_with_digits_lost!` in utilities/linear_algebra.jl
    # for how this is computed.
    digits_lost::Float64
    status::PathTrackerStatus.states
    patch::MaybePatchState
    accepted_steps::Int
    rejected_steps::Int
    last_step_failed::Bool
    consecutive_successfull_steps::Int
end

function PathTrackerState(x₁::AbstractVector, t₁, t₀, options::PathTrackerOptions,
                           patch::Union{Nothing, AbstractAffinePatchState}=nothing,
                           inner_product::AbstractInnerProduct=EuclideanIP())
    if isa(x₁, SVector)
        x = Vector(x₁)
    else
        x = copy(x₁)
    end
    x̂, x̄ = copy(x), copy(x)
    ẋ = Vector(x)
    η = ω = NaN
    segment = ComplexSegment(promote(complex(t₁), complex(t₀))...)
    s = 0.0
    Δs = convert(Float64, min(options.initial_step_size, length(segment), options.max_step_size ))
    Δs_prev = 0.0
    accuracy = 0.0
    digits_lost = 0.0
    accepted_steps = rejected_steps = 0
    status = PathTrackerStatus.tracking
    last_step_failed = false
    consecutive_successfull_steps = 0
    PathTrackerState(x, x̂, x̄, ẋ, inner_product, η, ω, segment, s, Δs, Δs_prev, accuracy, digits_lost, status, patch,
        accepted_steps, rejected_steps, last_step_failed, consecutive_successfull_steps)
end

Base.show(io::IO, state::PathTrackerState) = print_fieldnames(io, state)
Base.show(io::IO, ::MIME"application/prs.juno.inline", state::PathTrackerState) = state

function reset!(state::PathTrackerState, x₁::AbstractVector, t₁, t₀, options::PathTrackerOptions, setup_patch)
    state.segment = ComplexSegment(promote(t₁, t₀)...)
    state.s = 0.0
    state.Δs = min(options.initial_step_size, length(state.segment), options.max_step_size )
    state.Δs_prev = 0.0
    state.accuracy = 0.0
    state.accepted_steps = state.rejected_steps = 0
    state.status = PathTrackerStatus.tracking
    embed!(state.x, x₁)
    setup_patch && state.patch !== nothing && setup!(state.patch, state.x)
    state.η = state.ω = NaN
    if options.auto_scaling
        init_auto_scaling!(state.inner_product, state.x, options.auto_scaling_options)
    end
    state.last_step_failed = false
    state.consecutive_successfull_steps = 0
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
            wᵢ = max(opts.scale_min * point_norm, opts.scale_abs_min)
        elseif wᵢ > opts.scale_max * point_norm
            wᵢ = opts.scale_max * point_norm
        end
        w[i] = wᵢ
    end
    nothing
end
init_auto_scaling!(ip::EuclideanIP, x::AbstractVector, opts::AutoScalingOptions) = nothing

####################
# PathTrackerCache #
####################
mutable struct PathTrackerCache{H<:HomotopyWithCache, P<:AbstractPredictorCache,
             C<:AbstractCorrectorCache, T, F}
    homotopy::H
    predictor::P
    corrector::C
    Jac::Jacobian{T, F}
    out::Vector{T}
    r::Vector{T}
end
function PathTrackerCache(H::HomotopyWithCache, predictor, corrector, state::PathTrackerState)
    t = state.segment[state.s]
    pcache = cache(predictor, H, state.x, state.ẋ, t)
    ccache = cache(corrector, H, state.x, t)
    Jac = Jacobian(jacobian(H, state.x, t))
    out = H(state.x, t)
    r = copy(out)
    PathTrackerCache(H, pcache, ccache, Jac, out, r)
end


###############
# PathTracker #
###############
"""
     PathTracker(H::AbstractHomotopy, x₁, t₁, t₀; options...)::PathTracker

Create a `PathTracker` to track `x₁` from `t₁` to `t₀`. The homotopy `H`
needs to be homogenous. Note that a `PathTracker` is also a (mutable) iterator.

## PathTrackerOptions
* `corrector::AbstractCorrector`: The corrector used during in the predictor-corrector scheme. The default is [`NewtonCorrector`](@ref).
* `max_corrector_iters=3`: The maximal number of correction steps in a single step.
* `initial_step_size=0.1`: The step size of the first step.
* `max_steps=10_000`: The maximal number of iterations the path tracker has available.
* `min_step_size =1e-14`: The minimal step size.
* `max_step_size =Inf`: The maximal step size.
* `maximal_lost_digits::Real=-(log₁₀(eps) + 3)`: The tracking is terminated if we estimate that we loose more than `maximal_lost_digits` in the linear algebra steps.
* `predictor::AbstractPredictor`: The predictor used during in the predictor-corrector scheme. The default is [`Heun`](@ref)()`.
* `max_refinement_iters=max_corrector_iters`: The maximal number of correction steps used to refine the final value.
* `refinement_accuracy=1e-8`: The precision used to refine the final value.
* `accuracy=1e-7`: The precision used to track a value.
* `auto_scaling=true`: This only applies if we track in affine space. Automatically regauges the variables to effectively compute with a relative accuracy instead of an absolute one.
"""
struct PathTracker{H<:AbstractHomotopy,
    Predictor<:AbstractPredictor,
    Corrector<:AbstractCorrector,
    MaybePatch<:Union{Nothing, AbstractAffinePatch},
    S<:PathTrackerState,
    C<:PathTrackerCache}
    # these are fixed
    homotopy::H
    predictor::Predictor
    corrector::Corrector
    affine_patch::MaybePatch
    # these are mutable
    state::S
    options::PathTrackerOptions
    cache::C
end

"""
    PathTracker(problem::AbstractProblem, x₁, t₁, t₀; kwargs...)

Construct a [`PathTracker`](@ref) from the given `problem`.
"""
function PathTracker(prob::AbstractProblem, x₁, t₁, t₀; kwargs...)
    PathTracker(prob.homotopy, embed(prob, x₁), t₁, t₀; kwargs...)
end

# Tracking in Projective Space
function PathTracker(H::AbstractHomotopy, x₁::ProjectiveVectors.PVector, t₁, t₀;
    patch=OrthogonalPatch(),
    corrector::AbstractCorrector=NewtonCorrector(),
    predictor::AbstractPredictor=default_predictor(x₁), kwargs...)

    options = PathTrackerOptions(real(eltype(x₁)); kwargs...)
    # disable auto scaling always in projective space
    options.auto_scaling = false

    if H isa PatchedHomotopy
        error("You cannot pass a `PatchedHomotopy` to PathTracker. Instead pass the homotopy and patch separate.")
    end

    patch_state = state(patch, x₁)
    # We close over the patch state, the homotopy and its cache
    # to be able to pass things around more easily
    HC = HomotopyWithCache(PatchedHomotopy(H, patch_state), x₁, t₁)

    # We have to make sure that the element type of x is invariant under evaluation
    tracker_state = PathTrackerState(indempotent_x(HC, x₁, t₁), t₁, t₀, options, patch_state)
    cache = PathTrackerCache(HC, predictor, corrector, tracker_state)

    PathTracker(H, predictor, corrector, patch, tracker_state, options, cache)
end

# Tracking in affine space
function PathTracker(H::AbstractHomotopy, x₁::AbstractVector, t₁, t₀;
    corrector::AbstractCorrector=NewtonCorrector(),
    predictor::AbstractPredictor=default_predictor(x₁), kwargs...)

    options = PathTrackerOptions(real(eltype(x₁)); kwargs...)
    # We close over the patch state, the homotopy and its cache
    # to be able to pass things around more easily
    HC = HomotopyWithCache(H, x₁, t₁)

    inner_product = options.auto_scaling ? WeightedIP(x₁) : EuclideanIP()
    # We have to make sure that the element type of x is invariant under evaluation
    tracker_state = PathTrackerState(indempotent_x(HC, x₁, t₁), t₁, t₀, options, nothing, inner_product)
    cache = PathTrackerCache(HC, predictor, corrector, tracker_state)

    PathTracker(H, predictor, corrector, nothing, tracker_state, options, cache)
end

default_predictor(x::AbstractVector) = Heun()
# Do not really understand this but Heun doesn't work that great for multi-homogenous tracking
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

Base.show(io::IO, ::PathTracker) = print(io, "PathTracker()")
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::PathTracker) = x

##############
## TRACKING ##
##############

"""
    track(tracker, x₁, t₁=1.0, t₀=0.0; options...)::PathTrackerResult

Track a value `x₁` from `t₁` to `t₀` using the given `PathTracker` `tracker`.
This returns a `PathTrackerResult`. This modifies `tracker`.
See [`track!`](@ref) for the possible options.
"""
function track(tracker::PathTracker, x₁::AbstractVector, t₁=1.0, t₀=0.0; kwargs...)
     track!(tracker, x₁, t₁, t₀; kwargs...)
     PathTrackerResult(tracker)
end

"""
     track!(tracker, x₁, t₁=1.0, t₀=0.0; setup_patch=true, checkstartvalue=true, compute_ẋ=true)

Track a value `x₁` from `t₁` to `t₀` using the given `PathTracker` `tracker`.
Returns one of the enum values of `PathTrackerStatus.states` indicating the status.
If the tracking was successfull it is `PathTrackerStatus.success`.
If `setup_patch` is `true` then [`setup!`](@ref) is called at the beginning
of the tracking.

    track!(x₀, tracker, x₁, t₁=1.0, t₀=0.0; options...)

Additionally also stores the result in `x₀` if the tracking was successfull.
"""
function track!(x₀, tracker::PathTracker, x₁, t₁=1.0, t₀=0.0; setup_patch=tracker.options.update_patch, checkstartvalue=true, compute_ẋ=true)
     track!(tracker, x₁, t₁, t₀, setup_patch, checkstartvalue, compute_ẋ)
     retcode = currstatus(tracker)
     if retcode == PathTrackerStatus.success
         x₀ .= currx(tracker)
     end
     retcode
end
@inline function track!(tracker::PathTracker, x₁, t₁=1.0, t₀=0.0; setup_patch=tracker.options.update_patch, checkstartvalue=true, compute_ẋ=true)
    track!(tracker, x₁, t₁, t₀, setup_patch, checkstartvalue, compute_ẋ)
end
function track!(tracker::PathTracker, x₁, t₁, t₀, setup_patch, checkstartvalue=true, compute_ẋ=true)
    setup!(tracker, x₁, t₁, t₀, setup_patch, checkstartvalue, compute_ẋ)

    while tracker.state.status == PathTrackerStatus.tracking
        step!(tracker)
        check_terminated!(tracker)
    end
    if tracker.state.status == PathTrackerStatus.success
        refine!(tracker)
    end

    tracker.state.status
end

"""
    setup!(pathtracker, x₁, t₁=1.0, t₀=0.0, setup_patch=pathtracker.options.update_patch, checkstartvalue=true, compute_ẋ=true)

Setup `pathtracker` to track `x₁` from `t₁` to `t₀`. Use this if you want to use the
pathtracker as an iterator.
"""
function setup!(tracker::PathTracker, x₁::AbstractVector, t₁=1.0, t₀=0.0, setup_patch=tracker.options.update_patch, checkstartvalue=true, compute_ẋ=true)
    state, cache = tracker.state, tracker.cache

    try
        reset!(state, x₁, t₁, t₀, tracker.options, setup_patch)
        reset!(cache.predictor, state.x, t₁)
        checkstartvalue && checkstartvalue!(tracker)
        if compute_ẋ
            compute_ẋ!(state, cache, tracker.options)
            setup!(cache.predictor, cache.homotopy, state.x, state.ẋ, currt(state), cache.Jac)
        end
    catch err
        if !(err isa LinearAlgebra.SingularException)
            rethrow(err)
        end
        tracker.state.status = PathTrackerStatus.terminated_singularity
    end
    tracker
end

function checkstartvalue!(tracker::PathTracker)
    result = correct!(tracker.state.x̄, tracker)
    if isconverged(result)
        tracker.state.x .= tracker.state.x̄
    else
        tracker.state.status = PathTrackerStatus.terminated_invalid_startvalue
    end
    nothing
end

function compute_ẋ!(state, cache, options::PathTrackerOptions)
    @inbounds jacobian_and_dt!(cache.Jac.J, cache.out, cache.homotopy, state.x, currt(state))
    # apply row scaling to J and compute factorization
    updated_jacobian!(cache.Jac)

    @inbounds for i in eachindex(cache.out)
        cache.out[i] = -cache.out[i]
    end
    solve!(state.ẋ, cache.Jac, cache.out)
    nothing
end


function correct!(x̄, tracker::PathTracker,
        x=tracker.state.x,
        t=tracker.state.segment[tracker.state.s],
        norm=tracker.state.inner_product;
        accuracy=tracker.options.accuracy,
        max_iters=tracker.options.max_corrector_iters)

    correct!(x̄, tracker.corrector, tracker.cache.corrector,
             tracker.cache.homotopy, x, t, norm, accuracy, max_iters)
end

function step!(tracker::PathTracker)
    state, cache, options = tracker.state, tracker.cache, tracker.options
    H = cache.homotopy
    x, x̂, x̄, ẋ = state.x, state.x̂, state.x̄, state.ẋ

    try
        t, Δt = currt(state), currΔt(state)
        predict!(x̂, tracker.predictor, cache.predictor, H, x, t, Δt, ẋ)
        result = correct!(x̄, tracker, x̂, t + Δt)
        if isconverged(result)
            # Step is accepted, assign values
            state.accepted_steps += 1
            x .= x̄
            state.s += state.Δs
            state.accuracy = result.accuracy
            state.digits_lost = result.digits_lost
            # Step size change
            update_stepsize!(state, result, order(tracker.predictor), options)
            if state.patch !== nothing && options.update_patch
                changepatch!(state.patch, x)
            end
            options.auto_scaling && auto_scaling!(state, options)
            # update derivative
            compute_ẋ!(state, cache, options)
            # tell the predictors about the new derivative if they need to update something
            update!(cache.predictor, H, x, ẋ, t + Δt, cache.Jac)
        else
            # We have to reset the patch
            state.rejected_steps += 1
            state.digits_lost = result.digits_lost
            # Step failed, so we have to try with a new (smaller) step size
            update_stepsize!(state, result, order(tracker.predictor), options)
            Δt = currΔt(state)
            state.last_step_failed = true
            # Check termination criteria
            # 1) Step size get's too small:
            if state.Δs < options.min_step_size
                state.status = PathTrackerStatus.terminated_step_size_too_small
            end
            # 2) We became too ill-conditioned
            if !is_overdetermined_tracking(tracker) &&
               state.digits_lost > options.maximal_lost_digits
                state.status = PathTrackerStatus.terminated_ill_conditioned
            end
        end
    catch err
        if !(err isa LinearAlgebra.SingularException)
            rethrow(err)
        end
        tracker.state.status = PathTrackerStatus.terminated_singularity
    end
    nothing
end

function is_overdetermined_tracking(tracker)
    m, n = size(tracker.cache.homotopy)
    m > n
end

g(Θ) = sqrt(1+4Θ) - 1
# Choose 0.25 instead of 1.0 due to Newton-Kantorovich theorem
δ(opts::PathTrackerOptions, ω) = @fastmath min(√(ω/2) * τ(opts), 0.25)
τ(opts::PathTrackerOptions) = nthroot(opts.accuracy, 2 * opts.max_corrector_iters)

function update_stepsize!(state::PathTrackerState, result::CorrectorResult,
                          order::Int, options::PathTrackerOptions)

    if options.simple_step_size_alg
        simple_step_size_alg!(state, result, options)
        return nothing
    end

    # The step size control is described in https://arxiv.org/abs/1902.02968

    # we have to handle the special case that there is only 1 iteration
    # in this case we cannot estimate ω and therefore just assume ω = 2
    # Also note ||x̂-x|| = ||Δx₀||
    if result.iters == 1
        ω = 2.0
        d_x̂_x̄ = result.norm_Δx₀
    else
        ω = result.ω
        d_x̂_x̄ = distance(state.x̂, state.x, state.inner_product)
    end

    Δx₀ = result.norm_Δx₀
    if isconverged(result)
        # compute η and update
        η = d_x̂_x̄ / state.Δs^order

        # assume Δs′ = Δs
        ω′ = isnan(state.ω) ? ω : max(2ω - state.ω, ω)
        if isnan(state.η)
            Δs′ = state.Δs
        else
            d_x̂_x̄′ = max(2d_x̂_x̄ - state.η * state.Δs^(order), 0.75d_x̂_x̄)
            if state.last_step_failed
                d_x̂_x̄′ *= 2
            end
            λ = g(δ(options, ω′)) / (ω′ * d_x̂_x̄′)
            Δs′ = 0.8 * nthroot(λ, order) * state.Δs
        end
        if state.last_step_failed
            Δs′ = min(Δs′, state.Δs)
        end
        state.η = η
        state.ω = ω
        state.last_step_failed = false
    else
        ω = max(ω, state.ω)
        δ_N_ω = δ(options, ω)
        ω_η = ω / 2 * Δx₀
        if δ_N_ω < ω_η
            λ = g(δ_N_ω) / g(ω_η)
        elseif δ_N_ω < 2ω_η
            λ = g(δ_N_ω) / g(2ω_η)
        elseif δ_N_ω < 4ω_η
            λ = g(δ_N_ω) / g(4ω_η)
        elseif δ_N_ω < 8ω_η
            λ = g(δ_N_ω) / g(8ω_η)
        else
            λ = 0.5^order
        end
        Δs′ = 0.9 * nthroot(λ, order) * state.Δs
        state.last_step_failed = true
    end

    state.Δs = min(Δs′, length(state.segment) - state.s, options.max_step_size )

    if !isconverged(result) && state.Δs < options.min_step_size
        state.status = PathTrackerStatus.terminated_step_size_too_small
    end
    nothing
end

function simple_step_size_alg!(state::PathTrackerState, result::CorrectorResult, options::PathTrackerOptions)
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

    if !isconverged(result) && state.Δs < options.min_step_size
        state.status = PathTrackerStatus.terminated_step_size_too_small
    end
end

function auto_scaling!(state::PathTrackerState, options::PathTrackerOptions)
    auto_scaling!(state.inner_product, state.x, options.auto_scaling_options)
end
function auto_scaling!(ip::WeightedIP, x::AbstractVector, opts::AutoScalingOptions)
    norm_x = ip(x)
    for i in 1:length(x)
        wᵢ = (abs(x[i]) + ip.weight[i]) / 2
        if wᵢ < opts.scale_min * norm_x
            wᵢ = max(opts.scale_min * norm_x, opts.scale_abs_min)
        elseif wᵢ > opts.scale_max * norm_x
            wᵢ = opts.scale_max * norm_x
        end
        ip.weight[i] = wᵢ
    end
    nothing
end
auto_scaling!(ip::EuclideanIP, x::AbstractVector, opts::AutoScalingOptions) = nothing

function check_terminated!(tracker::PathTracker)
    if abs(tracker.state.s - length(tracker.state.segment)) < 2eps(length(tracker.state.segment))
        tracker.state.status = PathTrackerStatus.success
    elseif curriters(tracker) ≥ tracker.options.max_steps
        tracker.state.status = PathTrackerStatus.terminated_maximal_iterations
    end
    nothing
end

function refine!(tracker::PathTracker)
    if tracker.state.accuracy < tracker.options.refinement_accuracy
        return
    end
    result = correct!(tracker.state.x̄, tracker;
        accuracy=tracker.options.refinement_accuracy,
        max_iters=tracker.options.max_refinement_iters)
    if isconverged(result)
        tracker.state.x .= tracker.state.x̄
        tracker.state.accuracy = result.accuracy
    end
    nothing
end

# TODO: REMOVE THIS
function residual(tracker::PathTracker, x, t)
    evaluate!(tracker.cache.out, tracker.cache.homotopy, x, t)
    infinity_norm(tracker.cache.out)
end

function Base.iterate(tracker::PathTracker, state=nothing)
    state === nothing && return tracker, 1

    if tracker.state.status == PathTrackerStatus.tracking
        step!(tracker)
        check_terminated!(tracker)

        if tracker.state.status == PathTrackerStatus.success
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
affine_tracking(tracker::PathTracker) = !isa(tracker.state.x, ProjectiveVectors.PVector)

#################
## Query State ##
#################
"""
     currt(tracker::PathTracker)

Current `t`.
"""
currt(tracker::PathTracker) = currt(tracker.state)
currt(state::PathTrackerState) = state.segment[state.s]

"""
     currΔt(tracker::PathTracker)

Current step_size `Δt`.
"""
currΔt(tracker::PathTracker) = currΔt(tracker.state)
currΔt(state::PathTrackerState) = state.segment[state.Δs] - state.segment.start

"""
     curriters(tracker::PathTracker)

Current number of iterations.
"""
curriters(tracker::PathTracker) = curriters(tracker.state)
curriters(state::PathTrackerState) = state.accepted_steps + state.rejected_steps

"""
     currstatus(tracker::PathTracker)

Current status.
"""
currstatus(tracker::PathTracker) = currstatus(tracker.state)
currstatus(state::PathTrackerState) = state.status

"""
    currx(tracker::PathTracker)

Return the current value of `x`.
"""
currx(tracker::PathTracker) = currx(tracker.state)
currx(state::PathTrackerState) = state.x


##################
# Modify options #
##################
"""
     accuracy(tracker::PathTracker)

Current accuracy.
"""
accuracy(tracker::PathTracker) = tracker.options.accuracy

"""
     set_accuracy!(tracker::PathTracker, accuracy)

Set the current accuracy to `accuracy`.
"""
set_accuracy!(tracker::PathTracker, accuracy) = tracker.options.accuracy = accuracy

"""
     refinement_accuracy(tracker::PathTracker)

Current refinement accuracy.
"""
refinement_accuracy(tracker::PathTracker) = tracker.options.refinement_accuracy

"""
     set_max_refinement_iters!(tracker::PathTracker, accuracy)

Set the current refinement accuracy to `accuracy`.
"""
set_refinement_accuracy!(T::PathTracker, accuracy) = T.options.refinement_accuracy = accuracy

"""
     max_refinement_iters(tracker::PathTracker)

Current refinement max_steps.
"""
max_refinement_iters(T::PathTracker) = T.options.max_refinement_iters

"""
     set_max_refinement_iters!(tracker::PathTracker, n)

Set the current refinement max_steps to `n`.
"""
set_max_refinement_iters!(T::PathTracker, n) = T.options.max_refinement_iters = n

"""
     max_corrector_iters(tracker::PathTracker)

Current correction max_steps.
"""
max_corrector_iters(T::PathTracker) = T.options.max_corrector_iters

"""
     set_max_corrector_iters!(tracker::PathTracker, n)

Set the correction max_steps to `n`.
"""
set_max_corrector_iters!(T::PathTracker, n) = T.options.max_corrector_iters = n

"""
     max_step_size (tracker::PathTracker)

Current maximal step size.
"""
max_step_size(T::PathTracker) = T.options.max_step_size

"""
     set_max_corrector_iters!(tracker::PathTracker, Δs)

Set the maximal step size to `Δs`.
"""
set_max_step_size!(T::PathTracker, Δs) = T.options.max_step_size  = Δs

################
# PathIterator #
################
struct PathIterator{Tracker<:PathTracker}
    tracker::Tracker
    x_affine::Bool
    t_real::Bool
end
Base.IteratorSize(::Type{<:PathIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:PathIterator}) = Base.HasEltype()

"""
    iterator(tracker::PathTracker, x₁, t₁=1.0, t₀=0.0; affine=true)

Prepare a tracker to make it usable as a (stateful) iterator. Use this if you want to inspect a specific
path. In each iteration the tuple `(x,t)` is returned.
If `affine == true` then `x` is the affine solution (internally we compute in projective space).

## Example

Assume you have `PathTracker` `tracker` and you wan to track `x₁` from 1.0 to 0.25:
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
println("Success: ", currstatus(tracker) == PathTrackerStatus.success)
```
"""
function iterator(tracker::PathTracker, x₁, t₁=1.0, t₀=0.0; affine=true, kwargs...)
    setup!(tracker, x₁, t₁, t₀; kwargs...)
    PathIterator(tracker, affine, typeof(t₁ - t₀) <: Real)
end

function current_x_t(iter::PathIterator)
    x = currx(iter.tracker)
    t = currt(iter.tracker)
    (iter.x_affine ? ProjectiveVectors.affine_chart(x) : x,
     iter.t_real ? real(t) : t)
end

function Base.iterate(iter::PathIterator, state=nothing)
    state === nothing && return current_x_t(iter), 1
    iter.tracker.state.status != PathTrackerStatus.tracking && return nothing

    step_done = false
    while !step_done && (iter.tracker.state.status == PathTrackerStatus.tracking)

        step!(iter.tracker)
        check_terminated!(iter.tracker)

        if iter.tracker.state.status == PathTrackerStatus.success
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
    pathtracker_startsolutions(args...; kwargs...)

Construct a [`PathTracker`](@ref) and `startsolutions` in the same way `solve`
does it. This also takes the same input arguments as `solve`. This is convenient if you want
to investigate single paths.
"""
function pathtracker_startsolutions(args...; kwargs...)
    supported, rest = splitkwargs(kwargs,problem_startsolutions_supported_keywords)
    prob, startsolutions = problem_startsolutions(args...; supported...)
    tracker = PathTracker(prob, start_solution_sample(startsolutions), one(ComplexF64), zero(ComplexF64); rest...)

    (tracker=tracker, startsolutions=startsolutions)
end

"""
    pathtracker(args...; kwargs...)

Construct a [`PathTracker`](@ref) in the same way `solve`
does it. This also takes the same input arguments as `solve` with the exception that you do not need to specify startsolutions.
This is convenient if you want to investigate single paths.

## Examples

### Obtain single solution
We want to construct a path tracker to track a parameterized system `f` with parameters `p`
from the parameters `a` to `b`.
```julia
tracker = pathtracker(f, parameters=p, p₁=a, p₀=b)
```
You then can obtain a single solution at `b` by using
```julia
x_b = track(tracker, x_a).x
```

### Trace a path
To trace a path you can use the [`iterator`](@ref) method.

```julia
tracker = pathtracker(f, parameters=p, p₁=a, p₀=b, max_step_size =0.01)
for (x, t) in iterator(tracker, x₁)
    @show (x,t)
end
```

If we want to guarantee smooth traces we can limit the maximal step size.
```julia
tracker = pathtracker(f, parameters=p, p₁=a, p₀=b, max_step_size =0.01)
for (x, t) in iterator(tracker, x₁)
    @show (x,t)
end
```
"""
function pathtracker(args...; kwargs...)
    tracker, _ = pathtracker_startsolutions(args...; kwargs...)
    tracker
end
