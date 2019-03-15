export PathTracker, PathTrackerResult, PathTrackerStatus,
        pathtracker, pathtracker_startsolutions,
        track, track!, setup!, iterator,
        currx, currt, currΔt, curriters, currstatus, tol, corrector_maxiters,
        maximal_step_size, refinement_tol, refinement_maxiters,
        set_tol!, set_corrector_maxiters!, set_refinement_tol!,
        set_refinement_maxiters!, set_maximal_step_size!

const pathtracker_allowed_keywords = [:corrector, :predictor, :initial_steplength,
    :minimal_steplength, :maximal_steplength,
    :initial_step_size, :minimal_step_size, :maximal_step_size,
    :tol, :refinement_tol, :corrector_maxiters,  :refinement_maxiters,
    :maxiters, :simple_step_size]


####################
# PathTrackerState #
####################

module PathTrackerStatus
    """
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


######################
# PathTrackerOptions #
######################
mutable struct PathTrackerOptions
    tol::Float64
    corrector_maxiters::Int
    refinement_tol::Float64
    refinement_maxiters::Int
    maxiters::Int
    initial_step_size::Float64
    minimal_step_size::Float64
    maximal_step_size::Float64
    simple_step_size::Bool
    update_patch::Bool
    maximal_lost_digits::Float64
end

function PathTrackerOptions(::Type{Precision};tol=1e-7,
    refinement_tol=1e-8,
    corrector_maxiters::Int=2,
    refinement_maxiters=corrector_maxiters,
    maxiters=10_000,
    initial_steplength=0.1,
    initial_step_size=initial_steplength,
    minimal_steplength=1e-14,
    minimal_step_size=minimal_steplength,
    maximal_steplength=Inf,
    maximal_step_size=maximal_steplength,
    simple_step_size=false,
    update_patch=true,
    maximal_lost_digits=default_maximal_lost_digits(Precision)) where {Precision<:Real}

    PathTrackerOptions(tol, corrector_maxiters, refinement_tol, refinement_maxiters, maxiters,
            initial_step_size, minimal_step_size, maximal_step_size, simple_step_size, update_patch,
            float(maximal_lost_digits))
end

default_maximal_lost_digits(::Type{T}) where T = -log10(eps(T)) - 3

Base.show(io::IO, opts::PathTrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::PathTrackerOptions) = opts


####################
# PathTrackerState #
####################

mutable struct PathTrackerState{T, AV<:AbstractVector{T}, MaybePatchState <: Union{AbstractAffinePatchState, Nothing}}
    x::AV # current x
    x̂::AV # last prediction
    x̄::AV # canidate for new x
    ẋ::Vector{T} # derivative at current x
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

function PathTrackerState(x₁::AbstractVector, t₁, t₀, options::PathTrackerOptions, patch::Union{Nothing, AbstractAffinePatchState}=nothing)
    x, x̂, x̄ = copy(x₁), copy(x₁), copy(x₁)
    if x₁ isa PVector
        ẋ = copy(x₁.data)
    else
        ẋ = copy(x₁)
    end
    η = ω = NaN
    segment = ComplexSegment(promote(complex(t₁), complex(t₀))...)
    s = 0.0
    Δs = convert(Float64, min(options.initial_step_size, length(segment), options.maximal_step_size))
    Δs_prev = 0.0
    accuracy = 0.0
    digits_lost = 0.0
    accepted_steps = rejected_steps = 0
    status = PathTrackerStatus.tracking
    last_step_failed = false
    consecutive_successfull_steps = 0
    PathTrackerState(x, x̂, x̄, ẋ, η, ω, segment, s, Δs, Δs_prev, accuracy, digits_lost, status, patch,
        accepted_steps, rejected_steps, last_step_failed, consecutive_successfull_steps)
end

Base.show(io::IO, state::PathTrackerState) = print_fieldnames(io, state)
Base.show(io::IO, ::MIME"application/prs.juno.inline", state::PathTrackerState) = state

function reset!(state::PathTrackerState, x₁::AbstractVector, t₁, t₀, options::PathTrackerOptions, setup_patch)
    state.segment = ComplexSegment(promote(t₁, t₀)...)
    state.η = state.ω = NaN
    state.s = 0.0
    state.Δs = min(options.initial_step_size, length(state.segment), options.maximal_step_size)
    state.Δs_prev = 0.0
    state.accuracy = 0.0
    state.accepted_steps = state.rejected_steps = 0
    state.status = PathTrackerStatus.tracking
    embed!(state.x, x₁)
    setup_patch && state.patch !== nothing && setup!(state.patch, state.x)
    state.last_step_failed = false
    state.consecutive_successfull_steps = 0
    state
end

embed!(x::ProjectiveVectors.PVector, y) = ProjectiveVectors.embed!(x, y)
embed!(x::AbstractVector, y) = x .= y

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
* `corrector_maxiters=3`: The maximal number of correction steps in a single step.
* `initial_step_size=0.1`: The step size of the first step.
* `maxiters=10_000`: The maximal number of iterations the path tracker has available.
* `minimal_step_size=1e-14`: The minimal step size.
* `maximal_step_size=Inf`: The maximal step size.
* `maximal_lost_digits::Real=-(log₁₀(eps) + 3)`: The tracking is terminated if we estimate that we loose more than `maximal_lost_digits` in the linear algebra steps.
* `predictor::AbstractPredictor`: The predictor used during in the predictor-corrector scheme. The default is [`Heun`](@ref)()`.
* `refinement_maxiters=corrector_maxiters`: The maximal number of correction steps used to refine the final value.
* `refinement_tol=1e-8`: The precision used to refine the final value.
* `tol=1e-7`: The precision used to track a value.
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

    # We have to make sure that the element type of x is invariant under evaluation
    tracker_state = PathTrackerState(indempotent_x(HC, x₁, t₁), t₁, t₀, options)
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
    indem_x = similar(x₁, promote_type(typeof(u[1]), ComplexF64))
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


function correct!(x̄, tracker::PathTracker, x=tracker.state.x, t=tracker.state.segment[tracker.state.s];
    tol=tracker.options.tol,
    maxiters=tracker.options.corrector_maxiters)
    correct!(x̄, tracker.corrector, tracker.cache.corrector,
                        tracker.cache.homotopy, x, t, tol, maxiters)
end

function step!(tracker::PathTracker)
    state, cache, options = tracker.state, tracker.cache, tracker.options
    H = cache.homotopy
    x, x̂, x̄, ẋ = state.x, state.x̂, state.x̄, state.ẋ

    try
        t, Δt = currt(state), currΔt(state)
        predict!(x̂, tracker.predictor, cache.predictor, H, x, t, Δt, ẋ)
        result = correct!(x̄, tracker.corrector, cache.corrector, H, x̂, t + Δt, options.tol, options.corrector_maxiters)

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
            if state.Δs < options.minimal_step_size
                state.status = PathTrackerStatus.terminated_step_size_too_small
            end
            # 2) We became too ill-conditioned
            if state.digits_lost > options.maximal_lost_digits
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

g(Θ) = sqrt(1+4Θ) - 1
# Choose 0.25 instead of 1.0 due to Newton-Kantorovich theorem
δ(opts::PathTrackerOptions, ω) = @fastmath min(√(ω/2) * τ(opts), 0.25)
τ(opts::PathTrackerOptions) = nthroot(opts.tol, 2 * opts.corrector_maxiters)

function update_stepsize!(state::PathTrackerState, result::CorrectorResult,
                          order::Int, options::PathTrackerOptions)

    if options.simple_step_size
        simple_step_size!(state, result, options)
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
        d_x̂_x̄ = euclidean_distance(state.x̂, state.x)
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

    state.Δs = min(Δs′, length(state.segment) - state.s, options.maximal_step_size)

    if !isconverged(result) && state.Δs < options.minimal_step_size
        state.status = PathTrackerStatus.terminated_step_size_too_small
    end
    nothing
end

function simple_step_size!(state::PathTrackerState, result::CorrectorResult, options::PathTrackerOptions)
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

    if !isconverged(result) && state.Δs < options.minimal_step_size
        state.status = PathTrackerStatus.terminated_step_size_too_small
    end
end

function check_terminated!(tracker::PathTracker)
    if abs(tracker.state.s - length(tracker.state.segment)) < 2eps(length(tracker.state.segment))
        tracker.state.status = PathTrackerStatus.success
    elseif curriters(tracker) ≥ tracker.options.maxiters
        tracker.state.status = PathTrackerStatus.terminated_maximal_iterations
    end
    nothing
end

function refine!(tracker::PathTracker)
    if tracker.state.accuracy < tracker.options.refinement_tol
        return
    end
    result = correct!(tracker.state.x̄, tracker;
        tol=tracker.options.refinement_tol,
        maxiters=tracker.options.refinement_maxiters)
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


"""
    checkstart(H, x)

Check whether the `x` has the correct size.
"""
function checkstart(H, x)
    N = nvariables(H)
    N != length(x) && throw(error("Expected `x` to have length $(N) but `x` has length $(length(x))"))
    nothing
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
     tol(tracker::PathTracker)

Current tolerance.
"""
tol(tracker::PathTracker) = tracker.options.tol

"""
     set_tol!(tracker::PathTracker, tol)

Set the current tolerance to `tol`.
"""
set_tol!(tracker::PathTracker, tol) = tracker.options.tol = tol

"""
     refinement_tol(tracker::PathTracker)

Current refinement tolerance.
"""
refinement_tol(tracker::PathTracker) = tracker.options.refinement_tol

"""
     set_refinement_maxiters!(tracker::PathTracker, tol)

Set the current refinement tolerance to `tol`.
"""
set_refinement_tol!(T::PathTracker, tol) = T.options.refinement_tol = tol

"""
     refinement_maxiters(tracker::PathTracker)

Current refinement maxiters.
"""
refinement_maxiters(T::PathTracker) = T.options.refinement_maxiters

"""
     set_refinement_maxiters!(tracker::PathTracker, n)

Set the current refinement maxiters to `n`.
"""
set_refinement_maxiters!(T::PathTracker, n) = T.options.refinement_maxiters = n

"""
     corrector_maxiters(tracker::PathTracker)

Current correction maxiters.
"""
corrector_maxiters(T::PathTracker) = T.options.corrector_maxiters

"""
     set_corrector_maxiters!(tracker::PathTracker, n)

Set the correction maxiters to `n`.
"""
set_corrector_maxiters!(T::PathTracker, n) = T.options.corrector_maxiters = n

"""
     maximal_step_size(tracker::PathTracker)

Current maximal step size.
"""
maximal_step_size(T::PathTracker) = T.options.maximal_step_size

"""
     set_corrector_maxiters!(tracker::PathTracker, Δs)

Set the maximal step size to `Δs`.
"""
set_maximal_step_size!(T::PathTracker, Δs) = T.options.maximal_step_size = Δs

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
tracker = pathtracker(f, parameters=p, p₁=a, p₀=b, maximal_step_size=0.01)
for (x, t) in iterator(tracker, x₁)
    @show (x,t)
end
```

If we want to guarantee smooth traces we can limit the maximal step size.
```julia
tracker = pathtracker(f, parameters=p, p₁=a, p₀=b, maximal_step_size=0.01)
for (x, t) in iterator(tracker, x₁)
    @show (x,t)
end
```
"""
function pathtracker(args...; kwargs...)
    tracker, _ = pathtracker_startsolutions(args...; kwargs...)
    tracker
end
