module PathTracking

import ProjectiveVectors
import Random
import LinearAlgebra, TreeViews
import ..AffinePatches, ..Correctors, ..Homotopies,
       ..Predictors, ..Problems
using ..Utilities

export PathTracker, PathTrackerResult, pathtracker, pathtracker_startsolutions,
        allowed_keywords, track, track!, setup!, iterator!,
        currx, currt, currΔt, curriters, currstatus, tol, corrector_maxiters,
        refinement_tol, refinement_maxiters, set_tol!,
        set_corrector_maxiters!, set_refinement_tol!, set_refinement_maxiters!

const allowed_keywords = [:corrector, :predictor, :steplength,
    :tol, :refinement_tol, :corrector_maxiters,  :refinement_maxiters,
    :maxiters, :simple_step_size]

###########
# Options
##########
mutable struct Options
    tol::Float64
    corrector_maxiters::Int
    refinement_tol::Float64
    refinement_maxiters::Int
    maxiters::Int
    initial_steplength::Float64
    minimal_steplength::Float64
    simple_step_size::Bool
end

function Options(;tol=1e-7,
    refinement_tol=1e-8,
    corrector_maxiters::Int=2,
    refinement_maxiters=corrector_maxiters,
    maxiters=10_000,
    initial_steplength=0.1,
    minimal_steplength=1e-14,
    simple_step_size=false)

    Options(tol, corrector_maxiters, refinement_tol, refinement_maxiters, maxiters,
            initial_steplength, minimal_steplength, simple_step_size)
end
Base.show(io::IO, opts::Options) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::Options) = opts

###########
# State
##########
@moduleenum Status begin
    success
    tracking
    terminated_maximal_iterations
    terminated_invalid_startvalue
    terminated_steplength_too_small
    terminated_singularity
end

mutable struct State{T, N, PatchState <: AffinePatches.AbstractAffinePatchState}
    x::ProjectiveVectors.PVector{T,N} # current x
    x̂::ProjectiveVectors.PVector{T,N} # last prediction
    x̄::ProjectiveVectors.PVector{T,N} # canidate for new x
    ẋ::Vector{T} # derivative at current x
    η::Float64
    ω::Float64
    segment::ComplexSegment
    s::Float64 # current step length (0 ≤ s ≤ length(segment))
    Δs::Float64 # current step size
    Δs_prev::Float64 # previous step size
    accuracy::Float64
    status::Status.t
    patch::PatchState
    accepted_steps::Int
    rejected_steps::Int
    last_step_failed::Bool
    consecutive_successfull_steps::Int
end

function State(x₁::ProjectiveVectors.PVector, t₁, t₀, patch::AffinePatches.AbstractAffinePatchState, options::Options)
    x, x̂, x̄ = copy(x₁), copy(x₁), copy(x₁)
    ẋ = copy(x₁.data)
    η = ω = NaN
    segment = ComplexSegment(promote(complex(t₁), complex(t₀))...)
    s = 0.0
    Δs = convert(Float64, min(options.initial_steplength, length(segment)))
    Δs_prev = 0.0
    accuracy = 0.0
    accepted_steps = rejected_steps = 0
    status = Status.tracking
    last_step_failed = false
    consecutive_successfull_steps = 0
    State(x, x̂, x̄, ẋ, η, ω, segment, s, Δs, Δs_prev, accuracy, status, patch,
        accepted_steps, rejected_steps, last_step_failed, consecutive_successfull_steps)
end

function reset!(state::State, x₁::AbstractVector, t₁, t₀, options::Options, setup_patch)
    state.segment = ComplexSegment(promote(t₁, t₀)...)
    state.η = state.ω = NaN
    state.s = 0.0
    state.Δs = min(options.initial_steplength, length(state.segment))
    state.Δs_prev = 0.0
    state.accuracy = 0.0
    state.accepted_steps = state.rejected_steps = 0
    state.status = Status.tracking
    ProjectiveVectors.embed!(state.x, x₁)
    setup_patch && AffinePatches.setup!(state.patch, state.x)
    state.last_step_failed = false
    state.consecutive_successfull_steps = 0
    state
end

Base.show(io::IO, state::State) = print_fieldnames(io, state)
Base.show(io::IO, ::MIME"application/prs.juno.inline", state::State) = state


###########
# Cache
##########
mutable struct Cache{H<:Homotopies.HomotopyWithCache, P<:Predictors.AbstractPredictorCache,
             C<:Correctors.AbstractCorrectorCache, T, F}
    homotopy::H
    predictor::P
    corrector::C
    Jac::Jacobian{T, F}
    out::Vector{T}
    r::Vector{T}
end
function Cache(H::Homotopies.HomotopyWithCache, predictor, corrector, state::State)
    t = state.segment[state.s]
    pcache = Predictors.cache(predictor, H, state.x, state.ẋ, t)
    ccache = Correctors.cache(corrector, H, state.x, t)
    Jac = Jacobian(Homotopies.jacobian(H, state.x, t))
    out = H(state.x, t)
    r = copy(out)
    Cache(H, pcache, ccache, Jac, out, r)
end


##############
# PathTracker
##############
"""
     PathTracker(H::Homotopies.AbstractHomotopy, x₁, t₁, t₀; options...)::PathTracker

Create a `PathTracker` to track `x₁` from `t₁` to `t₀`. The homotopy `H`
needs to be homogenous. Note that a `PathTracker` is also a (mutable) iterator.

## Options
* `corrector::Correctors.AbstractCorrector`: The corrector used during in the predictor-corrector scheme. The default is [`Correctors.Newton`](@ref).
* `corrector_maxiters=3`: The maximal number of correction steps in a single step.
* `initial_steplength=0.1`: The step length of the first step.
* `maxiters=10_000`: The maximal number of iterations the path tracker has available.
* `minimal_steplength=1e-14`: The minimal step length.
* `predictor::Predictors.AbstractPredictor`: The predictor used during in the predictor-corrector scheme. The default is `[Predictors.RK4`](@ref)()`.
* `refinement_maxiters=corrector_maxiters`: The maximal number of correction steps used to refine the final value.
* `refinement_tol=1e-8`: The precision used to refine the final value.
* `tol=1e-7`: The precision used to track a value.
"""
struct PathTracker{H<:Homotopies.AbstractHomotopy,
    Predictor<:Predictors.AbstractPredictor,
    Corrector<:Correctors.AbstractCorrector,
    Patch<:AffinePatches.AbstractAffinePatch,
    S<:State,
    C<:Cache}
    # these are fixed
    homotopy::H
    predictor::Predictor
    corrector::Corrector
    affine_patch::Patch
    # these are mutable
    state::S
    options::Options
    cache::C
end

"""
    PathTracker(problem::Problems.AbstractProblem, x₁, t₁, t₀; kwargs...)

Construct a [`PathTracking.PathTracker`](@ref) from the given `problem`.
"""
function PathTracker(prob::Problems.Projective, x₁, t₁, t₀; kwargs...)
    y₁ = Problems.embed(prob, x₁)
    PathTracker(prob.homotopy, y₁, complex(t₁), complex(t₀); kwargs...)
end
function PathTracker(H::Homotopies.AbstractHomotopy, x₁::ProjectiveVectors.PVector, t₁, t₀;
    patch=AffinePatches.OrthogonalPatch(),
    corrector::Correctors.AbstractCorrector=Correctors.Newton(),
    predictor::Predictors.AbstractPredictor=Predictors.Heun(), kwargs...)

    options = Options(;kwargs...)

    if H isa Homotopies.PatchedHomotopy
        error("You cannot pass a `PatchedHomotopy` to PathTracker. Instead pass the homotopy and patch separate.")
    end

    patch_state = AffinePatches.state(patch, x₁)
    # We close over the patch state, the homotopy and its cache
    # to be able to pass things around more easily
    HC = Homotopies.HomotopyWithCache(Homotopies.PatchedHomotopy(H, patch_state), x₁, t₁)

    # We have to make sure that the element type of x is invariant under evaluation
    indempotent_x = begin
        u = Vector{Any}(undef, size(H)[1])
        Homotopies.evaluate!(u, HC, x₁, t₁)
        indem_x = similar(x₁, promote_type(typeof(u[1]), ComplexF64))
        indem_x .= x₁
    end
    state = State(indempotent_x, t₁, t₀, patch_state, options)
    cache = Cache(HC, predictor, corrector, state)

    PathTracker(H, predictor, corrector, patch, state, options, cache)
end

Base.show(io::IO, ::PathTracker) = print(io, "PathTracker()")
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::PathTracker) = x


"""
     PathTrackerResult(tracker)

Containing the result of a tracked path. The fields are
* `successfull::Bool` Indicating whether tracking was successfull.
* `returncode::PathTracking.Status.t` If the tracking was successfull then it is `PathTracking.Status.success`.
* `x::V` The result.
* `t::Float64` The `t` when the path tracker stopped.
"""
struct PathTrackerResult{T, N}
     returncode::Status.t
     x::ProjectiveVectors.PVector{T,N}
     t::ComplexF64
     accuracy::Float64
     accepted_steps::Int
     rejected_steps::Int
end

function PathTrackerResult(tracker::PathTracker)
    state = tracker.state
     PathTrackerResult(state.status,
          copy(state.x), state.segment[state.s],
          state.accuracy,
          state.accepted_steps,
          state.rejected_steps)
end

Base.show(io::IO, result::PathTrackerResult) = print_fieldnames(io, result)
Base.show(io::IO, ::MIME"application/prs.juno.inline", result::PathTrackerResult) = result

"""
    track(tracker, x₁, t₁, t₀; options...)::PathTrackerResult

Track a value `x₁` from `t₁` to `t₀` using the given `PathTracker` `tracker`.
This returns a `PathTrackerResult`. This modifies `tracker`.
See [`track!`](@ref) for the possible options.
"""
function track(tracker::PathTracker, x₁::AbstractVector, t₁, t₀; kwargs...)
     track!(tracker, x₁, t₁, t₀; kwargs...)
     PathTrackerResult(tracker)
end

"""
     track!(tracker, x₁, t₁, t₀; setup_patch=true, checkstartvalue=true, compute_ẋ=true)

Track a value `x₁` from `t₁` to `t₀` using the given `PathTracker` `tracker`.
Returns one of the enum values of `PathTracking.Status.t` indicating the status.
If the tracking was successfull it is `PathTracking.Status.success`.
If `setup_patch` is `true` then [`AffinePatches.setup!`](@ref) is called at the beginning
of the tracking.

    track!(x₀, tracker, x₁, t₁, t₀; options...)

Additionally also stores the result in `x₀` if the tracking was successfull.
"""
function track!(x₀, tracker::PathTracker, x₁, t₁, t₀; setup_patch=true, checkstartvalue=true, compute_ẋ=true)
     track!(tracker, x₁, t₁, t₀, setup_patch, checkstartvalue, compute_ẋ)
     retcode = currstatus(tracker)
     if retcode == Status.success
         x₀ .= currx(tracker)
     end
     retcode
end
function track!(tracker::PathTracker, x₁, t₁, t₀; setup_patch=true, checkstartvalue=true, compute_ẋ=true)
    track!(tracker, x₁, t₁, t₀, setup_patch, checkstartvalue, compute_ẋ)
end
function track!(tracker::PathTracker, x₁, t₁, t₀, setup_patch, checkstartvalue=true, compute_ẋ=true)
    setup!(tracker, x₁, t₁, t₀, setup_patch, checkstartvalue, compute_ẋ)

    while tracker.state.status == Status.tracking
        step!(tracker)
        check_terminated!(tracker)
    end
    if tracker.state.status == Status.success
        refine!(tracker)
    end

    tracker.state.status
end

"""
    setup!(pathtracker, x₁, t₁, t₀, setup_patch=true, checkstartvalue=true, compute_ẋ=true)

Setup `pathtracker` to track `x₁` from `t₁` to `t₀`. Use this if you want to use the
pathtracker as an iterator.
"""
function setup!(tracker::PathTracker, x₁::AbstractVector, t₁, t₀, setup_patch=true, checkstartvalue=true, compute_ẋ=true)
    state, cache = tracker.state, tracker.cache

    try
        reset!(state, x₁, t₁, t₀, tracker.options, setup_patch)
        Predictors.reset!(cache.predictor, state.x, t₁)
        checkstartvalue && checkstartvalue!(tracker)
        if compute_ẋ
            compute_ẋ!(state.ẋ, state.x, currt(state), cache)
            Predictors.setup!(cache.predictor, cache.homotopy, state.x, state.ẋ, currt(state), cache.Jac)
        end
    catch err
        if !(err isa LinearAlgebra.SingularException)
            rethrow(err)
        end
        tracker.state.status = Status.terminated_singularity
    end
    tracker
end

function checkstartvalue!(tracker)
    result = correct!(tracker.state.x̄, tracker)
    if Correctors.isconverged(result)
        tracker.state.x .= tracker.state.x̄
    else
        tracker.state.status = Status.terminated_invalid_startvalue
    end
    nothing
end

function compute_ẋ!(ẋ, x, t, cache)
    @inbounds Homotopies.jacobian_and_dt!(cache.Jac.J, cache.out, cache.homotopy, x, t)
    # applay row scaling to J and compute factorization
    Utilities.updated_jacobian!(cache.Jac)

    @inbounds for i in eachindex(cache.out)
        cache.out[i] = -cache.out[i]
    end

    Utilities.solve!(ẋ, cache.Jac, cache.out)
end

function correct!(x̄, tracker, x=tracker.state.x, t=tracker.state.segment[tracker.state.s];
    tol=tracker.options.tol,
    maxiters=tracker.options.corrector_maxiters)
    Correctors.correct!(x̄, tracker.corrector, tracker.cache.corrector,
                        tracker.cache.homotopy, x, t, tol=tol, maxiters=maxiters)
end

function step!(tracker)
    state, cache, options = tracker.state, tracker.cache, tracker.options
    H = cache.homotopy
    x, x̂, x̄, ẋ = state.x, state.x̂, state.x̄, state.ẋ

    try
        t, Δt = currt(state), currΔt(state)
        Predictors.predict!(x̂, tracker.predictor, cache.predictor, H, x, t, Δt, ẋ)
        result = Correctors.correct!(x̄, tracker.corrector, cache.corrector, H, x̂, t + Δt, tol=options.tol, maxiters=options.corrector_maxiters)
        if Correctors.isconverged(result)
            # Step is accepted, assign values
            state.accepted_steps += 1
            x .= x̄
            state.s += state.Δs
            state.accuracy = result.accuracy

            # Step size change
            update_stepsize!(state, result, Predictors.order(tracker.predictor), options)
            AffinePatches.changepatch!(state.patch, x)
            # update derivative
            compute_ẋ!(ẋ, x, t + Δt, cache)
            # tell the predictors about the new derivative if they need to update something
            Predictors.update!(cache.predictor, H, x, ẋ, t + Δt, cache.Jac)
        else
            # We have to reset the patch
            state.rejected_steps += 1

            # Step failed, so we have to try with a new (smaller) step size
            update_stepsize!(state, result, Predictors.order(tracker.predictor), options)
            Δt = currΔt(state)
            state.last_step_failed = true
            if state.Δs < options.minimal_steplength
                state.status = Status.terminated_steplength_too_small
            end
        end
    catch err
        if !(err isa LinearAlgebra.SingularException)
            rethrow(err)
        end
        tracker.state.status = Status.terminated_singularity
    end
    nothing
end

g(Θ) = sqrt(1+4Θ) - 1
# Choose 0.25 instead of 1.0 due to Newton-Kantorovich theorem
δ(opts::Options, ω) = @fastmath min(√(ω/2) * τ(opts), 0.25)
τ(opts::Options) = nthroot(opts.tol, 2 * opts.corrector_maxiters)

function update_stepsize!(state::State, result::Correctors.Result,
                          order::Int, options::Options)

    if options.simple_step_size
        simple_step_size!(state, result, options)
        return nothing
    end

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
    if Correctors.isconverged(result)
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

    state.Δs = min(Δs′, length(state.segment) - state.s)

    if !Correctors.isconverged(result) && state.Δs < options.minimal_steplength
        state.status = Status.terminated_steplength_too_small
    end
    nothing
end

function simple_step_size!(state::State, result::Correctors.Result, options::Options)
    if Correctors.isconverged(result)
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

    if !Correctors.isconverged(result) && state.Δs < options.minimal_steplength
        state.status = Status.terminated_steplength_too_small
    end
end

function check_terminated!(tracker)
    if abs(tracker.state.s - length(tracker.state.segment)) < 1e-15
        tracker.state.status = Status.success
    elseif curriters(tracker) ≥ tracker.options.maxiters
        tracker.state.status = Status.terminated_maximal_iterations
    end
    nothing
end

function refine!(tracker)
    if tracker.state.accuracy < tracker.options.refinement_tol
        return
    end
    result = correct!(tracker.state.x̄, tracker;
        tol=tracker.options.refinement_tol,
        maxiters=tracker.options.refinement_maxiters)
    if Correctors.isconverged(result)
        tracker.state.x .= tracker.state.x̄
        tracker.state.accuracy = result.accuracy
    end
    nothing
end

# TODO: REMOVE THIS
function residual(tracker, x, t)
    Homotopies.evaluate!(tracker.cache.out, tracker.cache.homotopy, x, t)
    infinity_norm(tracker.cache.out)
end


"""
    checkstart(H, x)

Check whether the `x` has the correct size.
"""
function checkstart(H, x)
    N = Homotopies.nvariables(H)
    N != length(x) && throw(error("Expected `x` to have length $(N) but `x` has length $(length(x))"))
    nothing
end

"""
     currt(tracker::PathTracker)

Current `t`.
"""
currt(tracker::PathTracker) = currt(tracker.state)
currt(state::State) = state.segment[state.s]

"""
     currΔt(tracker::PathTracker)

Current steplength `Δt`.
"""
currΔt(tracker::PathTracker) = currΔt(tracker.state)
currΔt(state::State) = state.segment[state.Δs] - state.segment.start

"""
     curriters(tracker::PathTracker)

Current number of iterations.
"""
curriters(tracker::PathTracker) = curriters(tracker.state)
curriters(state::State) = state.accepted_steps + state.rejected_steps

"""
     currstatus(tracker::PathTracker)

Current status.
"""
currstatus(tracker::PathTracker) = currstatus(tracker.state)
currstatus(state::State) = state.status

"""
    currx(tracker::PathTracker)

Return the current value of `x`.
"""
currx(tracker::PathTracker) = currx(tracker.state)
currx(state::State) = state.x

"""
     tol(tracker::PathTracker)

Current tolerance.
"""
tol(tracker::PathTracker) = tracker.options.tol

"""
     set_tol!(tracker::PathTracker, tol)

Set the current tolerance to `tol`.
"""
function set_tol!(tracker::PathTracker, tol)
     tracker.options.tol = tol
     tol
end

"""
     refinement_tol(tracker::PathTracker)

Current refinement tolerance.
"""
refinement_tol(tracker::PathTracker) = tracker.options.refinement_tol

"""
     set_refinement_maxiters!(tracker::PathTracker, tol)

Set the current refinement tolerance to `tol`.
"""
function set_refinement_tol!(tracker::PathTracker, tol)
     tracker.options.refinement_tol = tol
     tol
end

"""
     refinement_maxiters(tracker::PathTracker)

Current refinement maxiters.
"""
refinement_maxiters(tracker::PathTracker) = tracker.options.refinement_maxiters

"""
     set_refinement_maxiters!(tracker::PathTracker, n)

Set the current refinement maxiters to `n`.
"""
function set_refinement_maxiters!(tracker::PathTracker, n)
     tracker.options.refinement_maxiters = n
     n
end

"""
     corrector_maxiters(tracker::PathTracker)

Current correction maxiters.
"""
corrector_maxiters(tracker::PathTracker) = tracker.options.corrector_maxiters

"""
     set_corrector_maxiters!(tracker::PathTracker, n)

Set the current correction maxiters to `n`.
"""
function set_corrector_maxiters!(tracker::PathTracker, n)
     tracker.options.corrector_maxiters = n
     n
end


"""
    pathtracker_startsolutions(args...; kwargs...)

Construct a [`PathTracking.PathTracker`](@ref) and `startsolutions` in the same way `solve`
does it. This also takes the same input arguments as `solve`. This is convenient if you want
to investigate single paths.
"""
function pathtracker_startsolutions(args...; kwargs...)
    supported, rest = splitkwargs(kwargs, Problems.supported_keywords)
    prob, startsolutions = Problems.problem_startsolutions(args...; supported...)
    tracker = PathTracker(prob, Utilities.start_solution_sample(startsolutions), one(ComplexF64), zero(ComplexF64); rest...)

    (tracker=tracker, startsolutions=startsolutions)
end

"""
    pathtracker(args...; kwargs...)

Construct a [`PathTracking.PathTracker`](@ref) in the same way `solve`
does it. This als0 takes the same input arguments as `solve`. This is convenient if you want
to investigate single paths.
"""
function pathtracker(args...; kwargs...)
    tracker, _ = pathtracker_startsolutions(args...; kwargs...)
    tracker
end

"""
    iterator!(tracker::PathTracker, x₁, t₁, t₀)

Prepare a tracker to make it usable as a (stateful) iterator. Use this if you want to inspect
the state of the pathtracker at each iteration. In each iteration the PathTracker is returned
from the iterator.

## Example

Assume you have `PathTracker` `pathtracker` and you wan to track `x₁` from 1.0 to 0.25:
```julia
for tracker in iterator!(pathtracker, x₁, 1.0, 0.25)
    println("Current t: \$(real(PathTracking.currt(tracker)))") # The time `t` is always a complex number
end
```

Note that if you want to store the current value of `x` you have to create a **copy**.
`x` will also be a projective vector.
```julia
xs = []
for tracker in iterator!(pathtracker, x₁, 1.0, 0.25)
    x = PathTracking.currx(tracker)
     # We want to get the affine vector, this also creates a copy
    push!(xs, ProjectiveVectors.affine_chart(x))
end
```
"""
iterator!(tracker::PathTracker, x₁, t₁, t₀; kwargs...) = setup!(tracker, x₁, t₁, t₀; kwargs...)

function Base.iterate(tracker::PathTracker, state=1)
    if tracker.state.status == Status.tracking
        # return initial tracker once
        state == 1 && return tracker, state + 1
        step!(tracker)
        check_terminated!(tracker)

        if tracker.state.status == Status.success
            refine!(tracker)
        end
        tracker, state + 1
    else
        nothing
    end
end

end
