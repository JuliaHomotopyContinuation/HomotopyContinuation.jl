module PathTracking

import LinearAlgebra, TreeViews
import ..AffinePatches, ..Correctors, ..Homotopies,
       ..Predictors, ..Problems, ..ProjectiveVectors
using ..Utilities

export PathTracker, allowed_keywords

const allowed_keywords = [:corrector, :predictor, :steplength,
    :tol, :refinement_tol, :corrector_maxiters,  :refinement_maxiters,
    :maxiters]

########################
# Constructors and show
########################

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
end

function Options(;tol=1e-7,
    refinement_tol=1e-8,
    corrector_maxiters::Int=3,
    refinement_maxiters=corrector_maxiters,
    maxiters=500,
    initial_steplength=0.1,
    minimal_steplength=1e-14)

    Options(tol, corrector_maxiters, refinement_tol, refinement_maxiters, maxiters,
            initial_steplength, minimal_steplength)
end
Base.show(io::IO, opts::Options) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::Options) = opts

τ(opts::Options) = nthroot(opts.tol, 2 * (opts.corrector_maxiters - 1))
τ(opts::Options, λ) = nthroot(λ * opts.tol, 2 * (opts.corrector_maxiters - 1))
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

mutable struct State{T, PV<:ProjectiveVectors.AbstractProjectiveVector{T}, PatchState <: AffinePatches.AbstractAffinePatchState}
    x::PV # current x
    x̂::PV # last prediction
    x̄::PV # canidate for new x
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
end

function State(x₁::ProjectiveVectors.AbstractProjectiveVector, t₁, t₀, patch::AffinePatches.AbstractAffinePatchState, options::Options)
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
    State(x, x̂, x̄, ẋ, η, ω, segment, s, Δs, Δs_prev, accuracy, status, patch, accepted_steps, rejected_steps)
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
    J::Matrix{T}
    J_factorization::F
    out::Vector{T}
end
function Cache(H::Homotopies.HomotopyWithCache, predictor, corrector, state::State)
    t = state.segment[state.s]
    pcache = Predictors.cache(predictor, H, state.x, t)
    ccache = Correctors.cache(corrector, H, state.x, t)
    J = Homotopies.jacobian(H, state.x, t)
    J_factorization = factorization(J)
    out = H(state.x, t)
    Cache(H, pcache, ccache, J, J_factorization, out)
end


##############
# PathTracker
##############
"""
     PathTracker(H::Homotopies.AbstractHomotopy, x₁, t₁, t₀; options...)::PathTracker

Create a `PathTracker` to track `x₁` from `t₁` to `t₀`. The homotopy `H`
needs to be homogenous. Note that a `PathTracker` is also a (mutable) iterator.

## Options
* `corrector::Correctors.AbstractCorrector`:
The corrector used during in the predictor-corrector scheme. The default is
[`Correctors.Newton`](@ref).
* `corrector_maxiters=2`: The maximal number of correction steps in a single step.
* `predictor::Predictors.AbstractPredictor`:
The predictor used during in the predictor-corrector scheme. The default is
`[Predictors.RK4`](@ref)()`.
* `refinement_maxiters=corrector_maxiters`: The maximal number of correction steps used to refine the final value.
* `refinement_tol=1e-8`: The precision used to refine the final value.
* `steplength::StepLength.AbstractStepLength`
The step size logic used to determine changes of the step size. The default is
[`StepLength.HeuristicStepLength`](@ref).
* `tol=1e-6`: The precision used to track a value.
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
function PathTracker(H::Homotopies.AbstractHomotopy, x₁::ProjectiveVectors.AbstractProjectiveVector, t₁, t₀;
    patch=AffinePatches.OrthogonalPatch(),
    corrector::Correctors.AbstractCorrector=Correctors.Newton(),
    predictor::Predictors.AbstractPredictor=Predictors.RK4(), kwargs...)

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
        similar(x₁, promote_type(typeof(u[1]), ComplexF64))
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
* `returncode::Symbol` If the tracking was successfull then it is `:success`.
Otherwise the return code gives an indication what happened.
* `x::V` The result.
* `t::Float64` The `t` when the path tracker stopped.
"""
struct PathTrackerResult{V<:AbstractVector, T}
     returncode::Status.t
     x::V
     t::T
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
    track(tracker, x₁, t₁, t₀)::PathTrackerResult

Track a value `x₁` from `t₁` to `t₀` using the given `PathTracker` `tracker`.
This returns a `PathTrackerResult`. This modifies `tracker`.
"""
function track(tracker::PathTracker, x₁::AbstractVector, t₁, t₀; kwargs...)
     track!(tracker, x₁, t₁, t₀; kwargs...)
     PathTrackerResult(tracker)
end

"""
     track!(tracker, x₁, t₁, t₀; checkstartvalue=true, precondition=true)

Track a value `x₁` from `t₁` to `t₀` using the given `PathTracker` `tracker`.
Returns a `Symbol` indicating the status.
If the tracking was successfull it is `:success`.
If `predcondition` is `true` then [`Homotopies.precondition!`](@ref) is called at the beginning
of the tracking.

    track!(x₀, tracker, x₁, t₁, t₀)

Additionally also stores the result in `x₀` if the tracking was successfull.
"""
function track!(x₀, tracker::PathTracker, x₁, t₁, t₀)
     track!(tracker, x₁, t₁, t₀)
     retcode = currstatus(tracker)
     if retcode == Status.success
         x₀ .= currx(tracker)
     end
     retcode
end
function track!(tracker::PathTracker, x₁, t₁, t₀, args...)
    setup!(tracker, x₁, t₁, t₀, args...)

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
        compute_ẋ && compute_ẋ!(state.ẋ, state.x, currt(state), cache)
        Predictors.setup!(cache.predictor, state.x, state.ẋ, currt(state))
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
    @inbounds Homotopies.jacobian_and_dt!(cache.J, cache.out, cache.homotopy, x, t)
    @inbounds for i in eachindex(cache.out)
        cache.out[i] = -cache.out[i]
    end
    cache.J_factorization = Utilities.factorize!(cache.J_factorization, cache.J)
    Utilities.solve!(ẋ, cache.J_factorization, cache.out)
end


function correct!(x̄, tracker, x=tracker.state.x, t=tracker.state.segment[tracker.state.s];
    tol=tracker.options.tol,
    maxiters=tracker.options.corrector_maxiters)
    Correctors.correct!(x̄, tracker.corrector, tracker.cache.corrector,
                        tracker.cache.homotopy, x, t, tol, maxiters)
end

function step!(tracker)
    state, cache, options = tracker.state, tracker.cache, tracker.options
    H = cache.homotopy
    x, x̂, x̄, ẋ = state.x, state.x̂, state.x̄, state.ẋ

    try
        t, Δt = currt(state), currΔt(state)
        Predictors.predict!(x̂, tracker.predictor, cache.predictor, H, x, t, Δt, ẋ)
        result = Correctors.correct!(x̄, tracker.corrector, cache.corrector, H, x̂, t + Δt, options.tol, options.corrector_maxiters)
        if Correctors.isconverged(result)
            AffinePatches.changepatch!(state.patch, x̄)

            # Step is accepted, assign values
            state.accepted_steps += 1
            x .= x̄
            # state.Δs_prev = state.Δs
            state.s += state.Δs
            state.accuracy = result.accuracy

            # update derivative
            compute_ẋ!(ẋ, x̄, t + Δt, cache)
            # tell the predictors about the new derivative if they need to update something
            Predictors.update!(cache.predictor, H, x, ẋ, t + Δt, cache.J_factorization)

            # Step size change
            update_stepsize!(state, result, Predictors.order(tracker.predictor), options)

            return
        else
            # We have to reset the patch
            AffinePatches.changepatch!(state.patch, x)
            state.rejected_steps += 1

            # Step failed, so we have to try with a new (smaller) step size
            update_stepsize!(state, result, Predictors.order(tracker.predictor), options)
            Δt = currΔt(state)
            if state.Δs < options.minimal_steplength
                state.status = Status.terminated_steplength_too_small
                return
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
function update_stepsize!(state::State, result::Correctors.Result,
                          order::Int, options::Options)
    # we have to handle the special case that there is only 1 iteration
    # in this case we cannot estimate ω and therefore just assume ω = 2
    # Also note ||x̂-x|| = ||Δx₀||
    if result.iters == 1
        ω = 2.0
        d_x̂_x̄ = result.norm_Δx₀
    else
        ω = result.ω
        d_x̂_x̄ = euclidean_distance(state.x̂, state.x̄)
    end

    Δx₀ = result.norm_Δx₀
    τN = τ(options, 0.5)
    if Correctors.isconverged(result)
        # compute η and update
        η = d_x̂_x̄ / state.Δs^order
        if isnan(state.ω)
            ω′ = ω
        else
            # assume Δs′ = Δs
            ω′ = max(2ω - state.ω, ω)
        end
        if isnan(state.η)
            λ = g(√(ω′/2) * τN) / (ω′ * d_x̂_x̄)
            Δs′ = nthroot(λ, order) * state.Δs
        else
            d_x̂_x̄′ = max(2d_x̂_x̄ - state.η * state.Δs^(order), 0.25d_x̂_x̄)
            η′ = d_x̂_x̄′ / state.Δs^order
            λ = g(√(ω′/2) * τN) / (ω′ * d_x̂_x̄′)
            Δs′ = nthroot(λ, order) * state.Δs
        end

        Δs′ = max(Δs′, state.Δs)
        state.η = η
        state.ω = ω
    else
        ω = max(ω, state.ω)
        if √(ω / 2) * τ(options, 0.1) < ω / 2 * Δx₀
            λ = g(√(ω / 2) * τ(options, 0.1)) / g(ω / 2 * Δx₀)
        # require a slightly more restrictive τ
        elseif √(ω / 2) * τ(options, 0.01) < ω / 2 * Δx₀
            λ = g(√(ω / 2) * τ(options, 0.01)) / g(ω / 2 * Δx₀)
        elseif √(ω / 2) * τ(options, 0.001) < ω / 2 * Δx₀
            λ = g(√(ω / 2) * τ(options, 0.001)) / g(ω / 2 * Δx₀)
        else
            # fallback
            λ = 0.25
        end
        # To avoid too minimal changes, e.g. λ=0.9999 we require that λ is at most
        λ = min(λ, 0.95)
        Δs′ = nthroot(λ, order) * state.Δs
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
    push!(xs, ProjectiveVectors.affine(x))
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
