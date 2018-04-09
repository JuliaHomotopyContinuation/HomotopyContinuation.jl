export ProjectiveTracker

"""
    ProjectiveTracker(H::AbstractHomotopy; options...) <: AbstractPathTracker

Construct a path tracker for a projective homotopy `H`.
The options are
* `patch::AffinePatches.AbstractAffinePatch`:
The affine patch used to embed the homotopy into affine space (from projective space).
The default is [`AffinePatches.OrthogonalPatch`](@ref).
* `predictor::Predictors.AbstractPredictor`:
The predictor used during in the predictor-corrector scheme. The default is
`[Predictors.Euler`](@ref)()`.
* `corrector::Correctors.AbstractCorrector`:
The corrector used during in the predictor-corrector scheme. The default is
[`Correctors.Newton`](@ref).
* `step::StepSize.AbstractStepSize`
The step size logic used to determine changes of the step size. The default is
[`StepSize.HeuristicStepSize`](@ref).
* `options=Options()`: The options for the path tracker. See [Option](@ref) for details.
"""
struct ProjectiveTracker{H<:AbstractHomotopy, Patch<:AbstractAffinePatch, P, C, S<:AbstractStepSize} <: AbstractPathTracker{H}
    homotopy::H
    patch::Patch
    # randomize::Randomizer
    predictor_corrector::PredictorCorrector{P, C}
    step::S
    options::Options
end

function ProjectiveTracker(H::AbstractHomotopy;
    patch::AffinePatches.AbstractAffinePatch=AffinePatches.OrthogonalPatch(),
    predictor::Predictors.AbstractPredictor=Predictors.Euler(),
    corrector::Correctors.AbstractCorrector=Correctors.Newton(),
    step::StepSize.AbstractStepSize=StepSize.HeuristicStepSize(),
    options::Options=Options())

    ProjectiveTracker(
        H,
        patch,
        PredictorCorrector(predictor, corrector), step, options)
end

# STATE

mutable struct ProjectiveState{T<:Number, S<:AbstractStepSizeState} <: AbstractPathTrackerState
    # Our start time
    start::Float64
    # Our target time
    target::Float64
    # The relative progress. `t` always goes from 1.0 to 0.0
    t::Float64
    # The current step length. This is always relative.
    steplength::Float64
    steplength_state::S
    iter::Int
    status::Symbol

    x::Vector{T}
    xnext::Vector{T}
    patch::Vector{T} # This vector will also be used in the cache.
    #randomization::Matrix{T}
end

struct ProjectiveTrackerCache{M, N, H<:PatchedHomotopy{M, N}, HC<:PatchedHomotopyCache, P, C} <: AbstractPathTrackerCache
    homotopy::HomotopyWithCache{M, N, H, HC}
    predictor_corrector::PredictorCorrectorCache{P, C}
end

function state(tracker::ProjectiveTracker, x::Vector, start, target)
    @assert(NewHomotopies.nvariables(tracker.homotopy) == length(x),
        "Expected `x` to have length $(NewHomotopies.nvariables(tracker.homotopy)) but `x` has length $(length(x))")

    t = 1.0
    steplength = StepSize.initial_steplength(tracker.step)
    steplength_state = StepSize.state(tracker.step)
    iter = 0

    value = NewHomotopies.evaluate(tracker.homotopy, x, start)
    if norm(value) > 1e-4
        status = :invalid_startvalue
    else
        status = :ok
    end

    # We have to make sure that the element type of x is invariant under evaluation
    x = convert.(eltype(value), x)
    AffinePatches.precondition!(x, tracker.patch)
    xnext = copy(x)
    patch = AffinePatches.init_patch(tracker.patch, x)

    ProjectiveState(start, target, t, steplength, steplength_state, iter, status, x, xnext, patch)
end

function cache(tracker::ProjectiveTracker, state::ProjectiveState)
    PH = PatchedHomotopy(tracker.homotopy, state.patch)
    H = NewHomotopies.HomotopyWithCache(PH, state.x, get_t(state))
    pc_cache = cache(tracker.predictor_corrector, H, state.x, get_t(state))

    ProjectiveTrackerCache(H, pc_cache)
end

function reset!(state::ProjectiveState, tracker::ProjectiveTracker, cache::ProjectiveTrackerCache, x, start, target)
    NewHomotopies.nvariables(tracker.homotopy) != length(x) && throw(error("Expected `x` to have length $(NewHomotopies.nvariables(tracker.homotopy)) but `x` has length $(length(x))"))

    state.t = 1.0
    state.steplength = StepSize.initial_steplength(tracker.step)
    StepSize.reset!(state.steplength_state)
    state.iter = 0
    state.x .= x
    AffinePatches.precondition!(state.x, tracker.patch)
    state.xnext .= state.x
    AffinePatches.update_patch!(state.patch, tracker.patch, state.x,)

    if norm(cache.homotopy(state.x, start)) > 1e-4
        state.status = :invalid_startvalue
    else
        state.status = :ok
    end
    state
end


# ITERATOR
get_t(state::ProjectiveState) = (1-state.t) * state.target + state.t * state.start

function update_stepsize!(state::ProjectiveState, tracker::ProjectiveTracker, step_successfull)
    new_steplength, status = StepSize.update(state.steplength, tracker.step, state.steplength_state, step_successfull)
    if status != :ok
        tracker.status = status
    end
    if StepSize.isrelative(tracker.step)
        state.steplength = min(new_steplength, state.t)
    else
        state.steplength = min(new_steplength, state.t) / (state.start - state.target)
    end
    nothing
end

function step!(tracker::ProjectiveTracker, state::ProjectiveState, cache::ProjectiveTrackerCache)
    state.iter += 1

    H = cache.homotopy
    t = get_t(state)
    Δt = state.steplength * (state.start - state.target)

    x = state.x


    # try
    step_successfull = step!(state.xnext,
        tracker.predictor_corrector,
        cache.predictor_corrector,
        H, x, t, Δt, tracker.options.tolerance)
    # # catch err
    # #     warn(err)
    # #     state.status = :predictor_corrector_failed
    # #     return nothing
    # # end

    if step_successfull
        x .= state.xnext
        state.t = max(state.t - state.steplength, 0.0)
    end

    if step_successfull
        # Since x changed we should update the patch
        AffinePatches.update_patch!(state.patch, tracker.patch, x)
    end
    #  update steplength
    update_stepsize!(state, tracker, step_successfull)
    nothing
end

function isdone(tracker::ProjectiveTracker, state::ProjectiveState, cache::ProjectiveTrackerCache)
    state.t == 0.0 || state.status != :ok || state.iter ≥ tracker.options.maxiters
end
