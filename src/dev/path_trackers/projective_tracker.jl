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

mutable struct ProjectiveState{T<:Number, S<:AbstractStepSizeState} <: AbstractPathTrackerState
    start::Float64
    target::Float64
    t::Float64
    steplength::Float64
    steplength_state::S
    iter::Int
    status::Symbol

    x::Vector{T}
    xnext::Vector{T}
    patch::Vector{T} # This vector will also be used in the cache.
    #randomization::Matrix{T}
end

get_t(state::ProjectiveState) = (1-state.t) * state.target + state.t * state.start
get_stepsize(state::ProjectiveState) = state.steplength * (state.start - state.target)

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
        status = :default
    end

    # We have to make sure that the element type of x is invariant under evaluation
    x = convert.(eltype(value), x)
    AffinePatches.precondition!(x, tracker.patch)
    xnext = copy(x)
    patch = AffinePatches.init_patch(tracker.patch, x)

    ProjectiveState(start, target, t, steplength, steplength_state, iter, status, x, xnext, patch)
end

function reset!(state::ProjectiveState, tracker::ProjectiveTracker, x, start, target)
    @assert(NewHomotopies.nvariables(tracker.homotopy) == length(x),
        "Expected `x` to have length $(NewHomotopies.nvariables(tracker.homotopy)) but `x` has length $(length(x))")

    state.t = 1.0
    StepSize.reset!(state.steplength_state)
    state.iter = 0

    if norm(NewHomotopies.evaluate(tracker.homotopy, x, start)) > 1e-4
        status = :invalid_startvalue
    else
        status = :default
    end

    state.x .= x
    AffinePatches.precondition!(state.x, tracker.patch)
    state.xnext .= x
    AffinePatches.update_patch!(state.patch, tracker.patch, x)
end

struct ProjectiveTrackerCache{M, N, H<:PatchedHomotopy{M, N}, HC<:PatchedHomotopyCache, P, C} <: AbstractPathTrackerCache
    homotopy::HomotopyWithCache{M, N, H, HC}
    predictor_corrector::PredictorCorrectorCache{P, C}
end

function cache(tracker::ProjectiveTracker, state::ProjectiveState)
    PH = PatchedHomotopy(tracker.homotopy, state.patch)
    H = NewHomotopies.HomotopyWithCache(PH, state.x, get_t(state))
    pc_cache = cache(tracker.predictor_corrector, H, state.x, get_t(state))

    ProjectiveTrackerCache(H, pc_cache)
end
