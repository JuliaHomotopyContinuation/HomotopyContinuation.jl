import StaticArrays: SVector

import ..NewHomotopies
import ..NewHomotopies: AbstractHomotopy, AbstractStartTargetHomotopy, AbstractHomotopyCache, HomotopyWithCache
import ..Predictors
import ..Correctors
import ..PredictionCorrection: PredictorCorrector, PredictorCorrectorCache
import ..PredictionCorrection
import ..AffinePatches
import ..AffinePatches: AbstractAffinePatch
import ..StepLength
import ..StepLength: AbstractStepLength, AbstractStepLengthState
using ..Utilities


export Projective

include("projective/patched_homotopy.jl")

"""
    Projective(H::AbstractHomotopy; options...) <: AbstractPathTracker

Construct a path tracker for a projective homotopy `H`.
The options are
* `patch::AffinePatches.AbstractAffinePatch`:
The affine patch used to embed the homotopy into affine space (from projective space).
The default is [`AffinePatches.OrthogonalPatch`](@ref).
* `predictor::Predictors.AbstractPredictor`:
The predictor used during in the predictor-corrector scheme. The default is
`[Predictors.RK4`](@ref)()`.
* `corrector::Correctors.AbstractCorrector`:
The corrector used during in the predictor-corrector scheme. The default is
[`Correctors.Newton`](@ref).
* `step::StepLength.AbstractStepLength`
The step size logic used to determine changes of the step size. The default is
[`StepLength.HeuristicStepLength`](@ref).
"""
struct Projective{H<:AbstractHomotopy, Patch<:AbstractAffinePatch, P, C, S<:AbstractStepLength} <: AbstractPathTrackerMethod
    homotopy::H
    patch::Patch
    # randomize::Randomizer
    predictor_corrector::PredictorCorrector{P, C}
    step::S
end

function Projective(H::AbstractHomotopy;
    patch::AffinePatches.AbstractAffinePatch=AffinePatches.OrthogonalPatch(),
    predictor::Predictors.AbstractPredictor=Predictors.RK4(),
    corrector::Correctors.AbstractCorrector=Correctors.Newton(),
    step::StepLength.AbstractStepLength=StepLength.HeuristicStepLength())

    Projective(
        H,
        patch,
        PredictorCorrector(predictor, corrector), step)
end

# STATE

mutable struct ProjectiveState{T<:Number, S<:AbstractStepLengthState} <: AbstractPathTrackerState
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
    patch::Vector{T} # This vector will also be used in the cache.
    #randomization::Matrix{T}
end


current_t(state::ProjectiveState) = (1-state.t) * state.target + state.t * state.start
current_Δt(state::ProjectiveState) = state.steplength * (state.target - state.start)
current_x(state::ProjectiveState) = state.x
current_iter(state::ProjectiveState) = state.iter
current_status(state::ProjectiveState) = state.status

struct ProjectiveCache{M, N, H<:PatchedHomotopy{M, N}, HC<:PatchedHomotopyCache, P, C} <: AbstractPathTrackerCache
    homotopy::HomotopyWithCache{M, N, H, HC}
    predictor_corrector::PredictorCorrectorCache{P, C}
end


# INITIAL CONSTRUCTION

function state(method::Projective, x::Vector, start, target)
    checkstart(method.homotopy, x)

    t = 1.0
    steplength = StepLength.initial_steplength(method.step)
    steplength_state = StepLength.state(method.step)
    iter = 0

    value = NewHomotopies.evaluate(method.homotopy, x, start)
    if infinity_norm(value) > 1e-4
        status = :invalid_startvalue
    else
        status = :ok
    end

    # We have to make sure that the element type of x is invariant under evaluation
    x = convert.(eltype(value), x)
    patch = AffinePatches.init_patch(method.patch, x)
    AffinePatches.precondition!(patch, x, method.patch)

    ProjectiveState(start, target, t, steplength, steplength_state, iter, status, x, patch)
end

function checkstart(H, x)
    N = NewHomotopies.nvariables(H)
    N != length(x) && throw(error("Expected `x` to have length $(N) but `x` has length $(length(x))"))
end


function cache(method::Projective, state::ProjectiveState)
    PH = PatchedHomotopy(method.homotopy, state.patch)
    H = NewHomotopies.HomotopyWithCache(PH, state.x, current_t(state))
    pc_cache = PredictionCorrection.cache(method.predictor_corrector, H, state.x, current_t(state))

    ProjectiveCache(H, pc_cache)
end

# UPDATES

function reset!(state::ProjectiveState, method::Projective, cache::ProjectiveCache, x, start, target)
    checkstart(method.homotopy, x)

    state.t = 1.0
    state.steplength = StepLength.initial_steplength(method.step)
    StepLength.reset!(state.steplength_state)
    state.iter = 0
    state.x .= x
    AffinePatches.precondition!(state.patch, state.x, method.patch)

    if infinity_norm(cache.homotopy(state.x, start)) > 1e-4
        state.status = :invalid_startvalue
    else
        state.status = :ok
    end
    state
end

# ITERATION

function step!(state::ProjectiveState, method::Projective, cache::ProjectiveCache, options::Options)
    state.iter += 1

    step_successfull, status = PredictionCorrection.predict_correct!(state.x,
        method.predictor_corrector,
        cache.predictor_corrector,
        cache.homotopy, current_t(state), current_Δt(state), options.tol, options.corrector_maxiters)

    if step_successfull
        state.t = max(state.t - state.steplength, 0.0)
        # Since x changed we should update the patch
        AffinePatches.update_patch!(state.patch, method.patch, state.x)
    end
    #  update steplength
    update_steplength!(state, method.step, step_successfull)
    nothing
end

function update_steplength!(state::ProjectiveState, step::StepLength.AbstractStepLength, step_successfull)
    new_steplength, status = StepLength.update(state.steplength, step, state.steplength_state, step_successfull)
    if status != :ok
        state.status = status
    end
    if StepLength.isrelative(step)
        state.steplength = min(new_steplength, state.t)
    else
        state.steplength = min(new_steplength, state.t) / (state.start - state.target)
    end
    nothing
end

function check_terminated!(state::ProjectiveState, method::Projective, cache::ProjectiveCache, options::Options)
    if state.t == 0.0
        state.status = :success
    elseif state.iter ≥ options.maxiters
        state.status = :maxiters
    end
end

function refine!(state::ProjectiveState, method::Projective, cache::ProjectiveCache, options::Options)
    H, x, t = cache.homotopy, state.x, current_t(state)
    PredictionCorrection.refine!(x,
        method.predictor_corrector, cache.predictor_corrector,
        H, t, options.refinement_tol, options.corrector_maxiters)
end

normalize_x!(x, method::Projective) = normalize!(x)
