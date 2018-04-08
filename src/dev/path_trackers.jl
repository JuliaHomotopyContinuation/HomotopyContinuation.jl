module PathTrackers

import ..NewHomotopies
import ..NewHomotopies: AbstractStartTargetHomotopy, AbstractHomotopyCache, HomotopyWithCache
import ..Predictors
import ..Correctors
import ..AffinePatches
import ..AffinePatches: AbstractAffinePatch
import ..StepSize
import ..StepSize: AbstractStepSize, AbstractStepSizeState

include("path_trackers/patched_homotopy.jl")
include("path_trackers/predictor_corrector.jl")

abstract type AbstractPathTracker end
abstract type AbstractPathTrackerCache end
abstract type AbstractPathTrackerState end


struct Options
    tolerance::Float64
    maxiters::Int64
end
function Options(;tolerance=1e-7, maxiters=10_000)
    Options(tolerance, maxiters)
end



struct Iterator{Tracker<:AbstractPathTracker, C<:AbstractPathTrackerCache, S<:AbstractPathTrackerState}
    tracker::Tracker
    # TODO: The following actually depend on the current precision. So we would need to introduce
    # multiple of these and switch if necessary.
    cache::C
    state::S
end



# Pathtracker implementation
struct ProjectiveTracker{PH<:PatchedHomotopy, P, C, S<:AbstractStepSize} <: AbstractPathTracker
    homotopy::PH
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
        PatchedHomotopy(H, patch),
        PredictorCorrector(predictor, corrector), step, options)
end

mutable struct ProjectiveState{T<:Number, S<:AbstractStepSizeState}
    start::Float64
    target::Float64
    t::Float64
    steplength::Float64
    steplength_state::S
    iter::Int
    status::Symbol

    x::Vector{T}
    xnext::Vector{T}
    patch::Vector{T}
end

get_t(state::ProjectiveState) = (1-state.t) * state.target + state.t * state.start
get_stepsize(state::ProjectiveState) = state.steplength * (state.start - state.target)

function state(tracker::ProjectiveTracker, x::Vector; start=1.0, target=0.0)
    step = tracker.steplength
    t = 1.0
    steplength = StepSize.initial_steplength(step)
    steplength_state = StepSize.state(step)
    iter = 0
    status = :default

    x = copy(x)
    xnext = copy(x)
    affine_patch = AffinePatches.init_patch(patch(tracker.homotopy), x)

    ProjectiveState(start, target, t, steplength, steplength_state, iter, status, x, xnext, affine_patch)
end


struct ProjectiveTrackerCache{H<:AbstractHomotopyCache, P, C} <: AbstractPathTrackerCache
    homotopy::H
    predictor_corrector::PredictorCorrectorCache{P, C}
end

function cache(tracker::ProjectiveTracker, state::ProjectiveState)
    H = NewHomotopies.HomotopyWithCache(tracker.homotopy, state.x, get_t(state))
    homotopy_cache = H.cache
    pc_cache = cache(tracker.predictor_corrector, H, state.x, get_t(stae))

    ProjectiveTrackerCache(homotopy_cache, pc_cache)
end


end
