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
    steplength::S
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
    # Our start point in space
    start::Complex{Float64}
    # Our target point in space
    target::Complex{Float64}
    t::Float64 # The relative progress. `t` always goes from 1.0 to 0.0
    Δt::Float64 # Δt is the current relative step width
    steplength::S
    iter::Int
    status::Symbol

    x::Vector{T}
    patch::Vector{T} # This vector will also be used in the cache.
    #randomization::Matrix{T}
end

current_t(state::ProjectiveState) = (1-state.t) * state.target + state.t * state.start
current_Δt(state::ProjectiveState) = state.Δt * (state.target - state.start)
current_x(state::ProjectiveState) = state.x
current_iter(state::ProjectiveState) = state.iter
current_status(state::ProjectiveState) = state.status

struct ProjectiveCache{H<:HomotopyWithCache, PC <:PredictorCorrectorCache} <: AbstractPathTrackerCache
    homotopy::H
    predictor_corrector::PC
end


# INITIAL CONSTRUCTION

function state(method::Projective, x::Vector, t₁, t₀)
    checkstart(method.homotopy, x)

    start = convert(Complex{Float64}, t₁)
    target = convert(Complex{Float64}, t₀)
    steplength = StepLength.state(method.steplength, start, target)
    t = 1.0
    Δt = min(StepLength.relsteplength(steplength), 1.0)
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

    ProjectiveState(start, target, t, Δt, steplength, iter, status, x, patch)
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

    state.start = start
    state.target = target
    StepLength.reset!(state.steplength, method.steplength, start, target)
    state.t = 1.0
    state.Δt = min(StepLength.relsteplength(state.steplength), 1.0)

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

    update_t_and_steplength!(state, method.steplength, step_successfull)

    if step_successfull
        # state.t = max(state.t - state.steplength, 0.0)
        # Since x changed we should update the patch
        AffinePatches.update_patch!(state.patch, method.patch, state.x)
    end
    nothing
end

function update_t_and_steplength!(state::ProjectiveState, step::StepLength.AbstractStepLength, step_successfull)
    if step_successfull
        state.t -= state.Δt
    end

    status = StepLength.update!(state.steplength, step, step_successfull)
    if status != :ok
        state.status = status
        return nothing
    end

    relsteplength = StepLength.relsteplength(state.steplength)
    state.Δt = min(relsteplength, state.t)
    # Due to numerical errors we would sometimes be at t=0.1+2eps() and we would just
    # due a step of size 0.1 and then again one of size 2eps() which is quite wasteful.
    if relsteplength + 5e-16 > state.t
        state.Δt = state.t
    else
        state.Δt = relsteplength
    end
end

function check_terminated!(state::ProjectiveState, method::Projective, cache::ProjectiveCache, options::Options)
    if state.t == 0.0
        state.status = :success
    elseif state.iter ≥ options.maxiters
        state.status = :maxiters
    end
    nothing
end

function refine!(state::ProjectiveState, method::Projective, cache::ProjectiveCache, options::Options)
    H, x, t = cache.homotopy, state.x, current_t(state)
    PredictionCorrection.refine!(x,
        method.predictor_corrector, cache.predictor_corrector,
        H, t, options.refinement_tol, options.refinement_maxiters)
end

normalize_x!(x, method::Projective) = normalize!(x)
