import ..AffinePatches
import ..Correctors
import ..Homotopies
import ..PredictionCorrection
import ..Predictors
import ..Problems
import ..ProjectiveVectors
import ..StepLength

using ..Utilities

export PathTracker,
    t, Δt,
    status, value,
    tol,
    corrector_maxiters,
    refinement_tol,
    refinement_maxiters,
    set_tol!,
    set_corrector_maxiters!,
    set_refinement_tol!,
    set_refinement_maxiters!,
    fixpatch!



mutable struct Options
    tol::Float64
    corrector_maxiters::Int
    refinement_tol::Float64
    refinement_maxiters::Int
    maxiters::Int
    # we can disable an updating of the patch
    fixedpatch::Bool
end

mutable struct State{S<:StepLength.AbstractStepLengthState}
    # Our start point in space
    start::Complex{Float64}
    # Our target point in space
    target::Complex{Float64}
    steplength::S
    t::Float64 # The relative progress. `t` always goes from 1.0 to 0.0
    Δt::Float64 # Δt is the current relative step width
    iters::Int
    status::Symbol
end

function State(H::Homotopies.AbstractHomotopy, step,
    x₁::ProjectiveVectors.AbstractProjectiveVector,
    t₁, t₀)
    checkstart(H, x₁)

    start = convert(Complex{Float64}, t₁)
    target = convert(Complex{Float64}, t₀)
    steplength = StepLength.state(step, start, target)
    t = 1.0
    Δt = min(StepLength.relsteplength(steplength), 1.0)
    iters = 0

    status = :ok

    State(start, target, steplength, t, Δt, iters, status)
end
function checkstart(H, x)
    N = Homotopies.nvariables(H)
    N != length(x) && throw(error("Expected `x` to have length $(N) but `x` has length $(length(x))"))
end

struct Cache{H<:Homotopies.HomotopyWithCache, PC<:PredictionCorrection.PredictorCorrectorCache, T}
    homotopy::H
    predictor_corrector::PC
    out::Vector{T}
end
function Cache(homotopy, patch, predictor_corrector, state::State, x₁)
    patched = PatchedHomotopy(homotopy, copy(x₁))
    raw_x₁ = ProjectiveVectors.raw(x₁)
    H = Homotopies.HomotopyWithCache(patched, raw_x₁, state.start)
    pc_cache = PredictionCorrection.cache(predictor_corrector, H, raw_x₁, state.start)
    out = H(raw_x₁, state.start)
    Cache(H, pc_cache, out)
end

"""
     PathTracker(H::Homotopies.AbstractHomotopy, x₁, t₁, t₀; options...)::PathTracker

Create a `PathTracker` to track `x₁` from `t₁` to `t₀`. The homotopy `H`
needs to be homogenous.

## Options
* `corrector::Correctors.AbstractCorrector`:
The corrector used during in the predictor-corrector scheme. The default is
[`Correctors.Newton`](@ref).
* `corrector_maxiters=2`: The maximal number of correction steps in a single step.
* `patch::AffinePatches.AbstractAffinePatch`:
The affine patch used to embed the homotopy into affine space (from projective space).
The default is [`AffinePatches.OrthogonalPatch`](@ref).
* `predictor::Predictors.AbstractPredictor`:
The predictor used during in the predictor-corrector scheme. The default is
`[Predictors.RK4`](@ref)()`.
* `refinement_maxiters=corrector_maxiters`: The maximal number of correction steps used to refine the final value.
* `refinement_tol=1e-11`: The precision used to refine the final value.
* `steplength::StepLength.AbstractStepLength`
The step size logic used to determine changes of the step size. The default is
[`StepLength.HeuristicStepLength`](@ref).
* `tol=1e-5`: The precision used to track a value.

     PathTracker(method::AbstractPathTrackerMethod, x₁, t₁, t₀, options)

If a method is already assembled this constructor is beneficial.
"""
struct PathTracker{
    H<:Homotopies.AbstractHomotopy,
    AP<:AffinePatches.AbstractAffinePatch,
    P<:Predictors.AbstractPredictor,
    Corr<:Correctors.AbstractCorrector,
    SL<:StepLength.AbstractStepLength,
    S<:State,
    C<:Cache,
    V<:ProjectiveVectors.AbstractProjectiveVector}

    # these are fixed
    homotopy::H
    patch::AP
    # randomize::Randomizer
    predictor_corrector::PredictionCorrection.PredictorCorrector{P, Corr}
    steplength::SL

    # these are mutable
    state::S
    options::Options

    # TODO: The following actually depend on the current precision. So we would need to introduce
    # multiple of these and switch if necessary to allow multiple precision.
    x::V
    cache::C
end

function PathTracker(prob::Problems.AbstractProblem, x₁::AbstractVector, t₁, t₀; kwargs...)
     PathTracker(prob, Problems.embed(prob, x₁), t₁, t₀; kwargs...)
end
function PathTracker(prob::Problems.AbstractProblem, x₁::ProjectiveVectors.AbstractProjectiveVector, t₁, t₀; kwargs...)
     PathTracker(prob.homotopy, x₁, t₁, t₀; kwargs...)
end
function PathTracker(H::Homotopies.AbstractHomotopy, x₁::ProjectiveVectors.AbstractProjectiveVector, t₁, t₀;
    corrector::Correctors.AbstractCorrector=Correctors.Newton(),
    patch::AffinePatches.AbstractAffinePatch=AffinePatches.OrthogonalPatch(),
    predictor::Predictors.AbstractPredictor=Predictors.RK4(),
    steplength::StepLength.AbstractStepLength=StepLength.HeuristicStepLength(),
    tol=1e-7,
    refinement_tol=1e-11,
    corrector_maxiters::Int=2,
    refinement_maxiters=corrector_maxiters,
    maxiters=10_000,
    fixedpatch=false)

    predictor_corrector = PredictionCorrection.PredictorCorrector(predictor, corrector)
    # We have to make sure that the element type of x is invariant under evaluation
    x = ProjectiveVectors.converteltype(x₁, eltype(Homotopies.evaluate(H, ProjectiveVectors.raw(x₁), t₁)))

    state = State(H, steplength, x, t₁, t₀)
    cache = Cache(H, patch, predictor_corrector, state, x)
    options = Options(tol, corrector_maxiters, refinement_tol, refinement_maxiters, maxiters, fixedpatch)

    PathTracker(H, patch, predictor_corrector, steplength, state, options, x, cache)
end

function Base.show(io::IO, ::MIME"text/plain", tracker::PathTracker)
    print("PathTracker")
end

"""
     currt(tracker::PathTracker)

Current `t`.
"""
currt(tracker::PathTracker) = currt(tracker.state)
currt(state::State) = (1-state.t) * state.target + state.t * state.start

"""
     Δt(tracker::PathTracker)

Current steplength `Δt`.
"""
currΔt(tracker::PathTracker) = currΔt(tracker.state)
currΔt(state::State) = state.Δt * (state.target - state.start)

"""
     iters(tracker::PathTracker)

Current number of iterations.
"""
curriters(tracker::PathTracker) = curriters(tracker.state)
curriters(state::State) = state.iters

"""
     status(tracker::PathTracker)

Current status.
"""
currstatus(tracker::PathTracker) = currstatus(tracker.state)
currstatus(state::State) = state.status

"""
    currx(tracker::PathTracker)

Return the current value of `x`.
"""
currx(tracker::PathTracker) = tracker.x

"""
     tol(tracker::PathTracker)

Current tolerance.
"""
tol(tracker::PathTracker) = tracker.options.tol

"""
     refinement_tol(tracker::PathTracker)

Current refinement tolerance.
"""
refinement_tol(tracker::PathTracker) = tracker.options.refinement_tol

"""
     refinement_maxiters(tracker::PathTracker)

Current refinement maxiters.
"""
refinement_maxiters(tracker::PathTracker) = tracker.options.refinement_maxiters

"""
     corrector_maxiters(tracker::PathTracker)

Current correction maxiters.
"""
corrector_maxiters(tracker::PathTracker) = tracker.options.corrector_maxiters

"""
     set_tol!(tracker::PathTracker, tol)

Set the current tolerance to `tol`.
"""
function set_tol!(tracker::PathTracker, tol)
     tracker.options.tol = tol
     tol
end

"""
     set_refinement_maxiters!(tracker::PathTracker, tol)

Set the current refinement tolerance to `tol`.
"""
function set_refinement_tol!(tracker::PathTracker, tol)
     tracker.options.refinement_tol = tol
     tol
end

"""
     set_refinement_maxiters!(tracker::PathTracker, n)

Set the current refinement maxiters to `n`.
"""
function set_refinement_maxiters!(tracker::PathTracker, n)
     tracker.options.refinement_maxiters = n
     n
end

"""
     set_corrector_maxiters!(tracker::PathTracker, n)

Set the current correction maxiters to `n`.
"""
function set_corrector_maxiters!(tracker::PathTracker, n)
     tracker.options.corrector_maxiters = n
     n
end

"""
     fixedpatch!(tracker::PathTracker, flag::Bool)

Set whether the current (local) affine patch should be fixed and shouldn't be updatet anymore.
"""
function fixpatch!(tracker::PathTracker, flag)
     tracker.options.fixedpatch = flag
     nothing
end


"""
     PathTrackerResult(tracker)

Containing the result of a tracked path. The fields are
* `successfull::Bool` Indicating whether tracking was successfull.
* `returncode::Symbol` If the tracking was successfull then it is `:success`.
Otherwise the return code gives an indication what happened.
* `x::V` The result.
* `t::Float64` The `t` when the path tracker stopped.
* `res::Float64` The residual at `(x, t)`.
"""
struct PathTrackerResult{T, V<:AbstractVector}
     returncode::Symbol
     x::V
     t::T
     res::Float64
     iters::Int
end

function PathTrackerResult(tracker::PathTracker)
     PathTrackerResult(currstatus(tracker),
          copy(currx(tracker)), currt(tracker),
          currresidual(tracker), curriters(tracker))
end
