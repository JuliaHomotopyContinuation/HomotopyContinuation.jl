import ..Correctors
import ..Homotopies
import ..PredictionCorrection
import ..Predictors
import ..StepLength

using ..Utilities

export PathTracker

mutable struct Options
    tol::Float64
    corrector_maxiters::Int
    refinement_tol::Float64
    refinement_maxiters::Int
    maxiters::Int
    update_homotopies::Bool
end

mutable struct State{T, S<:StepLength.AbstractStepLengthState}
    # Our start point in space
    start::T
    # Our target point in space
    target::T
    steplength::S
    t::Float64 # The relative progress. `t` always goes from 1.0 to 0.0
    Δt::Float64 # Δt is the current relative step width
    iters::Int
    status::Symbol
end

function State(H::Homotopies.AbstractHomotopy, step, x₁::AbstractVector, t₁, t₀)
    checkstart(H, x₁)

    start, target = promote(t₁, t₀)
    steplength = StepLength.state(step, start, target)
    t = 1.0
    Δt = min(StepLength.relsteplength(steplength), 1.0)
    iters = 0

    status = :ok

    State(start, target, steplength, t, Δt, iters, status)
end

struct Cache{H<:Homotopies.HomotopyWithCache, PC<:PredictionCorrection.PredictorCorrectorCache, T}
    homotopy::H
    predictor_corrector::PC
    out::Vector{T}
end
function Cache(homotopy, predictor_corrector, state::State, x₁)
    H = Homotopies.HomotopyWithCache(homotopy, x₁, state.start)
    pc_cache = PredictionCorrection.cache(predictor_corrector, H, x₁, state.start)
    out = H(x₁, state.start)
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
    P<:Predictors.AbstractPredictor,
    Corr<:Correctors.AbstractCorrector,
    SL<:StepLength.AbstractStepLength,
    S<:State,
    C<:Cache,
    V<:AbstractVector}

    # these are fixed
    homotopy::H
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

function PathTracker(H::Homotopies.AbstractHomotopy, x₁::AbstractVector, t₁, t₀;
    corrector::Correctors.AbstractCorrector=Correctors.Newton(),
    predictor::Predictors.AbstractPredictor=Predictors.RK4(),
    steplength::StepLength.AbstractStepLength=StepLength.HeuristicStepLength(),
    tol=1e-7,
    refinement_tol=1e-11,
    corrector_maxiters::Int=2,
    refinement_maxiters=corrector_maxiters,
    maxiters=10_000,
    update_homotopies=true)

    predictor_corrector = PredictionCorrection.PredictorCorrector(predictor, corrector)
    # We have to make sure that the element type of x is invariant under evaluation
    x = similar(x₁, eltype(Homotopies.evaluate(H, x₁, t₁, Homotopies.cache(H, x₁, t₁))))

    state = State(H, steplength, x, t₁, t₀)
    cache = Cache(H, predictor_corrector, state, x)
    options = Options(tol, corrector_maxiters, refinement_tol, refinement_maxiters, maxiters, update_homotopies)

    PathTracker(H, predictor_corrector, steplength, state, options, x, cache)
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
struct PathTrackerResult{V<:AbstractVector, T}
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
