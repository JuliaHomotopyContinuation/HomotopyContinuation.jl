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

function Options(;tol=1e-6,
    refinement_tol=1e-8,
    corrector_maxiters::Int=2,
    refinement_maxiters=corrector_maxiters,
    maxiters=500,
    initial_steplength=0.1,
    minimal_steplength=1e-14)

    Options(tol, corrector_maxiters, refinement_tol, refinement_maxiters, maxiters,
            initial_steplength, minimal_steplength)
end
Base.show(io::IO, opts::Options) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/juno+inline", opts::Options) = opts


###########
# State
##########

@enum Status
    done
    tracking
    terminated_maximal_iterations
    terminated_invalid_startvalue
    terminated_steplength_too_small
end

mutable struct State{T, PV<:ProjectiveVectors.AbstractProjectiveVector{T}}
    x::PV # current x
    x̂::PV # last prediction
    x̄::PV # canidate for new x
    ẋ::Vector{T} # derivative at current x
    segment::ComplexSegment
    s::Float64 # current step length (0 ≤ s ≤ length(segment))
    Δs::Float64 # current step size
    accuracy::Float64
    iters::Int
    status::Status
end

function State(x₁::ProjectiveVectors.AbstractProjectiveVector, t₁, t₀, options::Options)
    x, x̂, x̄ = copy(x₁), copy(x₁), copy(x₁)
    ẋ = copy(x₁.data)
    segment = ComplexSegment(promote(t₁, t₀)...)
    s = 0.0
    Δs = convert(Float64, min(options.initial_steplength, length(segment)))
    accuracy = 0.0
    iters = 0
    status = tracking
    State(x, x̂, x̄, ẋ, segment, s, Δs, accuracy, iters, status)
end

Base.show(io::IO, state::State) = print_fieldnames(io, state)
Base.show(io::IO, ::MIME"application/juno+inline", state::State) = state


###########
# Cache
##########
struct Cache{H<:Homotopies.HomotopyWithCache, P<:Predictors.AbstractPredictorCache, C<:Correctors.AbstractCorrectorCache, T}
    homotopy::H
    predictor::P
    corrector::C
    J::Matrix{T}
    out::Vector{T}
end
function Cache(H::Homotopies.HomotopyWithCache, predictor, corrector, state::State)
    t = state.segment[state.s]
    pcache = Predictors.cache(predictor, H, state.x, t)
    ccache = Correctors.cache(corrector, H, state.x, t)
    J = Homotopies.jacobian(H, state.x, t)
    out = H(state.x, t)
    Cache(H, pcache, ccache, J, out)
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
    S<:State, C<:Cache}
    # these are fixed
    homotopy::H
    predictor::Predictor
    corrector::Corrector
    # these are mutable
    state::S
    options::Options
    cache::C
end

"""
    PathTracker(problem::Problems.AbstractProblem, x₁, t₁, t₀; kwargs...)

Construct a [`PathTracking.PathTracker`](@ref) from the given `problem`.
"""
function PathTracker(prob::Problems.Projective, x₁, t₁, t₀; patch=AffinePatches.OrthogonalPatch(), kwargs...)
    y₁ = Problems.embed(prob, x₁)
    H = Homotopies.PatchedHomotopy(prob.homotopy, AffinePatches.state(patch, y₁))
    PathTracker(H, y₁, complex(t₁), complex(t₀); kwargs...)
end
function PathTracker(H::Homotopies.AbstractHomotopy, x₁::ProjectiveVectors.AbstractProjectiveVector, t₁, t₀;
    corrector::Correctors.AbstractCorrector=Correctors.Newton(),
    predictor::Predictors.AbstractPredictor=Predictors.Euler(), kwargs...)

    options = Options(;kwargs...)
    HC = Homotopies.HomotopyWithCache(H, x₁, t₁)
    # We have to make sure that the element type of x is invariant under evaluation
    indempotent_x = begin
        u = Vector{Any}(undef, size(H)[1])
        Homotopies.evaluate!(u, HC, x₁, t₁)
        similar(x₁, promote_type(typeof(u[1]), ComplexF64))
    end
    state = State(indempotent_x, t₁, t₀, options)
    cache = Cache(HC, predictor, corrector, state)

    PathTracker(H, predictor, corrector, state, options, cache)
end

Base.show(io::IO, ::PathTracker) = print(io, "PathTracker()")
Base.show(io::IO, ::MIME"application/juno+inline", x::PathTracker) = x


include("path_tracking/getter_setter.jl")

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
     returncode::Symbol
     x::V
     t::T
     accuracy::Float64
     iters::Int
end

function PathTrackerResult(tracker::PathTracker)
    state = tracker.state
     PathTrackerResult(state.status,
          copy(state.x), state.segment[state.s],
          state.accuracy,
          state.iters)
end

Base.show(io::IO, result::PathTrackerResult) = print_fieldnames(io, result)
Base.show(io::IO, ::MIME"application/juno+inline", result::PathTrackerResult) = result

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
     if retcode == :success
         x₀ .= currx(tracker)
     end
     retcode
end
function track!(tracker::PathTracker, x₁, t₁, t₀; kwargs...)
    setup!(tracker, x₁, t₁, t₀; kwargs...)

    while tracker.state.status == tracking
        step!(tracker)
        check_terminated!(tracker)
    end
    if tracker.state.status == done
        refine!(tracker)
    end
    tracker.state.status
end

"""
    setup!(pathtracker, x₁, t₁, t₀, checkstartvalue=true))

Setup `pathtracker` to track `x₁` from `t₁` to `t₀`. Use this if you want to use the
pathtracker as an iterator.
"""
function setup!(tracker::PathTracker, x₁::AbstractVector, t₁, t₀; checkstartvalue=true, compute_ẋ=true)
    state = tracker.state
    reset!(tracker.state, x₁, t₁, t₀, tracker.options)
    Homotopies.update!(tracker.cache.homotopy, state.x, t₁)
    Predictors.reset!(tracker.cache.predictor, state.x, t₁)

    checkstartvalue && checkstartvalue!(tracker)
    compute_ẋ && compute_ẋ!(tracker)
    Predictors.setup!(tracker.cache.predictor, state.x, state.ẋ, currt(state))
    nothing
end

function reset!(state::State, x₁::AbstractVector, t₁, t₀, options::Options)
    state.segment = ComplexSegment(promote(t₁, t₀)...)
    state.s = 0.0
    state.Δs = min(options.initial_steplength, length(state.segment))
    state.accuracy = 0.0
    state.iters = 0
    state.status = tracking
    ProjectiveVectors.embed!(state.x, x₁)
end

patch(cache::Cache) = cache.homotopy.homotopy.patch # this is a little bit hacky...

function checkstartvalue!(tracker)
    result = correct!(tracker.state.x̄, tracker)
    if Correctors.converged(result)
        tracker.state.x .= tracker.state.x̄
    else
        tracker.state.status = terminated_invalid_startvalue
    end
    nothing
end

function compute_ẋ!(tracker)
    ẋ = tracker.state.ẋ
    @inbounds Homotopies.jacobian_and_dt!(tracker.cache.J, ẋ,
                        tracker.cache.homotopy, tracker.state.x, currt(tracker))
    @inbounds for i in eachindex(ẋ)
        ẋ[i] = -ẋ[i]
    end
    solve!(tracker.cache.J, ẋ)
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
    t = currt(state)
    try
        last_step_failed = false
        Δt = currΔt(state)
        while state.iters < options.maxiters
            state.iters += 1

            Predictors.predict!(x̂, tracker.predictor, cache.predictor, H, x, t, Δt, ẋ)
            result = Correctors.correct!(x̄, tracker.corrector, cache.corrector, H, x̂, t + Δt, options.tol, 100)
            # @show result
            if Correctors.converged(result)

                x .= x̄
                state.s += state.Δs
                Homotopies.update!(H, x, t + Δt)
                state.accuracy = result.accuracy

                # Compute the new derivate at t + Δt
                compute_ẋ!(tracker)
                Predictors.update!(cache.predictor, x, ẋ, t + Δt)

                # Step size change
                update_stepsize!(state, result,
                    tracker.predictor,
                    cache.predictor, last_step_failed)

                return
            else
                last_step_failed = true
                # Step failed, so we have to try with a new (smaller) step size
                update_stepsize!(state, result,
                    tracker.predictor,
                    cache.predictor, last_step_failed)
                Δt = currΔt(state)
                if state.Δs < options.minimal_steplength
                    state.status = terminated_steplength_too_small
                    return
                end
            end
        end
    catch err
        if !(err isa LinearAlgebra.SingularException)
            rethrow(err)
        end
    end
    nothing
end


const gθ̄ = sqrt(2) - 1
g(Θ) = sqrt(1+4Θ) - 1
function update_stepsize!(state::State, result::Correctors.Result,
                          predictor::Predictors.AbstractPredictor,
                          predictor_cache::Predictors.AbstractPredictorCache,
                          last_step_failed)
    # P.247 + 248
    correction_factor = 1.0
    if Correctors.converged(result)
        θ_min = 1e-6
        if result.Θ₀ < θ_min
            correction_factor = gθ̄ / (2g(θ_min))
        else
            correction_factor =
                (result.norm_Δx₀ * gθ̄ /
                 (euclidean_distance(state.x̂, state.x̄) * result.Θ₀))
        end
    else # Correct step size with a-posteriori estimate
        correction_factor = gθ̄ / g(result.Θᵢ)
    end


    Δs′ = 0.9 * Predictors.asymptotic_correction(predictor, predictor_cache, correction_factor, state.Δs)

    state.Δs = @fastmath min(Δs′, length(state.segment) - state.s)
end


function check_terminated!(tracker)
    if tracker.state.s ≈ length(tracker.state.segment)
        tracker.state.status = done
    elseif tracker.state.iters ≥ tracker.options.maxiters
        tracker.state.status = terminated_maximal_iterations
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
    if Correctors.converged(result)
        tracker.state.x .= tracker.state.x̄
        tracker.state.accuracy = result.accuracy
    end
    nothing
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
    tracker = PathTracking.PathTracker(prob, Utilities.start_solution_sample(startsolutions), one(ComplexF64), zero(ComplexF64); rest...)

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


end
