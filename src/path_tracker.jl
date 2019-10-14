export PathTracker,
       track!,
       track,
       PathTrackerStatus,
       is_at_infinity,
       is_success,
       is_tracking,
       is_failed,
       is_terminated_callback,
       is_invalid_startvalue,
       pathtracker,
       pathtracker_startsolutions,
       PathResult,
       solution,
       start_solution,
       accuracy,
       residual,
       winding_number,
       multiplicity,
       condition_jacobian,
       path_number,
       is_finite,
       is_singular,
       is_nonsingular,
       is_real,
       is_projective,
       is_affine

const pathtracker_supported_keywords = [
    :at_infinity_check,
    :endgame_start,
    :endgame_start_callback,
    :max_winding_number,
    :min_accuracy,
    :min_cond_eg,
    :min_step_size_before_eg,
    :min_step_size_eg,
    :s_always_consider_valuation,
    :samples_per_loop,
    :precision_strategy,
]


############
## STATUS ##
############

module PathTrackerStatus
"""
    enum PathTrackerStatus

The possible states a [`PathTracker`](@ref) can be in:

* `PathTrackerStatus.tracking`
* `PathTrackerStatus.success`
* `PathTrackerStatus.at_infinity`
* `PathTrackerStatus.excess_solution`
* `PathTrackerStatus.post_check_failed`
* `PathTrackerStatus.terminated_accuracy_limit`
* `PathTrackerStatus.terminated_ill_conditioned`
* `PathTrackerStatus.terminated_invalid_startvalue`
* `PathTrackerStatus.terminated_max_winding_number`
* `PathTrackerStatus.terminated_max_iters`
* `PathTrackerStatus.terminated_step_size_too_small`
"""
@enum states begin
    tracking
    success
    at_infinity
    terminated_accuracy_limit
    terminated_invalid_startvalue
    terminated_ill_conditioned
    terminated_max_iters
    terminated_max_winding_number
    terminated_step_size_too_small
    post_check_failed
    excess_solution
end
end

"""
    path_tracker_status(code::CoreTrackerStatus.states)

Construct a [`PathTrackerStatus.states`](@ref) from a [`CoreTrackerStatus.states`](@ref).
"""
function path_tracker_status(code::CoreTrackerStatus.states)
    if code == CoreTrackerStatus.success
        return PathTrackerStatus.success
    elseif code == CoreTrackerStatus.terminated_invalid_startvalue
        return PathTrackerStatus.terminated_invalid_startvalue
    elseif code == CoreTrackerStatus.terminated_maximal_iterations
        return PathTrackerStatus.terminated_max_iters
    elseif code == CoreTrackerStatus.terminated_step_size_too_small
        return PathTrackerStatus.terminated_step_size_too_small
    elseif code == CoreTrackerStatus.terminated_ill_conditioned
        return PathTrackerStatus.terminated_ill_conditioned
    else
        return PathTrackerStatus.tracking
    end
end

"""
    is_success(status::PathTrackerStatus.states)

Returns `true` if `status` indicates a success in tracking.
"""
is_success(status::PathTrackerStatus.states) = status == PathTrackerStatus.success

"""
    is_success(status::PathTrackerStatus.states)

Returns `true` if `status` indicates that a path diverged towards infinity.
"""
is_at_infinity(status::PathTrackerStatus.states) = status == PathTrackerStatus.at_infinity

"""
    is_tracking(status::PathTrackerStatus.states)

Returns `true` if `status` indicates the tracking is not going on.
"""
is_tracking(status::PathTrackerStatus.states) = status == PathTrackerStatus.tracking

"""
    is_invalid_startvalue(status::PathTrackerStatus.states)

Returns `true` if the provided start value was not valid.
"""
is_invalid_startvalue(status::PathTrackerStatus.states) =
    status == PathTrackerStatus.terminated_invalid_startvalue

"""
    is_terminated_callback(status::PathTrackerStatus.states)

Returns `true` if the provided callback indicated a termination of the path.
"""
is_terminated_callback(status::PathTrackerStatus.states) =
    status == PathTrackerStatus.terminated_callback

#############
## Options ##
#############
@enum PrecisionStrategy begin
    PREC_STRATEGY_NEVER
    PREC_STRATEGY_FINITE
    PREC_STRATEGY_ALWAYS
end

mutable struct PathTrackerOptions
    endgame_start::Float64
    at_infinity_check::Bool
    min_accuracy::Float64
    min_cond_eg::Float64
    min_step_size_before_eg::Float64
    min_step_size_eg::Float64
    precision_strategy::PrecisionStrategy
    s_always_consider_valuation::Float64
end

Base.show(io::IO, opts::PathTrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::PathTrackerOptions) = opts


###########
## STATE ##
###########

mutable struct PathTrackerState{AV<:AbstractVector}
    status::PathTrackerStatus.states
    valuation::Valuation
    prediction::AV
    solution::AV
    intermediate_solution::AV
    solution_accuracy::Float64
    solution_cond::Float64
    solution_residual::Float64
    winding_number::Union{Nothing,Int}
    max_winding_number_hit::Bool
    s::Float64
    eg_started::Bool
end

function PathTrackerState(x::AbstractVector; at_infinity_check::Bool = true)
    status = PathTrackerStatus.tracking
    valuation = Valuation(x; affine = at_infinity_check)
    prediction = copy(x)
    solution = copy(x)
    intermediate_solution = copy(x)
    solution_accuracy = NaN
    solution_cond = NaN
    solution_residual = NaN
    winding_number = nothing
    max_winding_number_hit = false
    s = 0.0
    eg_started = false
    PathTrackerState(
        status,
        valuation,
        prediction,
        solution,
        intermediate_solution,
        solution_accuracy,
        solution_cond,
        solution_residual,
        winding_number,
        max_winding_number_hit,
        s,
        eg_started,
    )
end

Base.show(io::IO, S::PathTrackerState) = print_fieldnames(io, S)
Base.show(io::IO, ::MIME"application/prs.juno.inline", S::PathTrackerState) = S

function init!(state::PathTrackerState, s::Float64)
    state.status = PathTrackerStatus.tracking
    init!(state.valuation)
    state.prediction .= 0.0
    state.solution .= 0.0
    state.intermediate_solution .= NaN
    state.solution_accuracy = NaN
    state.solution_cond = NaN
    state.solution_residual = NaN
    state.winding_number = nothing
    state.max_winding_number_hit = false
    state.s = s
    state.eg_started = false
end


##################
## PATH TRACKER ##
##################

abstract type AbstractPathTracker end
Base.broadcastable(T::AbstractPathTracker) = Ref(T)

"""
    PathTracker

The `PathTracker` combines the path tracking (with [`CoreTracker`](@ref)) with an *endgame*
routine. The *endgame*  is the name for special algorithms the end of the path tracking.
These enable to detect if a path is diverging or if it ends in a singular solution.
The `PathTracker` is more opinionated than the `CoreTracker` and implements additional logic
to handle numerical difficult solutions. In particular it always reparametrizes a solution
path to use a logarithmic time scale, i.e., ``x(t) → x(e^{-s})`` and we track for ``s`` from
``0`` to ``∞``.

In order to construct a `PathTracker` it is recommended to use the [`pathtracker`](@ref) and
[`pathtracker_startsolutions`](@ref) helper functions.
With a `PathTracker` constructed you can track a single path using the [`track`](@ref) method.
The result of this will be a [`PathResult`](@ref).

The `PathTracker` assumes that the provided homotopy `H` is defined in such a way that
`H(x,1)` is the start system and `H(x,0)` the target system.
During the path tracking an approximation of the valuation of a Puiseux series expansion of
the solution is computed. This is used to decide whether a path is diverging.
Before a certain treshold (`s_always_consider_valuation`) this approximation is only trusted
if a path gives some sign of ill-conditioning.
To compute singular solutions the *Cauchy endgame* is used which is based on
Cauchy's integral formula. For this we have to track solutions along a loop around the origin.
The number of loops necessary to arrive back at the start point is called the *winding number*.

## Options

The path tracker accepts all options which [`CoreTracker`](@ref) accepts. Furthermore the
following options are accepted:

* `at_infinity_check`: This is true if the provided start system is an affine polynomial system.
* `endgame_start` (default `2.0`): The value of `s` where the endgame is started earliest.
* `max_winding_number` (default `12`): The maximal winding number tried in the Cauchy endgame.
* `min_accuracy` (default `1e-5`): The `PathTracker` automatically lowers the desired
  accuracy automatically to `min_accuracy` if path tracking would fail othwerwise (since
  the desired cannot be reached)
* `min_cond_eg` (default `1e5`): The minimal condition number before the Cauchy endgame is
  started or a path is cut off.
* `min_step_size_before_eg` (default `exp2(-40)`): The minimal allowed step size before the
  endgame starts.
* `min_step_size_eg` (default `exp2(-120)`): The minimal allowed step size during the endgame.
  This is also control what is considered `∞` for the path tracking.
* `s_always_consider_valuation` (default `-log(1e-16)`) A threshold after which we always consider
  the valuation.
* `samples_per_loop (default `8`): The number of samples used during the endgame.
* `precision_strategy` (default `:adaptive_finite`): This controls whether `H(x,t)` is
  possibly evaluated with higher than machine precision during the endgame
  (the Jacobian is always computed with machine precision). The `:adaptive_finite` allows
  this only if we are optimistic that we can still obtain a finite solution.
  Other options are `:adaptive_never` where this is never allowed and `:adaptive_always`
  where it is always enabled.
  ``
"""
struct PathTracker{
    AV<:AbstractVector{Complex{Float64}},
    Prob<:AbstractProblem,
    CT<:CoreTracker{Float64,AV},
} <: AbstractPathTracker
    problem::Prob
    core_tracker::CT
    endgame::CauchyEndgame{AV}
    state::PathTrackerState{AV}
    options::PathTrackerOptions
    # We modify options of the core tracker dynamically. Therefore we need to be able
    # so start fresh for a new path
    default_ct_options::CoreTrackerOptions
end

function PathTracker(
    prob::AbstractProblem,
    x::AbstractVector{<:Number};
    min_step_size = 1e-30,
    kwargs...,
)

    core_tracker_kwargs, options = splitkwargs(kwargs, coretracker_supported_keywords)
    core_tracker = CoreTracker(
        prob,
        x;
        log_transform = true,
        predictor = Pade21(),
        min_step_size = min_step_size,
        core_tracker_kwargs...,
    )
    PathTracker(prob, core_tracker; options...)
end


function PathTracker(
    prob::Problem,
    core_tracker::CoreTracker;
    at_infinity_check = default_at_infinity_check(prob),
    endgame_start::Float64 = 2.0,
    max_winding_number::Int = 12,
    min_accuracy::Float64 = max(core_tracker.options.accuracy, 1e-5),
    min_cond_endgame::Float64 = 1e5,
    min_step_size_before_eg::Float64 = exp2(-40),
    min_step_size_eg = exp2(-120),
    s_always_consider_valuation::Float64 = -log(1e-16),
    samples_per_loop::Int = 8,
    precision_strategy::Symbol = :adaptive_finite,
)
    state = PathTrackerState(current_x(core_tracker); at_infinity_check = at_infinity_check)
    endgame = CauchyEndgame(
        current_x(core_tracker);
        samples_per_loop = samples_per_loop,
        max_winding_number = max_winding_number,
    )
    default_ct_options = copy(core_tracker.options)


    options = PathTrackerOptions(
        endgame_start,
        at_infinity_check,
        min_accuracy,
        min_cond_endgame,
        min_step_size_before_eg,
        min_step_size_eg,
        make_precision_strategy(precision_strategy),
        s_always_consider_valuation,
    )

    PathTracker(prob, core_tracker, endgame, state, options, default_ct_options)
end

default_at_infinity_check(prob::Problem{AffineTracking}) = true
default_at_infinity_check(prob::Problem{ProjectiveTracking}) = homvars(prob) !== nothing

always_true(_, _) = true

function make_precision_strategy(precision_strategy::Symbol)
    if precision_strategy == :adaptive_finite
        PREC_STRATEGY_FINITE
    elseif precision_strategy == :adaptive_never
        PREC_STRATEGY_NEVER
    elseif precision_strategy == :adaptive_always
        PREC_STRATEGY_ALWAYS
    else
        ArgumentError("Unknown argument `precision_strategy = $(precision_strategy)`." *
                      "Possible values are `:adaptive_finite`, `:adaptive_never`" *
                      " or `:adaptive_always`.") |> throw
    end
end

Base.show(io::IO, ::MIME"application/prs.juno.inline", PT::PathTracker) = PT
function Base.show(io::IO, S::PathTracker{AV}) where {AV}
    println(io, "PathTracker with solution type $AV")
end

seed(PT::PathTracker) = seed(PT.problem)

"""
    init!(tracker::PathTracker, x)

Prepare the `PathTracker`` `tracker` to track a path with start solution `x`.
"""
function init!(
    tracker::PathTracker,
    x;
    accuracy::Union{Nothing,Float64} = nothing,
    max_corrector_iters::Union{Nothing,Int} = nothing,
)
    options, ct_options = tracker.options, tracker.core_tracker.options
    copy!(ct_options, tracker.default_ct_options)
    if accuracy !== nothing
        ct_options.accuracy = accuracy
    end
    if max_corrector_iters !== nothing
        ct_options.max_corrector_iters = max_corrector_iters
    end
    if options.precision_strategy == PREC_STRATEGY_ALWAYS ||
       options.precision_strategy == PREC_STRATEGY_FINITE
        ct_options.precision = PRECISION_ADAPTIVE
    end

    init!(tracker.core_tracker, x, 0.0, options.endgame_start)
    init!(tracker.state, 0.0)

    # before endgame don't track the condition number
    # and we only update the limiting accuracy etc after a step failed
    ct_options.track_cond = false
    ct_options.steps_jacobian_info_update = -1
    ct_options.min_step_size = options.min_step_size_before_eg

    tracker
end

"""
    step!(tracker::PathTracker)

Perform a single step of the `PathTracker`.
"""
function step!(tracker::PathTracker)
    @unpack core_tracker, state, options, endgame = tracker
    step!(core_tracker)

    state.s = s = real(current_t(core_tracker))
    ct_status = status(core_tracker)
    if is_tracking(ct_status)
        state.eg_started || return state.status
        # If we didn't move forward there is nothing to do
        core_tracker.state.last_step_failed && return state.status

        # update the valuation
        update!(state.valuation, core_tracker)

        # We only care about the valuation if the path is starting to get worse
        # or we passed a certain threshold.
        cond_bad = cond(core_tracker) > options.min_cond_eg
        near_accuracy_limit = 1e4 * core_tracker.state.limit_accuracy > core_tracker.options.accuracy
        consider_always = s > options.s_always_consider_valuation
        if !(consider_always || cond_bad || near_accuracy_limit)
            return state.status
        end

        # Judge the current valuation to determine how to proceed
        verdict = judge(state.valuation; tol = 1e-2, tol_at_infinity = 1e-4)
        if verdict == VAL_AT_INFINITY &&
           (consider_always || cond_bad) && options.at_infinity_check
            state.status = PathTrackerStatus.at_infinity
        # If we expect a finite value and there is some ill conditioning let's do the
        # Cauchy endgame
        elseif verdict == VAL_FINITE && cond_bad
            # Perform endgame to estimate singular solution
            @label run_cauchy_eg
            retcode, m, p_accuracy = predict!(state.prediction, core_tracker, endgame)

            if retcode == CAUCHY_SUCCESS
                # We only accept a result of the Cauchy endgame only if two consecutive
                # loops resulted in the same number of loop
                m′ = state.winding_number
                if m′ === nothing
                    state.winding_number = m
                    state.solution .= state.prediction
                elseif m′ == m &&
                       norm(core_tracker)(state.solution, state.prediction) ≤ p_accuracy
                    @label cauchy_eg_success
                    state.status = PathTrackerStatus.success
                    state.solution_accuracy = p_accuracy
                    converged, s_accuracy, s_cond, s_res = check_converged!(
                        state.solution,
                        core_tracker,
                        state.prediction,
                        Inf;
                        tol = tracker.default_ct_options.accuracy,
                    )
                    state.solution_cond = s_cond
                    state.solution_residual = s_res
                    if converged && s_cond < 1e8 && s_accuracy < p_accuracy
                        state.solution_accuracy = s_accuracy
                    else
                        state.solution .= state.prediction
                    end
                else
                    state.solution .= state.prediction
                    state.winding_number = m
                end

            elseif retcode == CAUCHY_TERMINATED_ACCURACY_LIMIT
                if core_tracker.options.accuracy > options.min_accuracy
                    core_tracker.options.accuracy = options.min_accuracy
                    core_tracker.state.status = CoreTrackerStatus.tracking
                    @goto run_cauchy_eg
                elseif core_tracker.options.precision == PRECISION_FIXED_64
                    core_tracker.options.precision = PRECISION_ADAPTIVE
                    @goto run_cauchy_eg
                end

            elseif state.winding_number !== nothing
                @goto cauchy_eg_success

            elseif options.precision_strategy == PREC_STRATEGY_NEVER
                state.status = PathTrackerStatus.terminated_accuracy_limit
                # If we are here, we have a precision strategy which allows to use
                # higher precision. Therefore switch to it if this not yet happened
                # and retry to cauchy eg

            elseif retcode == CAUCHY_TERMINATED_MAX_WINDING_NUMBER
                # Terminate if we hit twice max winding number and make one step pause
                if state.max_winding_number_hit
                    state.status = PathTrackerStatus.terminated_max_winding_number
                else
                    state.max_winding_number_hit = true
                end

            elseif state.winding_number !== nothing
                @goto cauchy_eg_success

            else
                state.status = PathTrackerStatus.terminated_ill_conditioned
            end
        end

    # We split the tracking in two parts. Pre endgame and endgame. Therefore we have to
    # handle the success case twice
    elseif is_success(ct_status) && !state.eg_started
        state.intermediate_solution .= current_x(core_tracker)
        state.eg_started = true
            # set min_step size
        core_tracker.options.min_step_size = options.min_step_size_eg
            # setup path tracker to continue tracking
        t₀ = -log(4 * options.min_step_size_eg)
        init!(
            tracker.core_tracker,
            current_x(core_tracker),
            options.endgame_start,
            t₀;
            loop = true,
        )
            # start keeping track of the condition number during the endgame
        tracker.core_tracker.options.track_cond = true
            # update the limiting accuracy and cond every other step
        tracker.core_tracker.options.steps_jacobian_info_update = 2
            # update the valuation
        update!(state.valuation, core_tracker)

            # If our strategy is to allow adaptive precision only for finite values
            # then we start with not allowing adaptive precision during the endgame
            # and wait to hit the accuracy limit.
        if options.precision_strategy == PREC_STRATEGY_FINITE
            tracker.core_tracker.options.precision = PRECISION_FIXED_64
        end

    elseif is_success(ct_status) ||
           ct_status == CoreTrackerStatus.terminated_step_size_too_small
        converged, s_acc, s_cond, s_res = check_converged!(
            state.solution,
            core_tracker,
            current_x(core_tracker),
            Inf;
            tol = tracker.default_ct_options.accuracy,
        )
        if converged
            state.status = PathTrackerStatus.success
            state.solution_accuracy = s_acc
            state.solution_cond = s_cond
            state.solution_residual = s_res
            state.s = Inf
        elseif ct_status == CoreTrackerStatus.terminated_step_size_too_small
            state.status = PathTrackerStatus.terminated_step_size_too_small
        else
            state.status = PathTrackerStatus.post_check_failed
        end

    elseif ct_status == CoreTrackerStatus.terminated_accuracy_limit
        # First update the valuation  and check whether we are good
        update!(state.valuation, core_tracker)
        # we put less requirement on tol_at_infinity
        verdict = judge(state.valuation; tol = 1e-2, tol_at_infinity = 1e-4)
        # We also terminate at this point if the valuation indicates that we are
        # substantially less than 0
        if verdict == VAL_AT_INFINITY && options.at_infinity_check
            state.status = PathTrackerStatus.at_infinity
        # Now, we have to differentiate 3 different cases:
        # 1) Current accuracy is larger than the defined minimal accuracy
        #     -> decrease accuracy to minimal accuracy
        # 2) We are already at minimal accuracy and allow adaptive precision for finite
        #    values
        #  a) it looks like we could end up with something finite
        #     -> enable adaptive precision
        #  b) else
        #     -> terminate
        # 3) We are already at minimal accuracy and don't allow adaptive precision
        #     -> terminate
        # 1)
        elseif core_tracker.options.accuracy > options.min_accuracy
            # TODO: Currently we never increase this again
            core_tracker.options.accuracy = options.min_accuracy
            core_tracker.state.status = CoreTrackerStatus.tracking
        # 2)
        elseif options.precision_strategy == PREC_STRATEGY_FINITE
            # 2.a) We say that a path could still become finite if it is not classified
            #      as going to infinity for a loose tolerance
            #      Also if the current valuation is less than -1 (so significantly less than 0)
            #      we do not continue the tracking
            loose_verdict = judge(state.valuation; tol = 1e-2, tol_at_infinity = 1e-2)
            if loose_verdict != VAL_AT_INFINITY
                core_tracker.options.precision = PRECISION_ADAPTIVE
                core_tracker.state.status = CoreTrackerStatus.tracking
            # 2.b)
            else
                state.status = PathTrackerStatus.terminated_accuracy_limit
            end
        # 3)
        else
            state.status = PathTrackerStatus.terminated_accuracy_limit
        end
    elseif state.winding_number !== nothing
        p_accuracy = 1e-5
        @goto cauchy_eg_success
    else
        state.status = path_tracker_status(ct_status)
    end

    if !is_success(state.status) && !is_tracking(state.status)
        state.solution .= current_x(core_tracker)
        state.solution_cond = LA.cond(core_tracker)
    end


    state.status
end

"Update the valuation with the current state."
update!(val::Valuation, T::CoreTracker) =
    update!(val, T.state.x, T.state.ẋ, real(current_t(T)), T.predictor)


"""
    check_converged!(y, tracker::CoreTracker, x::AbstractVector, t)

Check that the CoreTracker is converged also converged at `(x,t)`. Returns a named tuple
`(converged, accuracy, cond)` and stores a (possibly) refined solution into `y`.
"""
function check_converged!(
    y,
    tracker::CoreTracker,
    x::AbstractVector,
    t::Number;
    tol::Float64 = 1e-6,
)
    init!(tracker.state.norm, x)
    # Make sure to evaluate the Jacobian only *once*. Otherwise it can happen at singuar
    # solutions that we bounce away from a good solution.
    result = correct!(y, tracker, x, t; tol = tol, max_iters = 3, full_steps = 1)
    if is_converged(result)
        tracker.state.residual .= abs.(tracker.corrector.r)
    else
        evaluate!(tracker.corrector.r, tracker.homotopy, x, t)
        tracker.state.residual .= abs.(tracker.corrector.r)
    end
    res = LA.norm(tracker.state.residual, InfNorm())
    # Fix the residual to the norm otherwise we can scale away zero rows.
    tracker.state.residual .= res
    # If the accuracy is NaN we assume that we dived by 0 -> Jacobian singular
    if isnan(result.accuracy)
        cond_jac = Inf
    else
        cond_jac = cond!(tracker.state.jacobian, tracker.state.norm, tracker.state.residual)
    end
    (
     converged = is_converged(result),
     accuracy = result.accuracy,
     cond = cond_jac,
     res = res,
    )
end

function Base.iterate(tracker::PathTracker, state::Int = 0)
    if is_tracking(tracker.state.status)
        step!(tracker)
        tracker, state + 1
    else
        nothing
    end
end

####################
## QUERYING STATE ##
####################

status(tracker::PathTracker) = tracker.state.status

"""
    solution(tracker::PathTracker)

Obtain the solution computed by the `PathTracker`.
"""
solution(tracker::PathTracker) = pull_back(tracker.problem, tracker.state.solution)

"""
    norm(tracker::PathTracker)

Obtain the norm used by the `PathTracker`.
"""
LA.norm(tracker::PathTracker) = LA.norm(tracker.core_tracker)


"""
    winding_number(tracker::PathTracker)

Obtain the estimate of the winding number computed by the `PathTracker`. Returns `nothing`
if the Cauchy endgame was not run.
"""
winding_number(tracker::PathTracker) = winding_number(tracker.state)
winding_number(state::PathTrackerState) = state.winding_number

###############################
## Convencience constructors ##
###############################

const pathtracker_startsolutions_supported_keywords = [
    problem_startsolutions_supported_keywords
    coretracker_supported_keywords
    pathtracker_supported_keywords
]


"""
    pathtracker_startsolutions(args...; kwargs...)

Construct a [`PathTracker`](@ref) and start solutions in the same way [`solve`](@ref) does it.
This also takes the same input arguments as `solve`. This is convenient if you want
to investigate single paths.
"""
function pathtracker_startsolutions(args...; kwargs...)
    invalid = invalid_kwargs(kwargs, pathtracker_startsolutions_supported_keywords)
    check_kwargs_empty(invalid, pathtracker_startsolutions_supported_keywords)
    supported, rest = splitkwargs(kwargs, problem_startsolutions_supported_keywords)
    prob, startsolutions = problem_startsolutions(args...; supported...)
    construct_tracker(prob, startsolutions; rest...), startsolutions
end

function construct_tracker(prob::Problem, startsolutions; kwargs...)
    PathTracker(prob, start_solution_sample(startsolutions); kwargs...)
end



"""
    pathtracker(args...; kwargs...)

Construct a [`PathTracker`](@ref) in the same way [`solve`](@ref) does it.
This also takes the same input arguments as `solve` with the exception that you do not need
to specify startsolutions.

## Examples

### Obtain single solution
We want to construct a path tracker to track a parametrized system `f` with parameters `p`
from the parameters `a` to `b`.
```julia
tracker = pathtracker(f, parameters=p, p₁=a, p₀=b)
```
You then can obtain a single solution at `b` by using
```julia
x_b = solution(track(tracker, x_a))
```
"""
pathtracker(args...; kwargs...) = first(pathtracker_startsolutions(args...; kwargs...))

#################
## PATH RESULT ##
#################

"""
    PathResult{V<:AbstractVector}

A `PathResult` is the result of tracking of a path with [`track`](@ref) using a
[`PathTracker`](@ref).

# Fields
* `return_code`: See the list of return codes below.
* `solution::V`: The solution vector.
* `t::Float64`: The value of `t` at which `solution` was computed. Note that if
  `return_code` is `:at_infinity`, then `t` is the value when this was decided.
* `accuracy::Float64`: An approximation of ``||x-x̄||`` where ``x`` is the computed solution
  and ``x̄`` is the true solution.
* `residual::Float64`: The value of the infinity-norm of `H(solution, 0)`.
* `multiplicity::Union{Nothing, Int}` is the multiplicity of the `solution`. This is
  only assigned by [`singular`](@ref).
* `condition_jacobian::Union{Nothing, Float64}`: This is the condition number of the
  row-equilibrated Jacobian at the solution. A high condition number indicates a singularity.
* `winding_number:Union{Nothing, Int}`: The estimated winding number. This is a lower bound
  on the multiplicity of the solution.
* `start_solution::Union{Nothing, Int}`: The start solution of the path.
* `accepted_steps::Int`: The number of accepted steps during the path tracking.
* `rejected_steps::Int`: The number of rejected steps during the path tracking.
* `valuation::Union{Nothing, Vector{Float64}}`: An approximation of the valuation of the
  Puiseux series expansion of `x(t)`.
* `valuation_accuracy::Union{Nothing, Vector{Float64}}`: An estimate of the accuracy of the
  valuation of the Puiseux series expansion of `x(t)`.

# Return codes

These is the list of possible return codes:

* `:success`: The `PathTracker` obtained a solution.
* `:at_infinity`: The `PathTracker` stopped the tracking of the path since it determined
  that that path is diverging towards infinity.
* `:terminated_callback`: One of the optional `PathTracker` callbacks terminated the tracking.
* `:terminated_max_iters`: The `PathTracker` terminated since it reached the limit accuracy.
* `:terminated_invalid_startvalue`: The `PathTracker` terminated since the provided start
  value is invalid.
* `:terminated_step_size_too_small`: The `PathTracker` terminated since the step size
  became smaller than the provided threshold.
* `:terminated_accuracy_limit`: The `PathTracker` terminated since the problem was too
  ill-conditioned to be tracked further with the desired minimal accuracy.
* `:terminated_ill_conditioned`: The `PathTracker` terminated since the Jacobian of the
  homotopy was too ill-conditioned.
* `:post_check_failed`: The verification of a non-singular solution failed.

# Constructors


     PathResult(tracker::PathTracker,
                start_solution=nothing,
                path_number::Union{Nothing,Int}=nothing;
                details=:default)

Construct a `PathResult` using the current state of the `PathTracker`.
Possible values for `details` are `:minimal` (minimal details), `:default` (default) and
`:extensive` (all information possible).
"""
struct PathResult{V<:AbstractVector}
    return_code::Symbol
    solution::V
    t::Float64
    accuracy::Float64
    residual::Float64
    condition_jacobian::Float64
    winding_number::Union{Nothing,Int}
    path_number::Base.RefValue{Union{Nothing,Int}}
    start_solution::Union{Nothing,V} # level 1+
    intermediate_solution::Union{Nothing,V}
    multiplicity::Base.RefValue{Union{Nothing,Int}} # only assigned by singular(Result)
    # performance stats
    accepted_steps::Int
    rejected_steps::Int
    valuation::Union{Nothing,Vector{Float64}} # level 2+
    valuation_accuracy::Union{Nothing,Vector{Float64}} # level 2+
end

function PathResult(
    tracker::PathTracker,
    start_solution::Union{Nothing,AbstractVector} = nothing,
    path_number::Union{Nothing,Int} = nothing;
    details::Symbol = :default,
)
    @unpack state, core_tracker = tracker
    details_level = detailslevel(details)

    return_code = Symbol(state.status)
    x = solution(tracker)
    t = return_code == :success ? 0.0 : exp(-state.s)
    accuracy = state.solution_accuracy
    residual = state.solution_residual
    condition_jac = state.solution_cond
    windingnumber = state.winding_number
    # this needs to be assigned manually
    multiplicity = Base.RefValue{Union{Nothing,Int}}(nothing)

    startsolution = nothing
    if details_level ≥ 1 && start_solution !== nothing
        # mimic the behaviour in track! to get a start solution of the same type as x
        embed!(core_tracker.state.x, tracker.problem, start_solution)
        startsolution = pull_back(tracker.problem, core_tracker.state.x; regauge = false)
    end
    intermediate_solution = nothing
    if !isnan(state.intermediate_solution[1])
        intermediate_solution = pull_back(tracker.problem, state.intermediate_solution)
    end

    accepted_steps = core_tracker.state.accepted_steps
    rejected_steps = core_tracker.state.rejected_steps

    valuation = valuation_accuracy = nothing
    if details_level == 2
        valuation = copy(tracker.state.valuation.ν)
        valuation_accuracy = abs.(tracker.state.valuation.ν̇)
    end

    PathResult(
        return_code,
        x,
        t,
        accuracy,
        residual,
        condition_jac,
        windingnumber,
        Base.RefValue{Union{Nothing,Int}}(path_number),
        startsolution,
        intermediate_solution,
        multiplicity,
        accepted_steps,
        rejected_steps,
        valuation,
        valuation_accuracy,
    )
end

function detailslevel(details::Symbol)
    if details == :minimal
        return 0
    elseif details == :extensive
        return 2
    else
        return 1
    end
end

"""
    result_type(tracker::PathTracker)

Returns the type of result `track` will return.
"""
result_type(tracker::PathTracker) = PathResult{typeof(solution(tracker))}

function Base.show(io::IO, r::PathResult{AV}) where {AV}
    iscompact = get(io, :compact, false)
    if iscompact || haskey(io, :typeinfo)
        println(io, "• return_code → :$(r.return_code)")
        if r.return_code ≠ :success
            println(io, " • t → $(r.t)")
        end
        println(io, " • solution → ", r.solution)
        r.accuracy !== nothing && println(
            io,
            " • accuracy → $(Printf.@sprintf "%.3e" r.accuracy)",
        )
        winding_number(r) !== nothing && println(
            io,
            " • winding_number → $(winding_number(r))",
        )
        if multiplicity(r) !== nothing
            println(io, " • multiplicity → $(multiplicity(r))")
        end
        path_number(r) !== nothing && println(io, " • path_number: ", path_number(r))
    else
        header = "PathResult{$AV}"
        compact_io = IOContext(io, :compact => true)
        println(io, "PathResult{$AV}")
        println(io, "="^length(header))
        path_number(r) !== nothing && println(io, " • path_number → ", path_number(r))
        println(io, " • return_code → :$(r.return_code)")
        if r.return_code ≠ :success
            println(compact_io, " • t → $(r.t)")
        end
        println(compact_io, " • solution → ", r.solution)
        r.accuracy !== nothing && println(compact_io, " • accuracy → ", r.accuracy)
        r.residual !== nothing && println(compact_io, " • residual → ", r.residual)
        if r.condition_jacobian !== nothing
            println(compact_io, " • condition_jacobian → ", r.condition_jacobian)
        end

        if multiplicity(r) !== nothing
            println(io, " • multiplicity → $(multiplicity(r))")
        end
        if winding_number(r) !== nothing
            println(io, " • winding_number → $(winding_number(r))")
        end
        if r.valuation !== nothing
            println(compact_io, " • valuation → ", r.valuation)
        end
        if r.valuation_accuracy !== nothing
            println(compact_io, " • valuation_accuracy → ", r.valuation_accuracy)
        end
    end
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::PathResult) = x

"""
    solution(r::PathResult)

Get the solution of the path.
"""
solution(r::PathResult) = r.solution


"""
    start_solution(r::PathResult)

Get the start solution of the path.
"""
start_solution(r::PathResult) = r.start_solution

"""
    accuracy(r::PathResult)

Get the accuracy of the solution. This is an estimate of the (relative) distance to the
true solution.
"""
accuracy(r::PathResult) = r.accuracy

"""
    residual(r::PathResult)

Get the residual of the solution ``x`` of the path, i.e., ``||H(x, 0)||₂``.
"""
residual(r::PathResult) = r.residual

"""
    winding_number(r::PathResult)

Get the winding number of the solution of the path. Returns `nothing` if it wasn't computed.
"""
winding_number(r::PathResult) = r.winding_number

"""
    multiplicity(P::PathResult)

Returns the multiplicity of `P`.
"""
multiplicity(P::PathResult) = P.multiplicity[]

"""
    condition_jacobian(r::PathResult)

Return the condition number of the Jacobian of the result.
"""
condition_jacobian(r::PathResult) = r.condition_jacobian

"""
    cond(r::PathResult)

Return the condition number of the Jacobian of the result.
"""
LA.cond(r::PathResult) = r.condition_jacobian

"""
    path_number(r::PathResult)

The number of the path.
"""
path_number(r::PathResult) = r.path_number[]

"""
    is_success(r::PathResult)

Checks whether the path is successfull.
"""
is_success(r::PathResult) = r.return_code == :success

"""
    is_failed(r::PathResult)

Checks whether the path failed.
"""
is_failed(r::PathResult) = !(r.return_code == :at_infinity || r.return_code == :success)

"""
    is_at_infinity(r::PathResult)

Checks whether the path goes to infinity.
"""
is_at_infinity(r::PathResult) = r.return_code == :at_infinity


"""
    is_finite(r::PathResult)

Checks whether the path result is finite.
"""
is_finite(r::PathResult) = r.return_code == :success

# Base fallback
Base.isfinite(r::PathResult) = is_finite(r)

"""
    is_singular(r::PathResult; tol=1e10)

Checks whether the path result is singular. This is true if
the multiplicity is larger than  1 or if the condition number of the Jacobian
is larger than `tol`.
"""
is_singular(r::PathResult; tol = 1e10) = is_singular(r, tol)
function is_singular(r::PathResult, tol::Real)
    (unpack(r.condition_jacobian, 1.0) > tol ||
     unpack(multiplicity(r), 1) > 1 || unpack(winding_number(r), 1) > 1) && is_success(r)
end

"""
    is_nonsingular(r::PathResult; tol=1e10)

Checks whether the path result is non-singular. This is true if
it is not singular.
"""
is_nonsingular(r::PathResult; kwargs...) = !is_singular(r; kwargs...) && is_success(r)
is_nonsingular(r::PathResult, tol::Real) = !is_singular(r, tol) && is_success(r)


"""
    is_real(r::PathResult; tol=1e-6)

We consider a result as `real` if the 2-norm of the imaginary part of the solution is at most `tol`.
"""
is_real(r::PathResult; tol = 1e-6) = is_real(r, tol)
is_real(r::PathResult, tol::Real) = is_real_vector(r.solution, tol)
# provide fallback since this in in Base
Base.isreal(r::PathResult, tol) = is_real(r, tol)
Base.isreal(r::PathResult; kwargs...) = is_real(r; kwargs...)

"""
    is_projective(r::PathResult)

Return`s true if the solution is a projective vector.
"""
is_projective(r::PathResult) = false
is_projective(r::PathResult{<:PVector}) = true

"""
    is_affine(pathresult)

Return`s true if the solution is an affine vector.
"""
is_affine(r::PathResult) = !is_projective(r)

###########
## track ##
###########

"""
    track(tracker::PathTracker, x₁;
            path_number=nothing,
            details::Symbol=:default,
            start_parameters = nothing,
            target_parameters = nothing)

Track the path `x(t)` with start solution `x₁` from ``1`` towards ``0``.
Returns a [`PathResult`](@ref).

The `details` options controls the level of details of the informations available
in the [`PathResult`](@ref).
If `tracker` uses a parameter homotopy you can set the start and target parameters
by setting the corresponding fields.
To investigate the behaviour of a particular take a look at [`path_info`](@ref).
"""
function track(
    tracker::PathTracker,
    x,
    t₁ = nothing,
    t₀ = nothing;
    path_number::Union{Int,Nothing} = nothing,
    details::Symbol = :default,
    accuracy::Union{Nothing,Float64} = nothing,
    max_corrector_iters::Union{Nothing,Int} = nothing,
    start_parameters = nothing,
    target_parameters = nothing,
)
    if t₁ !== nothing
        @warn(
            "track(tracker, x, t₁, t₀) is deprecated. `t₁` is always assumed to be `1` " *
            " and `t₀` set automatically.",
        )
    end

    start_parameters !== nothing && start_parameters!(tracker, start_parameters)
    target_parameters !== nothing && target_parameters!(tracker, target_parameters)

    track!(tracker, x; accuracy = accuracy, max_corrector_iters = max_corrector_iters)
    PathResult(tracker, x, path_number; details = details)
end

"""
    start_parameters!(tracker::PathTracker, p)

Set the start parameters of the homotopy in in `tracker` to `p`.
"""
start_parameters!(T::PathTracker, p) = start_parameters!(T.core_tracker, p)

"""
    target_parameters!(tracker::PathTracker, p)

Set the target parameters of the homotopy in in `tracker` to `p`.
"""
target_parameters!(T::PathTracker, p) = target_parameters!(T.core_tracker, p)


"""
    track!(tracker::PathTracker, x₁)::PathTrackerStatus.states

Track the path `x(t)` with start solution `x₁` from ``1`` towards ``0``.
Returns a [`PathTrackerStatus.states`](@ref).
"""
function track!(
    tracker::PathTracker,
    x;
    accuracy::Union{Nothing,Float64} = nothing,
    max_corrector_iters::Union{Nothing,Int} = nothing,
)
    init!(tracker, x; accuracy = accuracy, max_corrector_iters = max_corrector_iters)
    while is_tracking(tracker.state.status)
        step!(tracker)
    end
    tracker.state.status
end
