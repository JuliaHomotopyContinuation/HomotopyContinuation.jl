export PathTracker,
       PathTrackerStatus,
       is_at_infinity,
       is_success,
       is_tracking,
       PathResult,
       pathtracker,
       pathtracker_startsolutions,
       solution,
       winding_number

const pathtracker_supported_keywords = [
    :at_infinity_check,
    :samples_per_loop,
    :max_winding_number,
    :overdetermined_min_accuracy,
    :overdetermined_min_residual,
    :cond_eg_start,
    :min_cond_at_infinity,
    :t_eg_start,
    :tol_val_inf_accurate,
    :tol_val_finite_accurate,
    :accuracy,
    :accuracy_eg,
]


############
## STATUS ##
############

"""
    enum PathTrackerStatus

The possible states a [`PathTracker`](@ref) can be in:
* `PT_TRACKING`
* `PT_SUCCESS`
* `PT_AT_INFINITY`
* `PT_EXCESS_SOLUTION`
* `PT_POST_CHECK_FAILED`
* `PT_TERMINATED_ACCURACY_LIMIT`
* `PT_TERMINATED_CALLBACK`
* `PT_TERMINATED_ILL_CONDITIONED`
* `PT_TERMINATED_INVALID_STARTVALUE`
* `PT_TERMINATED_MAX_ITERS`
* `PT_TERMINATED_STEP_SIZE_TOO_SMALL`
"""
@enum PathTrackerStatus begin
    PT_TRACKING
    PT_SUCCESS
    PT_AT_INFINITY
    PT_TERMINATED_CALLBACK
    PT_TERMINATED_MAX_ITERS
    PT_TERMINATED_INVALID_STARTVALUE
    PT_TERMINATED_STEP_SIZE_TOO_SMALL
    PT_TERMINATED_ACCURACY_LIMIT
    PT_TERMINATED_ILL_CONDITIONED
    PT_POST_CHECK_FAILED
    PT_EXCESS_SOLUTION
end

"""
    path_tracker_status(code::CoreTrackerStatus)
Construct a [`PathTrackerStatus`](@ref) from a [`CoreTrackerStatus`](@ref).
"""
function path_tracker_status(code::CoreTrackerStatus)
    if code == CT_SUCCESS
        return PT_SUCCESS
    elseif code == CT_TERMINATED_INVALID_STARTVALUE
        return PT_TERMINATED_INVALID_STARTVALUE
    elseif code == CT_TERMINATED_MAX_ITERS
        return PT_TERMINATED_MAX_ITERS
    elseif code == CT_TERMINATED_STEP_SIZE_TOO_SMALL
        return PT_TERMINATED_STEP_SIZE_TOO_SMALL
    elseif code == CT_TERMINATED_ILL_CONDITIONED
        return PT_TERMINATED_ILL_CONDITIONED
    else
        return PT_TRACKING
    end
end

is_success(status::PathTrackerStatus) = status == PT_SUCCESS
is_at_infinity(status::PathTrackerStatus) = status == PT_AT_INFINITY
is_tracking(status::PathTrackerStatus) = status == PT_TRACKING


function symbol(status::PathTrackerStatus)
    if status == PT_SUCCESS
        :success
    elseif status == PT_AT_INFINITY
        :at_infinity
    elseif status == PT_TERMINATED_CALLBACK
        :terminated_callback
    elseif status == PT_TERMINATED_MAX_ITERS
        :terminated_max_iters
    elseif status == PT_TERMINATED_INVALID_STARTVALUE
        :terminated_invalid_startvalue
    elseif status == PT_TERMINATED_STEP_SIZE_TOO_SMALL
        :terminated_step_size_too_small
    elseif status == PT_TERMINATED_ACCURACY_LIMIT
        :terminated_accuracy_limit
    elseif status == PT_TERMINATED_ILL_CONDITIONED
        :terminated_ill_conditioned
    elseif status == PT_POST_CHECK_FAILED
        :post_check_failed
    elseif status == PT_EXCESS_SOLUTION
        :excess_solution
    else
        :unknown
    end
end


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
    status::PathTrackerStatus
    valuation::Valuation
    prediction::AV
    solution::AV
    solution_accuracy::Float64
    solution_cond::Float64
    solution_residual::Float64
    winding_number::Union{Nothing,Int}
    s::Float64
    eg_started::Bool
end

function PathTrackerState(x::AbstractVector; at_infinity_check::Bool = true)
    status = PT_TRACKING
    valuation = Valuation(x; affine = at_infinity_check)
    prediction = copy(x)
    solution = copy(x)
    solution_accuracy = NaN
    solution_cond = NaN
    solution_residual = NaN
    winding_number = nothing
    s = 0.0
    eg_started = false
    PathTrackerState(
        status,
        valuation,
        prediction,
        solution,
        solution_accuracy,
        solution_cond,
        solution_residual,
        winding_number,
        s,
        eg_started,
    )
end

Base.show(io::IO, S::PathTrackerState) = print_fieldnames(io, S)
Base.show(io::IO, ::MIME"application/prs.juno.inline", S::PathTrackerState) = S

function init!(state::PathTrackerState, s::Float64)
    state.status = PT_TRACKING
    init!(state.valuation)
    state.prediction .= 0.0
    state.solution .= 0.0
    state.solution_accuracy = NaN
    state.solution_cond = NaN
    state.solution_residual = NaN
    state.winding_number = nothing
    state.s = s
    state.eg_started = false
end


##################
## PATH TRACKER ##
##################

# The PathTracker sits on top of the CoreTracker
# it combines the path tracking (using a predictor-corrector) scheme with an endgame
# routine. It also implements additional logic to handle numerical different solutions.
# It is more opinionated than the CoreTracker and always compute in a logarithmic time scale.
#
# The PathTracker works in 2 stages:
# 1) Pre-Endgame
# 2) Endgame
#
# 1) Pre-Endgame
# Here we track from 0 to s_eg, the value where we start the endgame. In this phase a path
# is never considered going to infinity, nor that the Cauchy endgame already starts.
# Instead the focus is on getting solutions at s_eg.
# If the path tracker is hitting an accuracy limit it decreases (and increases again)
# the desired accuracy until in a prescribed range.
# If this is still not sufficient it passes to 128 the bit precision option.
# After arriving at s_eg a callback is called to check whether tracking should be continued.
#
# 2) Endgame
# Here we track from s_eg to s_t1. Now we cut paths off if we decide that they are diverging.
# If there is some sign of ill-conditioning (decreased limit accuracy or increasing condition number)
# and the valuation indicates that we are in the eg zone, then we start the Cauchy endgame.
# If we arrive at s_t1 we check that we also have a solution at s=Inf.
# If this is not the case we extend the tracking to s_t2 and check there again.


struct PathTracker{
    AV<:AbstractVector{Complex{Float64}},
    Prob<:AbstractProblem,
    CT<:CoreTracker{Float64,AV},
    F1<:Function,
}
    problem::Prob
    core_tracker::CT
    endgame::CauchyEndgame{AV}
    state::PathTrackerState{AV}
    options::PathTrackerOptions
    eg_start_callback::F1
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
    samples_per_loop::Int = 8,
    max_winding_number::Int = 12,
    endgame_start::Float64 = 2.0,
    endgame_start_callback = always_true,
    precision_strategy::PrecisionStrategy = PREC_STRATEGY_FINITE,
    min_accuracy::Float64 = max(core_tracker.options.accuracy, 1e-5),
    min_step_size_before_eg::Float64 = exp2(-40),
    min_step_size_eg = exp2(-100),
    s_always_consider_valuation::Float64 = -log(1e-16),
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
        min_step_size_before_eg,
        min_step_size_eg,
        precision_strategy,
        s_always_consider_valuation,
    )

    PathTracker(
        prob,
        core_tracker,
        endgame,
        state,
        options,
        endgame_start_callback,
        default_ct_options,
    )
end

default_at_infinity_check(prob::Problem{AffineTracking}) = true
default_at_infinity_check(prob::Problem{ProjectiveTracking}) = homvars(prob) !== nothing

always_true(_) = true

Base.show(io::IO, ::MIME"application/prs.juno.inline", PT::PathTracker) = PT
function Base.show(io::IO, S::PathTracker{AV}) where {AV}
    println(io, "PathTracker with solution type $AV")
end


function init!(tracker::PathTracker, x)
    options, ct_options = tracker.options, tracker.core_tracker.options
    copy!(ct_options, tracker.default_ct_options)
    init!(tracker.core_tracker, x, 0.0, options.endgame_start)
    init!(tracker.state, 0.0)

    if options.precision_strategy == PREC_STRATEGY_ALWAYS ||
       options.precision_strategy == PREC_STRATEGY_FINITE
        ct_options.precision = PRECISION_ADAPTIVE
    end
    # before endgame don't track the condition number
    # and we only update the limiting accuracy etc after a step failed
    ct_options.track_cond = false
    ct_options.steps_jacobian_info_update = -1
    ct_options.min_step_size = options.min_step_size_before_eg

    tracker
end

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
        cond_bad = cond(core_tracker) > 1e4
        near_accuracy_limit = 1e4 * core_tracker.state.limit_accuracy > core_tracker.options.accuracy
        consider_always = s > options.s_always_consider_valuation
        if !(consider_always || cond_bad || near_accuracy_limit)
            return state.status
        end

        # Judge the current valuation to determine how to proceed
        verdict = judge(state.valuation; tol = 1e-3, tol_at_infinity = 1e-4)
        if verdict == VAL_AT_INFINITY &&
           (consider_always || cond_bad) && options.at_infinity_check
            state.status = PT_AT_INFINITY
        # If we expect a finite value and there is some ill conditioning let's do the
        # Cauchy endgame
        elseif verdict == VAL_FINITE && (cond_bad || near_accuracy_limit)
            # Perform endgame to estimate singular solution
            @label run_cauchy_eg
            retcode, m, p_accuracy = predict!(state.prediction, core_tracker, endgame)

            if retcode == CAUCHY_SUCCESS
                state.winding_number = m
                state.status = PT_SUCCESS
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

            elseif retcode == CAUCHY_TERMINATED_ACCURACY_LIMIT
                if options.precision_strategy == PREC_STRATEGY_NEVER
                    state.status = PT_TERMINATED_ACCURACY_LIMIT
                # If we are here, we have a precision strategy which allows to use
                # higher precision. Therefore switch to it if this not yet happened
                # and retry to cauchy eg
                elseif core_tracker.options.precision == PRECISION_FIXED_64
                    core_tracker.options.precision = PRECISION_ADAPTIVE
                    @goto run_cauchy_eg
                end

            elseif retcode == CAUCHY_TERMINATED_MAX_WINDING_NUMBER
                # If we hit the max winding number, this is either due to the fact
                # that we are not yet in the eg zone or that the solution is just too singular
                # For now we just terminate here
                state.status = PT_TERMINATED_ILL_CONDITIONED
            else
                state.status = PT_TERMINATED_ILL_CONDITIONED
            end
        end

    # We split the tracking in two parts. Pre endgame and endgame. Therefore we have to
    # handle the success case twice
    elseif is_success(ct_status) && !state.eg_started
        # make callback, returns true if we can continue with the tracking
        if tracker.eg_start_callback(current_x(core_tracker))
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
        else
            state.status = PT_TERMINATED_CALLBACK
        end

    elseif is_success(ct_status)
        converged, s_acc, s_cond, s_res = check_converged!(
            state.solution,
            core_tracker,
            current_x(core_tracker),
            Inf;
            tol = tracker.default_ct_options.accuracy,
        )
        if converged
            state.status = PT_SUCCESS
            state.solution_accuracy = s_acc
            state.solution_cond = s_cond
            state.solution_residual = s_res
        else
            state.status = PT_POST_CHECK_FAILED
        end

    elseif ct_status == CT_TERMINATED_ACCURACY_LIMIT
        # First update the valuation  and check whether we are good
        update!(state.valuation, core_tracker)
        verdict = judge(state.valuation; tol = 1e-2, tol_at_infinity = 1e-4)
        # We also terminate at this point if the valuation indicates that we are substantially
        # less than 0
        if verdict == VAL_AT_INFINITY || minimum(state.valuation.ν) < -1
            state.status = PT_AT_INFINITY
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
            core_tracker.state.status = CT_TRACKING
        # 2)
        elseif options.precision_strategy == PREC_STRATEGY_FINITE
            # 2.a) We say that a path could still become finite if it is not classified
            #      as going to infinity for a loose tolerance
            #      Also if the current valuation is less than -1 (so significantly less than 0)
            #      we do not continue the tracking
            if verdict == VAL_FINITE ||
               judge(
                state.valuation;
                tol = 1e-2,
                tol_at_infinity = 1e-2,
            ) != VAL_AT_INFINITY
                core_tracker.options.precision = PRECISION_ADAPTIVE
                core_tracker.state.status = CT_TRACKING
            # 2.b)
            else
                state.status = PT_TERMINATED_ACCURACY_LIMIT
            end
        # 3)
        else
            state.status = PT_TERMINATED_ACCURACY_LIMIT
        end
    else
        state.status = path_tracker_status(ct_status)
    end

    state.status
end

"Update the valuation with the current state."
update!(val::Valuation, T::CoreTracker) =
    update!(val, T.state.x, T.state.ẋ, real(current_t(T)), T.predictor)

"""
    sign_of_ill_conditioning(CT::CoreTracker; safety_factor::Float64 = 1e4)

Returns `true` if the path has some signs of ill-conditioning. For this we check whether
the condition number of the Jacobian is larger than `safety_factor` or if
the safety factor to reaching limit accuracy is less than `safety_factor`.
"""
function sign_of_ill_conditioning(
    CT::CoreTracker;
    min_cond = 1e3,
    accuracy_safety_factor::Float64 = 1e5,
)
    # TODO: This yields wrong results if we track with adaptive precision
    # and have a too high accuracy.
    cond(CT) > min_cond || accuracy_safety_factor * CT.state.limit_accuracy > CT.options.accuracy
end

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

function track!(tracker::PathTracker, x)
    init!(tracker, x)
    while is_tracking(tracker.state.status)
        step!(tracker)
    end
    tracker.state.status
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

Obtain the solution computed by the `PathTracker`. This **does not** create a copy.
"""
solution(tracker::PathTracker) = solution(tracker.state)
solution(state::PathTrackerState) = state.solution


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
    tracker_startsolutions(prob, startsolutions; rest...)
end

function tracker_startsolutions(prob::Problem, startsolutions; kwargs...)
    tracker = PathTracker(prob, start_solution_sample(startsolutions); kwargs...)
    (tracker = tracker, startsolutions = startsolutions)
end

"""
    pathtracker(args...; kwargs...)

Construct a [`PathTracker`](@ref) in the same way [`solve`](@ref) does it.
This also takes the same input arguments as `solve` with the exception that you do not need to specify startsolutions.

## Examples

### Obtain single solution
We want to construct a path tracker to track a parameterized system `f` with parameters `p`
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

A `PathResult` is the result of tracking of a path using [`PathTracker`](@ref).
Its fields are

* `return_code`: One of `:success`, `:at_infinity` or any error code in [`PathTrackerStatus.states`](@ref) converted to a `Symbol`.
* `solution::V`: The solution vector.
* `t::Float64`: The value of `t` at which `solution` was computed. Note that if `return_code` is `:at_infinity`, then `t` is the value when this was decided.
* `accuracy::Union{Nothing, Float64}`: An approximation of ``||x-x^*||₂`` where ``x`` is the computed solution and ``x^*`` is the true solution.
* `residual::Union{Nothing, Float64}`: The value of the 2-norm of `H(solution, 0)`.
* `multiplicity::Union{Nothing, Int}` is the multiplicity of the `solution`. This is only assigned by. [`singular`](@ref).
* `condition_jacobian::Union{Nothing, Float64}`: This is the condition number of the row-equilibrated Jacobian at the solution. A high condition number indicates a singularity.
* `winding_number:Union{Nothing, Int}`: The estimated winding number. This is a lower bound on the multiplicity of the solution.
* `start_solution::Union{Nothing, Int}`: The start solution of the path.
* `accepted_steps::Int`: The number of accepted steps during the path tracking.
* `rejected_steps::Int`: The number of rejected steps during the path tracking.
* `valuation::Union{Nothing, Vector{Float64}}`: An approximation of the valuation of the Puiseux series expansion of `x(t)`.
* `valuation_accuracy::Union{Nothing, Vector{Float64}}`: An estimate of the accuracy of the valuation of the Puiseux series expansion of `x(t)`.

     PathResult(tracker::PathTracker, start_solution=nothing, path_number::Union{Nothing,Int}=nothing; details=:default)

Possible `details` values are `:minimal` (minimal details), `:default` (default) and `:extensive` (all information possible).
"""
struct PathResult{V<:AbstractVector}
    return_code::Symbol
    solution::V
    t::Float64
    accuracy::Union{Nothing,Float64}
    residual::Union{Nothing,Float64} # level 1+
    condition_jacobian::Union{Nothing,Float64}
    winding_number::Union{Nothing,Int}
    path_number::Union{Nothing,Int}
    start_solution::Union{Nothing,V} # level 1+
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
    return_code = symbol(state.status)
    windingnumber = winding_number(tracker)
    x = copy(solution(tracker))

    # accuracy
    if !isnan(state.solution_accuracy)
        accuracy = state.solution_accuracy
    else
        accuracy = nothing
    end

    if return_code == :success
        t = 0.0
    else
        t = exp(-state.s)
    end
    # condition
    if isnan(state.solution_cond)
        condition_jac = cond(core_tracker)
    else
        condition_jac = state.solution_cond
    end
    # residual
    if return_code == :success && details_level ≥ 1
        res = residual(tracker)
    else
        res = nothing
    end

    if details_level ≥ 1
        # mimic the behaviour in track! to get a start solution of the same type as x
        embed!(cache.base_point, tracker.problem, start_solution)
        startsolution = pull_back(tracker.problem, cache.base_point; regauge = false)
    else
        startsolution = nothing
    end

    accepted_steps = core_tracker.state.accepted_steps
    rejected_steps = core_tracker.state.rejected_steps

    if details_level == 2
        valuation = copy(tracker.state.val.v)
        valuation_accuracy = abs.(tracker.state.val.v̇)
    else
        valuation = nothing
        valuation_accuracy = nothing
    end

    # this is only assigned by using the singular() function
    multiplicity = nothing

    PathResult(
        return_code,
        x,
        t,
        accuracy,
        res,
        multiplicity,
        condition_jac,
        windingnumber,
        path_number,
        startsolution,
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

function residual(tracker::PathTracker, x = tracker.state.solution)
    evaluate!(
        tracker.cache.target_residual,
        tracker.cache.target_system,
        x,
        tracker.cache.target_newton_cache.system_cache,
    )
    euclidean_norm(tracker.cache.target_residual)
end

function condition_jacobian(tracker::PathTracker, x = tracker.state.solution)
    jac = tracker.cache.target_jacobian
    jacobian!(
        jac.J,
        tracker.cache.target_system,
        x,
        tracker.cache.target_newton_cache.system_cache,
    )
    updated_jacobian!(jac; update_infos = true)
    jac.cond
end

"""
    multiplicity(P::PathResult{T})

Returns the multiplicity of `P`.
"""
multiplicity(P::PathResult{T}) where {T} = P.multiplicity

"""
    result_type(tracker::PathTracker)

Returns the type of result `track` will return.
"""
result_type(tracker::PathTracker) = PathResult{typeof(solution(tracker))}

function Base.show(io::IO, r::PathResult)
    iscompact = get(io, :compact, false)
    if iscompact || haskey(io, :typeinfo)
        println(io, "• return_code: $(r.return_code)")
        if r.return_code != PathTrackerStatus.success
            println(io, " • t: $(r.t)")
        end
        println(io, " • solution: ", r.solution)
        r.accuracy !== nothing && println(
            io,
            " • accuracy: $(Printf.@sprintf "%.3e" r.accuracy)",
        )
        r.winding_number !== nothing && println(
            io,
            " • winding_number: $(r.winding_number)",
        )
        r.path_number !== nothing && println(io, " • path_number: ", r.path_number)
    else
        println(io, "PathResult")
        println(io, "=================")
        println(io, " • return_code: $(r.return_code)")
        if r.return_code != PathTrackerStatus.success
            println(io, " • t: $(r.t)")
        end
        println(io, " • solution: ", r.solution)
        r.accuracy !== nothing && println(
            io,
            " • accuracy: $(Printf.@sprintf "%.3e" r.accuracy)",
        )
        r.residual !== nothing && println(
            io,
            " • residual: $(Printf.@sprintf "%.3e" r.residual)",
        )
        r.winding_number !== nothing && println(
            io,
            " • winding_number: $(r.winding_number)",
        )
        r.condition_jacobian !== nothing && println(
            io,
            " • condition_jacobian: $(Printf.@sprintf "%.3e" r.condition_jacobian)",
        )
        r.path_number !== nothing && println(io, " • path_number: ", r.path_number)
    end
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::PathResult) = x

"""
    solution(pathresult)

Get the solution of the path.
"""
solution(r::PathResult) = r.solution


"""
    accuracy(pathresult)

Get the accuracy of the solution ``x`` of the path, i.e., ``||H(x, 0)||₂``.
"""
accuracy(r::PathResult) = r.accuracy


"""
    residual(pathresult)

Get the residual of the solution ``x`` of the path, i.e., ``||H(x, 0)||₂``.
"""
residual(r::PathResult) = r.residual

"""
    start_solution(pathresult)

Get the start solution of the solution ``x`` of the path.
"""
start_solution(r::PathResult) = r.start_solution

"""
    is_success(pathresult)

Checks whether the path is successfull.
"""
is_success(r::PathResult) = r.return_code == :success

"""
    is_failed(pathresult)

Checks whether the path failed.
"""
is_failed(r::PathResult) = !(r.return_code == :at_infinity || r.return_code == :success)


"""
    is_at_infinity(pathresult)

Checks whether the path goes to infinity.
"""
is_at_infinity(r::PathResult) = r.return_code == :at_infinity

"""
    isfinite(pathresult)

Checks whether the path result is finite.
"""
Base.isfinite(r::PathResult) = r.return_code == :success # we don't check is_affine to make other code easier

"""
    is_singular(pathresult; tol=1e10)

Checks whether the path result is singular. This is true if
the multiplicity is larger than  1 or if the condition number of the Jacobian
is larger than `tol`.
"""
is_singular(r::PathResult; tol = 1e10) = is_singular(r, tol)
function is_singular(r::PathResult, tol::Real)
    (unpack(r.condition_jacobian, 1.0) > tol || unpack(r.multiplicity, 1) > 1) &&
    is_success(r)
end

"""
    is_nonsingular(pathresult; tol=1e10)

Checks whether the path result is non-singular. This is true if
it is not singular.
"""
is_nonsingular(r::PathResult; kwargs...) = !is_singular(r; kwargs...) && is_success(r)
is_nonsingular(r::PathResult, tol::Real) = !is_singular(r, tol) && is_success(r)


"""
    is_real(pathresult; tol=1e-6)

We consider a result as `real` if the 2-norm of the imaginary part of the solution is at most `tol`.
"""
is_real(r::PathResult; tol = 1e-6) = is_real(r, tol)
is_real(r::PathResult, tol::Real) = is_real_vector(r.solution, tol)
# provide fallback since this in in Base
Base.isreal(r::PathResult, tol) = is_real(r, tol)
Base.isreal(r::PathResult; kwargs...) = is_real(r; kwargs...)

"""
    is_projective(pathresult)

Return`s true if the solution is a projective vector.
"""
is_projective(r::PathResult{<:PVector}) = true
is_projective(r::PathResult) = false

"""
    is_affine(pathresult)

Return`s true if the solution is an affine vector.
"""
is_affine(r::PathResult) = !is_projective(r)
