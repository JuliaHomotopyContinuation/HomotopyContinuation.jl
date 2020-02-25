Base.@kwdef mutable struct EndgameTrackerOptions
    endgame_start::Float64 = 0.1
    min_cond_at_infinity::Float64 = 1e4
    max_winding_number::Int = 12
end


module EGTrackerReturnCode
import ..TrackerReturnCode
"""
    enum EGTrackerReturnCode

The possible states a `PathTracker` can be in:

* `EGTrackerReturnCode.tracking`
* `EGTrackerReturnCode.success`
* `EGTrackerReturnCode.at_infinity`
* `EGTrackerReturnCode.excess_solution`
* `EGTrackerReturnCode.post_check_failed`
* `EGTrackerReturnCode.terminated_accuracy_limit`
* `EGTrackerReturnCode.terminated_ill_conditioned`
* `EGTrackerReturnCode.terminated_invalid_startvalue`
* `EGTrackerReturnCode.terminated_max_winding_number`
* `EGTrackerReturnCode.terminated_max_iters`
* `EGTrackerReturnCode.terminated_step_size_too_small`
"""
@enum codes begin
    tracking
    success
    at_infinity
    terminated_accuracy_limit
    terminated_invalid_startvalue
    terminated_ill_conditioned
    terminated_max_iters
    terminated_max_winding_number
    terminated_step_size_too_small
    terminated_unknown
    post_check_failed
    excess_solution
end

function Base.convert(::Type{codes}, code::TrackerReturnCode.codes)
    if code == TrackerReturnCode.success
        return success
    elseif code == TrackerReturnCode.terminated_max_iters
        return terminated_max_iters
    elseif code == TrackerReturnCode.terminated_accuracy_limit
        return terminated_accuracy_limit
    elseif code == TrackerReturnCode.terminated_ill_conditioned
        return terminated_ill_conditioned
    elseif code == TrackerReturnCode.terminated_invalid_startvalue
        return terminated_invalid_startvalue
    elseif code == TrackerReturnCode.terminated_step_size_too_small
        return terminated_step_size_too_small
    elseif code == TrackerReturnCode.terminated_unknown
        return terminated_unknown
    else
        return tracking
    end
end

end

# """
#     endgame_tracker_code(code::TrackerReturnCode.codes)
#
# Construct a [`EGTrackerReturnCode.codes`](@ref) from a [`TrackerReturnCode.codes`](@ref).
# """
# endgame_tracker_code(code::TrackerReturnCode.codes) =
#     EGTrackerReturnCode.code(code)

"""
    is_success(code::EGTrackerReturnCode.codes)

Returns `true` if `status` indicates a success in tracking.
"""
is_success(code::EGTrackerReturnCode.codes) =
    code == EGTrackerReturnCode.success

"""
    is_success(code::EGTrackerReturnCode.codes)

Returns `true` if `status` indicates that a path diverged towards infinity.
"""
is_at_infinity(code::EGTrackerReturnCode.codes) =
    code == EGTrackerReturnCode.at_infinity

"""
    is_tracking(code::EGTrackerReturnCode.codes)

Returns `true` if `status` indicates the tracking is not going on.
"""
is_tracking(code::EGTrackerReturnCode.codes) =
    code == EGTrackerReturnCode.tracking

"""
    is_invalid_startvalue(code::EGTrackerReturnCode.codes)

Returns `true` if the provided start value was not valid.
"""
is_invalid_startvalue(code::EGTrackerReturnCode.codes) =
    code == EGTrackerReturnCode.terminated_invalid_startvalue

"""
    is_terminated_callback(code::EGTrackerReturnCode.codes)

Returns `true` if the provided callback indicated a termination of the path.
"""
is_terminated_callback(code::EGTrackerReturnCode.codes) =
    code == EGTrackerReturnCode.terminated_callback


# State

mutable struct EndgameTrackerState{V<:AbstractVector}
    code::EGTrackerReturnCode.codes
    val::Valuation
    solution::V
    winding_number::Union{Nothing,Int}
    accuracy::Float64
    prediction::V
    last_point::V
    last_t::Float64
end

function EndgameTrackerState(x::AbstractVector)
    code = EGTrackerReturnCode.tracking
    val = Valuation(length(x))
    solution = copy(x)
    winding_number = nothing
    accuracy = NaN
    prediction = copy(x)
    last_point = copy(x)
    last_t = NaN
    EndgameTrackerState(
        code,
        val,
        solution,
        winding_number,
        accuracy,
        prediction,
        last_point,
        last_t,
    )
end


struct EndgameTracker{
    H<:AbstractHomotopy,
    # V and V̄ need to have the same container type
    V<:AbstractVector{ComplexF64},
    V̄<:AbstractVector{ComplexDF64},
    M<:AbstractMatrix{ComplexF64},
}
    tracker::Tracker{H,V,V̄,M}
    state::EndgameTrackerState
    options::EndgameTrackerOptions
end

function EndgameTracker(tracker::Tracker; kwargs...)
    options = EndgameTrackerOptions(; kwargs...)
    state = EndgameTrackerState(tracker.state.x)
    EndgameTracker(tracker, state, options)
end

function init!(eg_tracker::EndgameTracker, x, t₁::Real)
    @unpack tracker, state, options = eg_tracker

    init!(tracker, x, t₁, 0.0)

    state.code = status(tracker)
    init!(state.val)
    state.solution .= NaN
    state.accuracy = NaN
    state.winding_number = nothing

    eg_tracker
end

"""
    CauchyEndgameResult

An enum indicating the result of the [`cauchy!`](@ref) computation.

# Cases
* `CAUCHY_SUCCESS`: The endgame was successfull.
* `CAUCHY_TERMINATED_MAX_WINDING_NUMBER`: The endgame was terminated since the winding
  number is larger than the provided threshold.
* `CAUCHY_TERMINATED`: The endgame was terminated due to some other error in the path
  tracking.
"""
@enum CauchyEndgameResult begin
    CAUCHY_SUCCESS
    CAUCHY_TERMINATED_MAX_WINDING_NUMBER
    CAUCHY_TERMINATED
end

"""
    cauchy!(p, tracker::CoreTracker, ::CauchyEndgame)

Try to predict the value of `x(0)` using the [`CauchyEndgame`](@ref).
For this we track the polygon defined by ``te^{i2πk/n}`` until we end again at ``x``.
Here ``n`` is the number of samples we take per loop, `samples_per_loop`.
The computation gives up if we have a winding number larger than `max_winding_number`.
It returns a tuple denoting the success ([`CauchyEndgameResult`](@ref)) the computed
winding number `m::Int` and th expected accuracy of the solution.

[Cauchy's integral formula]: https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula
"""
function cauchy!(
    state::EndgameTrackerState,
    tracker::Tracker,
    options::EndgameTrackerOptions,
)
    @unpack last_point, prediction = state

    t = real(tracker.state.t)
    n₀ = max(ceil(Int, log(eps()) / log(t)), 3)
    @unpack x, μ, ω = tracker.state

    last_point .= x
    state.last_t = t
    prediction .= 0.0

    m = 1
    Δθ = 2π / n₀
    point_acc = tracker.state.accuracy
    result = CAUCHY_TERMINATED_MAX_WINDING_NUMBER
    while m ≤ options.max_winding_number
        θⱼ = 0.0
        tⱼ = t
        for j = 1:n₀
            θⱼ += Δθ
            tⱼ = j == n₀ ? t : t * cis(θⱼ)
            res = track!(tracker, tⱼ)
            # TODO: Refine to guarantee high accuracy
            point_acc = max(point_acc, tracker.state.accuracy)

            if !is_success(res.code)
                result = CAUCHY_TERMINATED
                @goto _return
            end

            prediction .+= x
        end
        # Check that loop is closed
        if tracker.state.norm(last_point, x) <
           4 * max(point_acc, tracker.state.accuracy)
            n = n₀ * m
            prediction .= prediction ./ n

            result = CAUCHY_SUCCESS
            break
        end

        m += 1
    end

    @label _return

    if result == CAUCHY_TERMINATED
        init!(tracker, base_point, t, 0.0, ω, μ; keep_steps = true)
    else
        init!(tracker, 0.0)
    end

    result, m, point_acc
end


function step!(eg_tracker::EndgameTracker)
    @unpack tracker, state, options = eg_tracker

    step!(tracker)
    t = real(tracker.state.t)
    state.code = status(tracker)

    is_tracking(state.code) || return state.code
    t < options.endgame_start || return state.code
    !tracker.state.last_step_failed || return state.code

    # update valution
    @unpack x, x¹, x², x³, x⁴ = tracker.state
    update!(state.val, x, x¹, x², x³, x⁴, t)

    println(state.val)

    verdict = judge(state.val, 1e-2)
    @show t verdict LA.cond(tracker)

    if verdict == VAL_AT_INFINITY
        # && LA.cond(tracker) > options.min_cond_at_infinity
        return (state.code = EGTrackerReturnCode.at_infinity)
    end

    # if valuation indicates singular solution
    # or finite but bad conditioned, start cauchy endgame
    if verdict == VAL_FINITE_SINGULAR || (
        verdict == VAL_FINITE &&
        (LA.cond(tracker) > 1e6 || tracker.state.use_extended_prec)
    )
        res, m, acc_est = cauchy!(state, tracker, options)
        if res == CAUCHY_SUCCESS
            state.winding_number = m
            state.solution .= state.prediction
            state.accuracy = acc_est
            state.code = EGTrackerReturnCode.success
        end
    end

    state.code
end

function track!(eg_tracker::EndgameTracker, x, t₁::Real)
    init!(eg_tracker, x, t₁)

    while is_tracking(eg_tracker.state.code)
        step!(eg_tracker)
    end

    eg_tracker.state.code
end
