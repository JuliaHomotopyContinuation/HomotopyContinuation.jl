export PathTracker,
    PathResult,
    solution,
    accuracy,
    is_success,
    is_at_infinity,
    is_failed,
    is_finite,
    winding_number,
    last_path_point,
    steps,
    accepted_steps,
    rejected_steps,
    valuation

Base.@kwdef mutable struct PathTrackerOptions
    endgame_start::Float64 = 0.1
    # eg parameter
    min_cond_eg::Float64 = 1e6
    # valuation etc
    at_zero_check::Bool = false
    at_infinity_check::Bool = true
    val_finite_tol::Float64 = 1e-2
    val_at_infinity_tol::Float64 = 1e-3
    # singular solutions parameters
    max_winding_number::Int = 20
end

Base.show(io::IO, opts::PathTrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::PathTrackerOptions) = opts

module PathTrackerCode
import ..TrackerCode
"""
    enum PathTrackerCode

The possible states a `PathTracker` can be in:

* `PathTrackerCode.tracking`
* `PathTrackerCode.success`
* `PathTrackerCode.at_infinity`
* `PathTrackerCode.at_zero`
* `PathTrackerCode.excess_solution`
* `PathTrackerCode.post_check_failed`
* `PathTrackerCode.terminated_accuracy_limit`
* `PathTrackerCode.terminated_ill_conditioned`
* `PathTrackerCode.terminated_invalid_startvalue`
* `PathTrackerCode.terminated_max_winding_number`
* `PathTrackerCode.terminated_max_iters`
* `PathTrackerCode.terminated_step_size_too_small`
"""
@enum codes begin
    tracking
    success
    at_infinity
    at_zero
    terminated_accuracy_limit
    terminated_invalid_startvalue
    terminated_ill_conditioned
    terminated_max_iters
    terminated_max_winding_number
    terminated_step_size_too_small
    terminated_unknown
    post_check_failed
    excess_solution
    polyhedral_failed
end

function Base.convert(::Type{codes}, code::TrackerCode.codes)
    if code == TrackerCode.success
        return success
    elseif code == TrackerCode.terminated_max_iters
        return terminated_max_iters
    elseif code == TrackerCode.terminated_accuracy_limit
        return terminated_accuracy_limit
    elseif code == TrackerCode.terminated_ill_conditioned
        return terminated_ill_conditioned
    elseif code == TrackerCode.terminated_invalid_startvalue
        return terminated_invalid_startvalue
    elseif code == TrackerCode.terminated_step_size_too_small
        return terminated_step_size_too_small
    elseif code == TrackerCode.terminated_unknown
        return terminated_unknown
    else
        return tracking
    end
end

end


"""
    is_success(code::PathTrackerCode.codes)

Returns `true` if `status` indicates a success in tracking.
"""
is_success(code::PathTrackerCode.codes) = code == PathTrackerCode.success

"""
    is_success(code::PathTrackerCode.codes)

Returns `true` if `status` indicates that a path diverged towards infinity.
"""
is_at_infinity(code::PathTrackerCode.codes) =
    code == PathTrackerCode.at_infinity

"""
    is_tracking(code::PathTrackerCode.codes)

Returns `true` if `status` indicates the tracking is not going on.
"""
is_tracking(code::PathTrackerCode.codes) = code == PathTrackerCode.tracking

"""
    is_invalid_startvalue(code::PathTrackerCode.codes)

Returns `true` if the provided start value was not valid.
"""
is_invalid_startvalue(code::PathTrackerCode.codes) =
    code == PathTrackerCode.terminated_invalid_startvalue

"""
    is_terminated_callback(code::PathTrackerCode.codes)

Returns `true` if the provided callback indicated a termination of the path.
"""
is_terminated_callback(code::PathTrackerCode.codes) =
    code == PathTrackerCode.terminated_callback


# State
Base.@kwdef mutable struct PathTrackerState{V<:AbstractVector}
    code::PathTrackerCode.codes = PathTrackerCode.tracking
    val::Valuation
    endgame_started::Bool = false
    row_scaling::Vector{Float64}
    col_scaling::Vector{Float64}
    cond_eg_start::Float64 = 1.0
    solution::V
    winding_number::Union{Nothing,Int} = nothing
    accuracy::Float64 = NaN
    cond::Float64 = 1.0
    prediction::V
    last_point::V
    last_t::Float64 = NaN
    max_winding_number_hit::Bool = false
    cauchy_failures::Int = 0
    jump_to_zero_failed::Tuple{Bool,Bool} = (false, false)
end

PathTrackerState(npolynomials::Integer, x::AbstractVector) = PathTrackerState(;
    val = Valuation(length(x)),
    row_scaling = zeros(npolynomials),
    col_scaling = zeros(length(x)),
    solution = copy(x),
    prediction = copy(x),
    last_point = copy(x),
)


struct PathTracker{
    H<:AbstractHomotopy,
    N, # AutomaticDifferentiation
    # V and V̄ need to have the same container type
    V<:AbstractVector{ComplexF64},
    V̄<:AbstractVector{ComplexDF64},
    M<:AbstractMatrix{ComplexF64},
}
    tracker::Tracker{H,N,V,V̄,M}
    state::PathTrackerState{V}
    options::PathTrackerOptions
end

PathTracker(H::AbstractHomotopy; kwargs...) = PathTracker(Tracker(H); kwargs...)
function PathTracker(tracker::Tracker; kwargs...)
    options = PathTrackerOptions(; kwargs...)
    state = PathTrackerState(size(tracker.homotopy, 1), tracker.state.x)
    PathTracker(tracker, state, options)
end

Base.broadcastable(T::PathTracker) = Ref(T)

function init!(eg_tracker::PathTracker, x, t₁::Real; ω::Float64 = NaN, μ::Float64 = NaN)
    @unpack tracker, state, options = eg_tracker

    init!(tracker, x, t₁, 0.0; ω = ω, μ = μ)

    state.code = status(tracker)
    init!(state.val)
    state.endgame_started = false
    state.row_scaling .= 1
    state.col_scaling .= weights(tracker.state.norm)
    state.cond_eg_start = 1.0
    state.cond = 1.0
    state.solution .= NaN
    state.accuracy = NaN
    state.winding_number = nothing
    state.max_winding_number_hit = false
    state.cauchy_failures = 0
    state.jump_to_zero_failed = (false, false)

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
function cauchy!(state::PathTrackerState, tracker::Tracker, options::PathTrackerOptions)
    @unpack last_point, prediction = state

    t = real(tracker.state.t)
    n₀ = max(ceil(Int, log(eps()) / log(t)), 3)
    @unpack x, μ, ω = tracker.state

    point_acc = refine_current_solution!(tracker)
    state.last_point .= tracker.state.x
    state.last_t = tracker.state.t
    prediction .= 0.0

    m = 1
    Δθ = 2π / n₀
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
        d = tracker.state.norm(last_point, x)
        if d < 10 * max(point_acc, tracker.state.accuracy)
            n = n₀ * m
            prediction .= prediction ./ n

            result = CAUCHY_SUCCESS
            break
        end
        m += 1
    end

    @label _return

    if result == CAUCHY_TERMINATED
        init!(tracker, last_point, t, 0.0; ω = ω, μ = μ, keep_steps = true)
    else
        init!(tracker, 0.0)
    end

    result, m, point_acc
end

function step!(eg_tracker::PathTracker, debug::Bool = false)
    @unpack tracker, state, options = eg_tracker

    proposed_t′ = real(tracker.state.t′)
    state.last_point .= tracker.state.x
    state.last_t = tracker.state.t
    step!(tracker)
    # @show tracker.state.t

    state.code = status(tracker)

    if !is_tracking(state.code)
        state.accuracy = tracker.state.accuracy
        state.solution .= tracker.state.x
        state.winding_number = nothing
        state.cond = LA.cond(tracker.state.jacobian, state.row_scaling, state.col_scaling)

        # If a path got terminated, let's be more generous to classify them as at_infinity
        if state.code == PathTrackerCode.terminated_accuracy_limit ||
           state.code == PathTrackerCode.terminated_ill_conditioned

            verdict = analyze(
                state.val;
                finite_tol = 0.0,
                at_infinity_tol = sqrt(options.val_at_infinity_tol),
                zero_is_finite = !options.at_zero_check,
            )

            if options.at_infinity_check && verdict.at_infinity
                return (state.code = PathTrackerCode.at_infinity)
            elseif options.at_zero_check && verdict.at_zero
                return (state.code = PathTrackerCode.at_zero)
            end
        end
        return state.code
    end


    t = real(tracker.state.t)
    t < options.endgame_start || return state.code

    if !state.endgame_started
        #=
         We enter the endgame.
         We monitor the condition number with a fixed scaling of the columns and rows,
         which is recorded now.
        =#
        state.col_scaling .= weights(tracker.state.norm)
        row_scaling!(
            state.row_scaling,
            jacobian(tracker.state.jacobian),
            state.col_scaling,
        )
        κ = LA.cond(tracker.state.jacobian, state.row_scaling, state.col_scaling)
        state.cond_eg_start = κ

        tracker.state.use_strict_β_τ = true
        state.endgame_started = true
    end

    if tracker.state.last_step_failed
        state.jump_to_zero_failed = (last(state.jump_to_zero_failed), iszero(proposed_t′))
        return state.code
    else
        state.jump_to_zero_failed = (last(state.jump_to_zero_failed), false)
    end


    update!(state.val, tracker.state.tx³, t)
    (finite, at_infinity, at_zero) = analyze(
        state.val;
        finite_tol = options.val_finite_tol * exp10(-state.cauchy_failures),
        at_infinity_tol = options.val_at_infinity_tol,
        zero_is_finite = !options.at_zero_check,
    )

    if debug
        color = tracker.state.extended_prec ? :blue : :yellow
        printstyled("t = ", t, "  t′ = ", real(tracker.state.t′), "\n", color = color)
        κ = LA.cond(tracker.state.jacobian, state.row_scaling, state.col_scaling)
        print("κ = ")
        Printf.@printf("%.3e (%.3e)\n", κ, κ / state.cond_eg_start)
        println(state.val)

        printstyled("Val: ", bold = true)
        if finite
            printstyled("finite ", color = :blue, bold = true)
        end
        if at_infinity
            printstyled("at_infinity ", color = :blue, bold = true)
        end
        if at_zero
            printstyled("at_zero ", color = :blue, bold = true)
        end
        println("\n")
    end

    if at_infinity || finite || at_zero
        state.cond = LA.cond(tracker.state.jacobian, state.row_scaling, state.col_scaling)
    end

    if options.at_infinity_check && at_infinity && state.cond > options.min_cond_eg
        state.accuracy = tracker.state.accuracy
        state.solution .= tracker.state.x
        return (state.code = PathTrackerCode.at_infinity)

    elseif options.at_zero_check && at_zero && state.cond > options.min_cond_eg
        state.accuracy = tracker.state.accuracy
        state.solution .= tracker.state.x
        return (state.code = PathTrackerCode.at_zero)

    elseif finite && state.cond > options.min_cond_eg
        res, m, acc_est = cauchy!(state, tracker, options)
        if debug
            printstyled("Cauchy result: ", res, " ", m, " ", acc_est, "\n"; color = :blue)
        end
        if res == CAUCHY_SUCCESS
            if state.winding_number === nothing
                @label save_cauchy_result
                state.winding_number = m
                state.solution .= state.prediction
                state.accuracy = acc_est
            elseif state.winding_number != m
                state.cauchy_failures += 1
                @goto save_cauchy_result
            elseif state.winding_number == m
                d = tracker.state.norm(state.prediction, state.solution)
                if d < 100 * max(state.accuracy, acc_est)
                    state.solution .= state.prediction
                    state.accuracy = d
                    state.cond = LA.cond(
                        tracker,
                        state.solution,
                        0.0,
                        state.row_scaling,
                        state.col_scaling,
                    )
                    return (state.code = PathTrackerCode.success)
                else
                    state.cauchy_failures += 1
                    @goto save_cauchy_result
                end
            end
        elseif res == CAUCHY_TERMINATED_MAX_WINDING_NUMBER
            if state.max_winding_number_hit
                state.code = PathTrackerCode.terminated_max_winding_number
            else
                state.cauchy_failures += 1
                state.max_winding_number_hit = true
            end
        elseif res == CAUCHY_TERMINATED
            state.code = tracker.state.code
        end
    end

    state.code
end

function track!(
    eg_tracker::PathTracker,
    x,
    t₁::Real;
    ω::Float64 = NaN,
    μ::Float64 = NaN,
    debug::Bool = false,
)
    init!(eg_tracker, x, t₁; ω = ω, μ = μ)

    while is_tracking(eg_tracker.state.code)
        step!(eg_tracker, debug)
    end

    eg_tracker.state.code
end

Base.@kwdef mutable struct PathResult{V<:AbstractVector}
    return_code::PathTrackerCode.codes
    solution::V
    t::Float64
    accuracy::Float64
    cond::Float64
    winding_number::Union{Nothing,Int}
    last_path_point::Union{Nothing,Tuple{V,Float64}}
    valuation::Union{Nothing,Vector{Float64}}
    ω::Float64
    μ::Float64
    extended_precision::Bool
    # performance stats
    accepted_steps::Int
    rejected_steps::Int
    extended_precision_used::Bool
end

function PathResult(egtracker::PathTracker)
    @unpack tracker, state, options = egtracker
    t = real(tracker.state.t)
    PathResult(
        return_code = state.code,
        solution = is_success(state.code) ? copy(state.solution) : copy(tracker.state.x),
        t = is_success(state.code) ? 0.0 : t,
        accuracy = state.accuracy,
        cond = state.cond,
        winding_number = state.winding_number,
        last_path_point = (copy(state.last_point), state.last_t),
        valuation = t > options.endgame_start ? nothing : copy(state.val.val_x),
        ω = tracker.state.ω,
        μ = tracker.state.μ,
        extended_precision = tracker.state.extended_prec,
        accepted_steps = tracker.state.accepted_steps,
        rejected_steps = tracker.state.rejected_steps,
        extended_precision_used = tracker.state.used_extended_prec,
    )
end

Base.show(io::IO, r::PathResult) = print_fieldnames(io, r)
Base.show(io::IO, ::MIME"application/prs.juno.inline", r::PathResult) = r


function track(eg_tracker::PathTracker, x, t₁::Real = 1.0; kwargs...)
    track!(eg_tracker, x, t₁; kwargs...)
    PathResult(eg_tracker)
end


"""
    solution(r::PathResult)

Get the solution of the path.
"""
solution(r::PathResult) = r.solution


"""
    accuracy(r::PathResult)

Get the accuracy of the solution. This is an estimate of the (relative) distance to the
true solution.
"""
accuracy(r::PathResult) = r.accuracy

"""
    steps(r::PathResult)

Total number of steps the path tracker performed.
"""
steps(r::PathResult) = accepted_steps(r) + rejected_steps(r)

"""
    accepted_steps(r::PathResult)

Total number of steps the path tracker accepted.
"""
accepted_steps(r::PathResult) = r.accepted_steps


"""
    rejected_steps(r::PathResult)

Total number of steps the path tracker rejected.
"""
rejected_steps(r::PathResult) = r.rejected_steps


"""
    winding_number(r::PathResult)

Get the winding number of the solution of the path. Returns `nothing` if it wasn't computed.
"""
winding_number(r::PathResult) = r.winding_number

"""
    last_path_point(r::PathResult)

Returns a tuple `(x,t)` containing the last zero of `H(x, t)` before the Cauchy endgame was used.
Returns `nothing` if the endgame strategy was not invoked.
"""
last_path_point(r::PathResult) = r.last_path_point

"""
    is_success(r::PathResult)

Checks whether the path is successfull.
"""
is_success(r::PathResult) = is_success(r.return_code)

"""
    is_failed(r::PathResult)

Checks whether the path failed.
"""
is_failed(r::PathResult) = !(is_at_infinity(r) || is_success(r))

"""
    is_at_infinity(r::PathResult)

Checks whether the path goes to infinity.
"""
is_at_infinity(r::PathResult) = is_at_infinity(r.return_code)


"""
    is_finite(r::PathResult)

Checks whether the path result is finite.
"""
is_finite(r::PathResult) = is_success(r.return_code)

"""
    valuation(r::PathResult)

Get the computed valuation of the path.
"""
valuation(r::PathResult) = r.valuation
