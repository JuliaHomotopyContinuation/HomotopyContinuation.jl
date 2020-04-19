export PathTracker,
    PathTrackerOptions,
    PathResult,
    track,
    solution,
    accuracy,
    residual,
    steps,
    accepted_steps,
    rejected_steps,
    winding_number,
    path_number,
    start_solution,
    multiplicity,
    last_path_point,
    valuation,
    is_success,
    is_at_infinity,
    is_excess_solution,
    is_failed,
    is_finite,
    is_singular,
    is_nonsingular,
    is_real



"""
    PathTrackerOptions(; options...)

Options controlling the behaviour of a [`PathTracker`](@ref).

## Options

* `endgame_start::Float64 = 0.1`: The point `t` in time where the endgame starts.

### Endgame parameters
These parameters control the behaviour during the endgame. See [^BT20] for details.

* `min_cond_eg::Float64 = 1e6`: The minimal condition number for which an endgame strategy
  is applied.
* `at_infinity_check::Bool = true`: Whether divering paths should be truncated.
* `zero_is_at_infinity::Bool = false`: Whether paths going to a solution where at least one
  coordinates is zero should also be considered diverging.
* `val_finite_tol::Float64 = 1e-2`: Tolerance on the valuation which has to be satisfied
  before the Cauchy endgame is started.
* `max_winding_number::Int = 20`: The maximal winding number which is attempted in the
  Cauchy endgame.
* `val_at_infinity_tol::Float64 = 1e-3`: Tolerance on the valuation which has to be
  satisfied before a path is considered to diverge / go to infinity.

[^BT20]: Breiding, P. and Timme, S. "Tropical Endgame", In preparation (2020)
"""
Base.@kwdef mutable struct PathTrackerOptions
    endgame_start::Float64 = 0.1
    # eg parameter
    min_cond_eg::Float64 = 1e6
    # valuation etc
    zero_is_at_infinity::Bool = false
    at_infinity_check::Bool = true
    val_finite_tol::Float64 = 1e-2
    val_at_infinity_tol::Float64 = 1e-3
    # singular solutions parameters
    max_winding_number::Int = 20
end

Base.show(io::IO, opts::PathTrackerOptions) = print_fieldnames(io, opts)
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::PathTrackerOptions) = opts

"""
    PathTrackerCode

The possible states a `PathTracker` can be in:

* `PathTrackerCode.tracking`
* `PathTrackerCode.success`
* `PathTrackerCode.at_infinity`
* `PathTrackerCode.at_zero`
* `PathTrackerCode.excess_solution`
* `PathTrackerCode.post_check_failed`
* `PathTrackerCode.polyhedral_failed`
* `PathTrackerCode.terminated_accuracy_limit`
* `PathTrackerCode.terminated_ill_conditioned`
* `PathTrackerCode.terminated_invalid_startvalue`
* `PathTrackerCode.terminated_max_winding_number`
* `PathTrackerCode.terminated_max_steps`
* `PathTrackerCode.terminated_step_size_too_small`
"""
module PathTrackerCode
import ..TrackerCode

@enum codes begin
    tracking
    success
    at_infinity
    at_zero
    terminated_accuracy_limit
    terminated_invalid_startvalue
    terminated_ill_conditioned
    terminated_max_steps
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
    elseif code == TrackerCode.terminated_max_steps
        return terminated_max_steps
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
is_at_infinity(code::PathTrackerCode.codes) = code == PathTrackerCode.at_infinity

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

###
### State
###
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

"""
    PathTracker(tracker::Tracker; options = PathTrackerOptions())
    PathTracker(H::AbstractHomotopy; options = PathTrackerOptions())

A `PathTracker` combines a [`Tracker`](@ref) with an endgame. That is,
while a [`Tracker`](@ref) assumes that the solution path is non-singular and convergent, the endgame
allows to handle singular endpoints as well as diverging paths.
To compute singular solutions the *Cauchy endgame* used, for divering paths a strategy
based on the valuation of local Puiseux series expansion of the path is used.
See [^BT20] for a detailed description.
By convention, a `PathTracker` always tracks from ``t=1`` to ``t = 0``.
See [`PathTrackerOptions`](@ref) for the possible options.

[^BT20]: Breiding, P. and Timme, S. "Tropical Endgame", In preparation (2020)
"""
struct PathTracker{
    H<:AbstractHomotopy,
    N, # AutomaticDifferentiation
    # V and V̄ need to have the same container type
    V<:AbstractVector{ComplexF64},
    V̄<:AbstractVector{ComplexDF64},
    M<:AbstractMatrix{ComplexF64},
} <: AbstractTracker
    tracker::Tracker{H,N,V,V̄,M}
    state::PathTrackerState{V}
    options::PathTrackerOptions
end

PathTracker(H::Homotopy; kwargs...) = PathTracker(ModelKitHomotopy(H); kwargs...)
function PathTracker(
    H::AbstractHomotopy;
    automatic_differentiation = 3,
    tracker_options = TrackerOptions(),
    kwargs...,
)
    tracker = Tracker(
        H;
        automatic_differentiation = automatic_differentiation,
        options = tracker_options,
    )
    PathTracker(tracker; kwargs...)
end
function PathTracker(tracker::Tracker; options = PathTrackerOptions())
    state = PathTrackerState(size(tracker.homotopy, 1), tracker.state.x)
    PathTracker(tracker, state, options)
end

Base.broadcastable(T::PathTracker) = Ref(T)

function init!(path_tracker::PathTracker, x, t₁::Real; ω::Float64 = NaN, μ::Float64 = NaN)
    @unpack tracker, state, options = path_tracker

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

    path_tracker
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
    cauchy!(state::PathTrackerState, tracker::Tracker, options::PathTrackerOptions)

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

            if !is_success(res)
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

function step!(path_tracker::PathTracker, debug::Bool = false)
    @unpack tracker, state, options = path_tracker

    proposed_t′ = real(tracker.state.t′)
    state.last_point .= tracker.state.x
    state.last_t = tracker.state.t
    step!(tracker)

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
                zero_is_finite = !options.zero_is_at_infinity,
            )

            if options.at_infinity_check && verdict.at_infinity
                return (state.code = PathTrackerCode.at_infinity)
            elseif options.zero_is_at_infinity && verdict.at_zero
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
        row_scaling!(state.row_scaling, jacobian(tracker.state.jacobian), state.col_scaling)
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
        zero_is_finite = !options.zero_is_at_infinity,
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

    elseif options.zero_is_at_infinity && at_zero && state.cond > options.min_cond_eg
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
    path_tracker::PathTracker,
    x,
    t₁::Real;
    ω::Float64 = NaN,
    μ::Float64 = NaN,
    debug::Bool = false,
)
    init!(path_tracker, x, t₁; ω = ω, μ = μ)

    while is_tracking(path_tracker.state.code)
        step!(path_tracker, debug)
    end

    path_tracker.state.code
end

function solution(path_tracker::PathTracker)
    get_solution(path_tracker.tracker.homotopy, path_tracker.state.solution, 0.0)
end

"""
    PathResult{V<:AbstractVector}

A `PathResult` is the result of tracking of a path with [`track`](@ref) using a
[`PathTracker`](@ref).

## Fields

General solution information:
* `return_code`: See the list of return codes below.
* `solution::V`: The solution vector.
* `t::Float64`: The value of `t` at which `solution` was computed. Note that if
  `return_code` is `:at_infinity`, then `t` is the value when this was decided.
* `accuracy::Float64`: An estimate the (relative) accuracy of the computed solution.
* `residual::Float64`: The infinity norm of `H(solution,t)`.
* `condition_jacobian::Float64`: This is the condition number of the Jacobian at the
  solution. A high condition number indicates a singular solution or a
  solution on a positive dimensional component.
* `winding_number:Union{Nothing, Int}`: The computed winding number. This is a lower bound
  on the multiplicity of the solution. It is ``nothing`` if the Cauchy endgame was not used.
* `extended_precision::Bool`: Indicate whether extended precision is necessary to achieve
  the accuracy of the `solution`.
* `path_number::Union{Nothing,Int}`: The number of the path (optional).
* `start_solution::Union{Nothing,V}`: The start solution of the path (optional).

Performance information:
* `accepted_steps::Int`: The number of accepted steps during the path tracking.
* `rejected_steps::Int`: The number of rejected steps during the path tracking.
* `extended_precision_used::Bool`: Indicates whether extended precision was
  necessary to track the path.

Additional path and solution informations
* `valuation::Vector{Float64}`: An approximation of the valuation of the
  Puiseux series expansion of ``x(t)``.
* `last_path_point::Tuple{V,Float64}`: The last pair ``(x,t)`` before the solution was
  computed. If the solution was computed with the Cauchy endgame, then the pair ``(x,t)``
  can be used to rerun the endgame.

## Return codes

Possible return codes are:
* `:success`: The `PathTracker` obtained a solution.
* `:at_infinity`: The `PathTracker` stopped the tracking of the path since it determined
  that that path is diverging towards infinity.
* `:at_zero`: The `PathTracker` stopped the tracking of the path since it determined
  that that path has a solution where at least one coordinate is 0. This only happens if
  the option `zero_is_at_infinity` is `true`.
* `:excess_solution`: For the solution of the system, the system had to be
  modified which introduced artificial solutions and this solution is one of them.
* various return codes indicating termination of the tracking
"""
Base.@kwdef mutable struct PathResult{V<:AbstractVector}
    return_code::Symbol
    solution::V
    t::Float64
    accuracy::Float64
    residual::Float64
    multiplicity::Union{Nothing,Int}
    condition_jacobian::Float64
    winding_number::Union{Nothing,Int}
    extended_precision::Bool
    path_number::Union{Nothing,Int}
    start_solution
    last_path_point::Tuple{V,Float64}
    valuation::Union{Nothing,Vector{Float64}}
    ω::Float64
    μ::Float64
    # performance stats
    accepted_steps::Int
    rejected_steps::Int
    extended_precision_used::Bool
end

function PathResult(
    path_tracker::PathTracker,
    start_solution = nothing,
    path_number = nothing,
)
    @unpack tracker, state, options = path_tracker
    H = tracker.homotopy


    if is_success(state.code)
        t = 0.0
        solution = get_solution(H, state.solution, 0.0)
        evaluate!(tracker.corrector.r, H, state.solution, complex(0.0))
        residual = LA.norm(tracker.corrector.r, InfNorm())
    else
        t = real(tracker.state.t)
        solution = get_solution(H, tracker.state.x, t)
        evaluate!(tracker.corrector.r, H, tracker.state.x, complex(t))
        residual = LA.norm(tracker.corrector.r, InfNorm())
    end

    PathResult(
        return_code = Symbol(state.code),
        solution = solution,
        t = t,
        accuracy = state.accuracy,
        residual = residual,
        condition_jacobian = state.cond,
        winding_number = state.winding_number,
        multiplicity = nothing,
        last_path_point = (get_solution(H, state.last_point, state.last_t), state.last_t),
        valuation = t > options.endgame_start ? nothing : copy(state.val.val_x),
        start_solution = start_solution,
        path_number = path_number,
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

"""
    result_type(tracker::PathTracker)

Returns the type of result `track` will return.
"""
result_type(tracker::PathTracker) = PathResult{typeof(solution(tracker))}

"""
    track(path_tracker::PathTracker, x::AbstractVector, t::Real = 1.0;
          path_number = nothing, debug = false)

Track the given start solution `x` from `t` towards `0` using the given `path_tracker`.
Returns a [`PathResult`](@ref).

    track(path_tracker::PathTracker, r::PathResult, t::Real = 1.0;
          path_number = nothing, debug = false)

Track `solution(r)` from `t` towards `0` using the given `path_tracker`.
"""
function track(
    path_tracker::PathTracker,
    x::AbstractVector,
    t₁::Real = 1.0;
    path_number::Union{Nothing,Int} = nothing,
    kwargs...,
)
    track!(path_tracker, x, t₁; kwargs...)
    PathResult(path_tracker, x, path_number)
end
function track(path_tracker::PathTracker, r::PathResult, t₁::Real; kwargs...)
    track(path_tracker, solution(r), t₁; ω = r.ω, μ = r.μ, kwargs...)
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
    residual(r::PathResult)

Get the residual of the solution.
"""
residual(r::PathResult) = r.residual

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
    cond(r::PathResult)

Return the condition number of the Jacobian of the result. A high condition number indicates
that the solution is singular or on a positive dimensional component.
"""
LA.cond(r::PathResult) = r.condition_jacobian

"""
    path_number(r::PathResult)

Get the number of the path. Returns `nothing` if it wasn't provided.
"""
path_number(r::PathResult) = r.path_number

"""
    start_solution(r::PathResult)

Get the start solution of the path.
"""
start_solution(r::PathResult) = r.start_solution

"""
    multiplicity(r::PathResult)

Get the multiplicity of the solution of the path. Returns `nothing` if it wasn't computed.
"""
multiplicity(r::PathResult) = r.multiplicity

"""
    last_path_point(r::PathResult)

Returns a tuple `(x,t)` containing the last zero of `H(x, t)` before the Cauchy endgame was used.
Returns `nothing` if the endgame strategy was not invoked.
"""
last_path_point(r::PathResult) = r.last_path_point

"""
    valuation(r::PathResult)

Get the computed valuation of the path.
"""
valuation(r::PathResult) = r.valuation

"""
    is_success(r::PathResult)

Checks whether the path is successfull.
"""
is_success(r::PathResult) = r.return_code == :success

"""
    is_at_infinity(r::PathResult)

Checks whether the path goes to infinity.
"""
is_at_infinity(r::PathResult) = r.return_code == :at_infinity || r.return_code == :at_zero

"""
    is_excess_solution(r::PathResult)

Checks whether the path is successfull.
"""
is_excess_solution(r::PathResult) = r.return_code == :excess_solution

"""
    is_failed(r::PathResult)

Checks whether the path failed.
"""
is_failed(r::PathResult) = !(is_at_infinity(r) || is_success(r) || is_excess_solution(r))

"""
    is_finite(r::PathResult)

Checks whether the path result is finite.
"""
is_finite(r::PathResult) = is_success(r)
Base.isfinite(r::PathResult) = is_finite(r)

"""
    is_singular(r::PathResult; tol::Float64 = 1e10)

Checks whether the path result `r` is singular. This is true if
the winding number is larger than  1 or if the condition number of the Jacobian
is larger than `tol`.
"""
is_singular(r::PathResult; tol::Float64 = 1e10) = is_singular(r, tol)
function is_singular(r::PathResult, tol::Real)
    (unpack(r.condition_jacobian, 1.0) > tol || unpack(winding_number(r), 1) > 1) &&
    is_success(r)
end

"""
    is_nonsingular(r::PathResult; tol::Float64 = 1e10)

Checks whether the path result is non-singular. This is true if
it is not singular.
"""
is_nonsingular(r::PathResult; kwargs...) = !is_singular(r; kwargs...) && is_success(r)
is_nonsingular(r::PathResult, tol::Real) = !is_singular(r, tol) && is_success(r)


"""
    is_real(r::PathResult; tol::Float64 = 1e-6)

We consider a result as `real` if the infinity-norm of the imaginary part of the solution
is at most `tol`.
"""
is_real(r::PathResult; tol::Float64 = 1e-6) = is_real(r, tol)
is_real(r::PathResult, tol::Real) = maximum(abs ∘ imag, r.solution) < tol
# provide fallback since this in in Base
Base.isreal(r::PathResult, tol) = is_real(r, tol)
Base.isreal(r::PathResult; kwargs...) = is_real(r; kwargs...)
