###
### TrackerCode
###

@doc """
    TrackerCode

The possible states a `CoreTracker` can have are of type `TrackerCode.codes` and can be

* `TrackerCode.success`: Indicates a successfull tracking.
* `TrackerCode.tracking`: The tracking is still in progress.
* `TrackerCode.terminated_accuracy_limit`: Tracking terminaed since the accuracy was insufficient.
* `TrackerCode.terminated_invalid_startvalue`: Tracking terminated since the provided start value was invalid.
* `TrackerCode.terminated_ill_conditioned`: Tracking terminated since the path was too ill-conditioned.
* `TrackerCode.terminated_max_steps`: Tracking terminated since maximal number of steps is reached.
* `TrackerCode.terminated_step_size_too_small`: Trackint terminated since the step size was too small.
* `TrackerCode.terminated_unknown`: An unintended error occured. Please consider reporting an issue.
"""
module TrackerCode

@enum codes begin
    tracking
    success
    terminated_max_steps
    terminated_accuracy_limit
    terminated_ill_conditioned
    terminated_invalid_startvalue
    terminated_invalid_startvalue_singular_jacobian
    terminated_step_size_too_small
    terminated_unknown
end

end

"""
    is_success(code::TrackerCode.codes)

Returns `true` if `code` indicates a success in the path tracking.
"""
is_success(S::TrackerCode.codes) = S == TrackerCode.success

"""
    is_terminated(code::TrackerCode.codes)

Returns `true` if `code` indicates that the path tracking got terminated.
"""
is_terminated(S::TrackerCode.codes) = S ≠ TrackerCode.tracking && S ≠ TrackerCode.success

"""
    is_invalid_startvalue(code::TrackerCode.codes)

Returns `true` if `code` indicates that the path tracking got terminated since the start
value was not a regular zero.
"""
function is_invalid_startvalue(S::TrackerCode.codes)
    S == TrackerCode.terminated_invalid_startvalue ||
        S == TrackerCode.terminated_invalid_startvalue_singular_jacobian
end

"""
    is_tracking(code::TrackerCode.codes)

Returns `true` if `code` indicates that the path tracking is not yet finished.
"""
is_tracking(S::TrackerCode.codes) = S == TrackerCode.tracking

###
### RESULT
###

"""
    TrackerResult

Containing the result of tracking a path with a [`Tracker`](@ref).

## Fields

* `return_code::Symbol`: A code indicating whether the tracking was successfull (`:success`).
  See [`TrackerCode`](@ref) for all possible values.
* `solution::V`: The solution when the tracking stopped.
* `t::ComplexF64`: The value of `t` when the tracking stopped.
* `accuracy::Float64`: Estimate of the relative accuracy of the `solution`.
* `accepted_steps::Int`: Number of steps that got accepted.
* `rejected_steps::Int`: Number of steps that got rejected.
* `extended_precision::Bool`: Indicate whether extended precision is necessary to achieve
  the accuracy of the `solution`.
* `extended_precision_used::Bool`: This is `true` if during the tracking at any point
  extended precision was used.
"""
struct TrackerResult
    return_code::Symbol
    solution::Vector{ComplexF64}
    t::ComplexF64
    accuracy::Float64
    accepted_steps::Int
    rejected_steps::Int
    extended_precision::Bool
    extended_precision_used::Bool
    ω::Float64
    μ::Float64
    τ::Float64
end

Base.show(io::IO, result::TrackerResult) = print_fieldnames(io, result)
Base.show(io::IO, ::MIME"application/prs.juno.inline", result::TrackerResult) = result

"""
    is_success(result::TrackerResult)

Returns `true` if the path tracking was successfull.
"""
is_success(R::TrackerResult) = R.return_code == :success

"""
    is_invalid_startvalue(result::TrackerResult)

Returns `true` if the path tracking failed since the start value was invalid.
You can inspect `result.return_code` to get the exact return code. Possible values
if `is_invalid_startvalue(result) == true` are
* `:terminated_invalid_startvalue_singular_jacobian` indicates that the Jacobian of the homotopy at
  the provided start value is singular, i.e., if it has not full-column rank.
* `:terminated_invalid_startvalue` indicates that the the provided start value is not sufficiently
  close to a solution of the homotopy.
"""
function is_invalid_startvalue(R::TrackerResult)
    R.return_code == :terminated_invalid_startvalue ||
        R.return_code == :terminated_invalid_startvalue_singular_jacobian
end

"""
    solution(result::TrackerResult)

Returns the solutions obtained by the `Tracker`.
"""
solution(result::TrackerResult) = result.solution

"""
    steps(result::TrackerResult)

Returns the number of steps done.
"""
steps(result::TrackerResult) = accepted_steps(result) + rejected_steps(result)

"""
    accepted_steps(result::TrackerResult)

Returns the number of accepted steps.
"""
accepted_steps(result::TrackerResult) = result.accepted_steps

"""
    rejected_steps(result::TrackerResult)

Returns the number of rejected_steps steps.
"""
rejected_steps(result::TrackerResult) = result.rejected_steps
