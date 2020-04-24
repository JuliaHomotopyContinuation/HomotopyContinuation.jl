export PathResult,
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
    PathResult

A `PathResult` is the result of tracking of a path with [`track`](@ref) using an [`AbstractPathTracker`](@ref) (
e.g. [`PathTracker`](@ref))

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
Base.@kwdef mutable struct PathResult
    return_code::Symbol
    solution::Vector{ComplexF64}
    t::Float64
    accuracy::Float64
    residual::Float64
    multiplicity::Union{Nothing,Int}
    condition_jacobian::Float64
    winding_number::Union{Nothing,Int}
    extended_precision::Bool
    path_number::Union{Nothing,Int}
    start_solution
    last_path_point::Tuple{Vector{ComplexF64},Float64}
    valuation::Union{Nothing,Vector{Float64}}
    ω::Float64
    μ::Float64
    # performance stats
    accepted_steps::Int
    rejected_steps::Int
    extended_precision_used::Bool
end

Base.show(io::IO, r::PathResult) = print_fieldnames(io, r)
Base.show(io::IO, ::MIME"application/prs.juno.inline", r::PathResult) = r

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
