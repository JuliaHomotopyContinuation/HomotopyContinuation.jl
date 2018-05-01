import ..Homotopies
import ..ProjectiveVectors: raw
using ..Utilities


struct PathResultCache{Hom, T}
    H::Hom
    v::Vector{T}
    J::Matrix{T}
end

function PathResultCache(prob::Problems.AbstractProblem, x)
    H = Homotopies.HomotopyWithCache(prob.homotopy, raw(x), rand())
    v, J = Homotopies.evaluate_and_jacobian(H, raw(x), rand())
    PathResultCache(H, v, J)
end


"""
    PathResult(startvalue, pathtracker_result, endgamer_result, solver)

Construct a `PathResult` for a given `startvalue`. `pathtracker_result` is the
[`PathtrackerResult`](@ref) until the endgame radius is reached. `endgamer_result`
is the [`EndgamerResult`](@ref) resulting from the corresponding endgame.

A `PathResult` contains:
* `returncode`: One of `:success`, `:at_infinity` or any error code from the `EndgamerResult`
* `solution::Vector{T}`: The solution vector. If the algorithm computed in projective space
and the solution is at infinity then the projective solution is given. Otherwise
an affine solution is given if the startvalue was affine and a projective solution
is given if the startvalue was projective.
* `residual::Float64`: The value of the infinity norm of `H(solution, 0)`.
* `newton_residual`: The value of the 2-norm of ``J_H(\\text{solution})^{-1}H(\\text{solution}, 0)``
* `log10_condition_number`: A high condition number indicates singularty. See [`Homotopies.κ`](@ref) for details.
    The value is the logarithmic condition number (with base 10).
* `windingnumber`: The estimated winding number
* `angle_to_infinity`: The angle to infinity is the angle of the solution to the hyperplane where the homogenizing coordinate is ``0``.
* `real_solution`: Indicates whether the solution is real given the defined tolerance `at_infinity_tol` (from the solver options).
* `startvalue`: The startvalue of the path
* `iterations`: The number of iterations the pathtracker needed.
* `endgame_iterations`: The number of steps in the geometric series the endgamer did.
* `npredictions`: The number of predictions the endgamer did.
* `predictions`: The predictions of the endgamer.
"""
struct PathResult{T1, T2, T3}
    returncode::Symbol
    returncode_detail::Symbol

    solution::Vector{T1}
    t::T2

    residual::Float64
    condition_number::Float64
    windingnumber::Int

    start_solution::Vector{T3}

    iterations::Int
    npredictions::Int
end


function PathResult(prob::Problems.AbstractProblem, x₁, t₀, r, cache::PathResultCache)
    PathResult(prob.homogenization_strategy, x₁, t₀, r, cache)
end
function PathResult(::Problems.NullHomogenization, x₁, t₀, r, cache::PathResultCache)
    returncode, returncode_detail = makereturncode(r.returncode)
    x = raw(r.x)
    Homotopies.evaluate_and_jacobian!(cache.v, cache.J, cache.H, x, t₀)
    res = infinity_norm(v)

    if returncode != :success
        condition = 0.0
    else
        condition = cond(J)
    end

    windingnumber, npredictions = windingnumber_npredictions(r)

    PathResult(returncode, returncode_detail, x, real(r.t), res, condition, windingnumber, x₁, r.iters, npredictions)
end
function PathResult(::Problems.DefaultHomogenization, x₁, t₀, r, cache::PathResultCache)
    returncode = r.returncode

    ProjectiveVectors.affine!(r.x)
    x = raw(r.x)

    Homotopies.evaluate_and_jacobian!(cache.v, cache.J, cache.H, x, t₀)
    res = infinity_norm(cache.v)
    if res > 0.1 && r.t == t₀ && returncode == :success
        returncode = :at_infinity
    end

    returncode, returncode_detail = makereturncode(returncode)

    if returncode != :success
        solution = x
        condition = 0.0
    else
        solution = x[2:end]
        condition = cond(@view cache.J[:,2:end])
    end

    windingnumber, npredictions = windingnumber_npredictions(r)

    PathResult(returncode, returncode_detail, solution, real(r.t), res, condition, windingnumber, x₁, r.iters, npredictions)
end

function makereturncode(retcode)
    if retcode != :at_infinity && retcode != :success
        :path_failed, retcode
    else
        retcode, :none
    end
end
windingnumber_npredictions(r::Endgame.EndgamerResult) = (r.windingnumber, r.npredictions)
windingnumber_npredictions(r::PathTracking.PathTrackerResult) = (0, 0)

function Base.show(io::IO, ::MIME"text/plain", r::PathResult)
    println(io, typeof(r), ":")
    println(io, "* returncode: $(r.returncode)")
    if r.returncode_detail != :none
        println(io, "* returncode_detail: $(r.returncode_detail)")
    end
    println(io, "* solution: $(r.solution)")
    println(io, "* t: $(r.t)")
    println(io, "---------------------------------------------")
    println(io, "* iterations: $(r.iterations)")
    println(io, "* npredictions: $(r.npredictions)")
    println(io, "---------------------------------------------")
    println(io, "* residual: $(@sprintf "%.3e" r.residual)")
    println(io, "* log10 of the condition_number: $(@sprintf "%.3e" log10(r.condition_number))")
    println(io, "* windingnumber: $(r.windingnumber)")
end
