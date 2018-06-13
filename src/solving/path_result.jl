export PathResult, solution,
    residual, start_solution, issuccess,
    isfailed, isaffine, isprojective,
    isatinfinity, issingular, issmooth


using Compat

import ..Homotopies
import ..ProjectiveVectors
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
* `condition_number`: A high condition number indicates singularty. See [`Homotopies.κ`](@ref) for details.
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
    solution_type::Symbol
    t::T2

    residual::Float64
    condition_number::Float64
    windingnumber::Int

    pathnumber::Int
    start_solution::Vector{T3}
    endgame_start_solution::Vector{T1}

    iterations::Int
    npredictions::Int
end

function PathResult(prob::Problems.AbstractProblem, k, x₁, x_e, t₀, r, cache::PathResultCache, patchswitcher)
    PathResult(prob.homogenization, k, x₁, x_e, t₀, r, cache, patchswitcher)
end
function PathResult(::Problems.NullHomogenization, k, x₁, x_e, t₀, r, cache::PathResultCache, ::Nothing)
    returncode, returncode_detail = makereturncode(r.returncode)
    x = raw(align_axis!(copy(r.x)))
    Homotopies.evaluate_and_jacobian!(cache.v, cache.J, cache.H, x, t₀)
    res = infinity_norm(cache.v)

    if returncode != :success
        condition = 0.0
    else
        condition = cond(cache.J)
    end

    windingnumber, npredictions = windingnumber_npredictions(r)


    PathResult(returncode, returncode_detail, x, :projective, real(r.t), res, condition,
        windingnumber, k, x₁, raw(x_e), r.iters, npredictions)
end
function PathResult(::Problems.Homogenization, k, x₁, x_e, t₀, r, cache::PathResultCache, patchswitcher::PatchSwitching.PatchSwitcher)
    returncode = r.returncode
    windingnumber, npredictions = windingnumber_npredictions(r)

    switch_to_affine!(r.x, returncode, windingnumber, patchswitcher)

    Homotopies.evaluate_and_jacobian!(cache.v, cache.J, cache.H, raw(r.x), t₀)
    res = infinity_norm(cache.v)

    returncode, returncode_detail = makereturncode(returncode)

    intermediate_sol = ProjectiveVectors.affine(x_e)

    if returncode != :success
        solution = copy(raw(r.x))
        condition = 0.0
    else
        solution = ProjectiveVectors.affine(r.x)
        condition = cond(@view cache.J[:,2:end])
    end


    PathResult(returncode, returncode_detail, solution, :affine, real(r.t), res,
        condition, windingnumber, k, x₁, intermediate_sol, r.iters, npredictions)
end

function makereturncode(retcode)
    if retcode != :at_infinity && retcode != :success
        :path_failed, retcode
    else
        retcode, :none
    end
end

function switch_to_affine!(x::ProjectiveVectors.PVector{<:Complex, Int}, returncode, windingnumber, patchswitcher)
    if returncode == :success && windingnumber == 1 && abs2(x[x.homvar]) < 1.0
        PatchSwitching.switch!(x, patchswitcher)
    elseif returncode != :at_infinity
        ProjectiveVectors.affine!(x)
    end
    x
end

"""
    align_axis!(x)

Multiplies `x` with a complex number of norm 1 such that the largest
entry is real.
"""
function align_axis!(x)
    maxval = abs2(x[1])
    maxind = 1
    for i=2:length(x)
        val = abs2(x[i])
        if val > maxval
            maxval = val
            maxind = i
        end
    end

    v = conj(x[maxind]) / sqrt(maxval)
    for i=1:length(x)
        x[i] *= v
    end
    x
end

windingnumber_npredictions(r::Endgaming.Result) = (r.windingnumber, r.npredictions)
windingnumber_npredictions(r::PathTracking.PathTrackerResult) = (0, 0)

function Base.show(io::IO, r::PathResult)
    iscompact = get(io, :compact, false)

    if iscompact
        println(io, "* returncode: $(r.returncode)")
        println(io, " * solution: ", r.solution)
        println(io, " * residual: $(@sprintf "%.3e" r.residual)")
    else
        println(io, " ---------------------------------------------")
        println(io, "* returncode: $(r.returncode)")
        if r.returncode_detail != :none
            println(io, " * returncode_detail: $(r.returncode_detail)")
        end
        println(io, " * solution: ", r.solution)
        println(io, " * residual: $(@sprintf "%.3e" r.residual)")
        println(io, " * condition_number: $(@sprintf "%.3e" r.condition_number)")
        println(io, " * windingnumber: $(r.windingnumber)")
        println(io, " ----")
        println(io, " * path number: ", r.pathnumber)
        println(io, " * start_solution: ", r.start_solution)
        println(io, " ----")
        println(io, " * t: ", r.t)
        println(io, " * iterations: ", r.iterations)
        println(io, " * npredictions: $(r.npredictions)")


    end

end

function Juno.render(i::Juno.Inline, x::PathResult{T1, T2, T3}) where {T1, T2, T3}
    t = Juno.render(i, Juno.defaultrepr(x, true))
    t[:head] = Juno.render(i, Text("PathResult"))
    t
end

"""
    solution(pathresult)

Get the solution of the path.
"""
solution(r::PathResult) = r.solution

"""
    residual(pathresult)

Get the residual of the solution ``x`` of the path, i.e., ``||H(x, t)||_{\\infty}``.
"""
residual(r::PathResult) = r.residual

"""
    residual(pathresult)

Get the residual of the solution ``x`` of the path, i.e., ``||H(x, t)||_{\\infty}``.
"""
start_solution(r::PathResult) = r.start_solution

"""
    issuccess(pathresult)

Checks whether the path is successfull.
"""
issuccess(r::PathResult) = r.returncode == :success

"""
    isfailed(pathresult)

Checks whether the path failed.
"""
isfailed(r::PathResult) = r.returncode == :path_failed

"""
    isaffine(pathresult)

Checks whether the path result is affine.
"""
isaffine(r::PathResult) = r.solution_type == :affine

"""
    isprojective(pathresult)

Checks whether the path result is projective.
"""
isprojective(r::PathResult) = r.solution_type == :projective

"""
    isatinfinity(pathresult)

Checks whether the path goes to infinity.
"""
isatinfinity(r::PathResult) = (r.returncode == :at_infinity && isaffine(r))


"""
    isfinite(pathresult)

Checks whether the path result is finite.
"""
Base.isfinite(r::PathResult) = (r.returncode == :success && isaffine(r))

"""
    issingular(pathresult; tol=1e10)
    issingular(pathresult, tol)

Checks whether the path result is singular. This is true if
the winding number > 1 or if the condition number of the Jacobian
is larger than `tol`.
"""
issingular(r::PathResult; tol=1e10) = issingular(r, tol)
function issingular(r::PathResult, tol::Real)
    if isprojective(r)
        (r.windingnumber > 1 || r.condition_number > tol) && issuccess(r)
    else
        (r.windingnumber > 1 || r.condition_number > tol) && isfinite(r) && issuccess(r)
    end
end

"""
    issmooth(pathresult; tol=1e10)

Checks whether the path result is smooth. This is true if
it is not singular.
"""
issmooth(r::PathResult; tol=1e10) = issmooth(r, tol)
function issmooth(r::PathResult, tol::Real)
    if isprojective(r)
        r.windingnumber ≤ 1 && r.condition_number ≤ tol
    else
        r.windingnumber ≤ 1 && r.condition_number ≤ tol && isfinite(r)
    end
end


"""
    isreal(pathresult; tol=1e-6)
    isreal(pathresult, tol)

Determine whether infinity norm of the imaginary part of the so
"""
Base.isreal(r::PathResult; tol=1e-6) = isreal(r, tol)
function Base.isreal(r::PathResult, tol::Real)
    sqrt(maximum(_abs2imag, r.solution)) < tol
end
_abs2imag(x) = abs2(imag(x))
