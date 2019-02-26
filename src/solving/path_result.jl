export PathResult, solution,
    residual, startsolution, #issuccess,
    isfailed, isaffine, isprojective,
    isatinfinity, issingular, isnonsingular

struct PathResultCache{Hom, T}
    H::Hom
    v::Vector{T}
    J::Matrix{T}
end

function PathResultCache(prob::AbstractProblem, x)
    H = HomotopyWithCache(prob.homotopy, x, rand())
    v, J = evaluate_and_jacobian(H, x, rand())
    PathResultCache(H, v, J)
end


"""
    PathResult(startvalue, pathtracker_result, endgamer_result, solver)

A `PathResult` is the result of the tracking of a path (inclusive endgame).
Its fields are
* `returncode`: One of `:success`, `:at_infinity` or any error code from the `EndgamerResult`
* `solution::Vector{T}`: The solution vector. If the algorithm computed in projective space and the solution is at infinity then the projective solution is given. Otherwise an affine solution is given if the startvalue was affine and a projective solution is given if the startvalue was projective.
* `residual::Float64`: The value of the infinity norm of `H(solution, 0)`.
* `newton_residual`: The value of the 2-norm of ``J_H(\\text{solution})^{-1}H(\\text{solution}, 0)``
* `condition_number`: This is the condition number of the row-equilibrated Jacobian at the solution. A high condition number indicates a singularity.
* `windingnumber`: The estimated winding number
* `angle_to_infinity`: The angle to infinity is the angle of the solution to the hyperplane where the homogenizing coordinate is ``0``.
* `real_solution`: Indicates whether the solution is real given the defined tolerance `at_infinity_tol` (from the solver options).
* `startvalue`: The startvalue of the path
* `iterations`: The number of iterations the pathtracker needed.
* `endgame_iterations`: The number of steps in the geometric series the endgamer did.
* `npredictions`: The number of predictions the endgamer did.
* `predictions`: The predictions of the endgamer.
"""
struct PathResult{T1, T2, V<:AbstractVector}
    returncode::Symbol
    returncode_detail::Symbol

    solution::Vector{T1}
    solution_type::Symbol
    t::T2

    residual::Float64
    condition_number::Float64
    windingnumber::Int

    pathnumber::Int
    start_solution::V
    endgame_start_solution::Vector{T1}

    iterations::Int
    npredictions::Int
end

function PathResult(prob::ProjectiveProblem, k, x₁, x_e, t₀, r, cache::PathResultCache)
    PathResult(homvars(prob), prob.vargroups, k, x₁, x_e, t₀, r, cache)
end
function PathResult(homvars::Nothing, vargroups, k, x₁, x_e, t₀, r, cache::PathResultCache)
    returncode, returncode_detail = makereturncode(r.returncode)
    x = align_axis!(copy(r.x.data))
    evaluate_and_jacobian!(cache.v, cache.J, cache.H, x, t₀)
    res = infinity_norm(cache.v)


    if returncode != :success
        condition = Inf
    else
        # Before we compute the condition number we row-equilibrate the Jacobian matrix
        for i=1:size(cache.J, 1)
            rᵢ = abs2(cache.J[i, 1])
            for j=2:size(cache.J, 2)
                rᵢ += abs2(cache.J[i, j])
            end
            rᵢ = inv(√rᵢ)
            for j=1:size(cache.J, 2)
                cache.J[i, j] *= rᵢ
            end
        end
        σ = LinearAlgebra.svdvals(cache.J)
        if size(cache.J, 1) > ngroups(vargroups) + 1
            condition = σ[1]/σ[end-ngroups(vargroups)]
        else
            condition = inv(σ[1])
        end
    end

    windingnumber, npredictions = windingnumber_npredictions(r)

    PathResult(returncode, returncode_detail, x, :projective, real(r.t), res, condition,
        windingnumber, k, x₁, x_e.data, iters(r), npredictions)
end
function PathResult(homvars::NTuple{M,Int}, vargroups, k, x₁, x_e, t₀, r, cache::PathResultCache) where {M}
    returncode = Symbol(r.returncode)
    windingnumber, npredictions = windingnumber_npredictions(r)

    evaluate_and_jacobian!(cache.v, cache.J, cache.H, r.x, t₀)
    res = infinity_norm(cache.v)

    returncode, returncode_detail = makereturncode(returncode)

    intermediate_sol = ProjectiveVectors.affine_chart(x_e)

    if returncode != :success
        solution = copy(r.x.data)
        condition = Inf
    else
        solution = ProjectiveVectors.affine_chart(r.x)

        # Before we compute the condition number we row-equilibrate the Jacobian matrix
        for i=1:size(cache.J, 1)
            rᵢ = abs2(cache.J[i, 1])
            for j=2:size(cache.J, 2)-1
                rᵢ += abs2(cache.J[i, j])
            end
            rᵢ = inv(√rᵢ)
            for j=1:size(cache.J, 2)-1
                cache.J[i, j] *= rᵢ
            end
        end
        condition = LinearAlgebra.cond(cache.J)
    end

    if x₁ isa ProjectiveVectors.PVector
        start = ProjectiveVectors.affine_chart(x₁)
    else
        start = x₁
    end

    PathResult(returncode, returncode_detail, solution, :affine, real(r.t), res,
        condition, windingnumber, k, start, intermediate_sol, iters(r), npredictions)
end

iters(R::EndgameResult) = R.iters
iters(R::PathTrackerResult) = R.accepted_steps + R.rejected_steps


function makereturncode(retcode)
    if retcode != :at_infinity && retcode != :success
        :path_failed, retcode
    else
        retcode, :none
    end
end
makereturncode(retcode::PathTrackerStatus.states) = :path_failed, :none


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

windingnumber_npredictions(r::EndgameResult) = (r.windingnumber, r.npredictions)
windingnumber_npredictions(r::PathTrackerResult) = (0, 0)

function Base.show(io::IO, r::PathResult)
    iscompact = get(io, :compact, false)
    if iscompact || haskey(io, :typeinfo)
        println(io, "• returncode: $(r.returncode)")
        println(io, " • solution: ", r.solution)
        println(io, " • residual: $(Printf.@sprintf "%.3e" r.residual)")
        println(io, " • pathnumber: ", r.pathnumber)
    else
        println(io, "PathResult")
        println(io, "==========")
        println(io, " • returncode: $(r.returncode)")
        if r.returncode_detail != :none
            println(io, " • returncode_detail: $(r.returncode_detail)")
        end
        println(io, " • solution: ", r.solution)
        println(io, " • residual: $(Printf.@sprintf "%.3e" r.residual)")
        println(io, " • condition_number: $(Printf.@sprintf "%.3e" r.condition_number)")
        println(io, " • windingnumber: $(r.windingnumber)")
        println(io, "")
        println(io, " • pathnumber: ", r.pathnumber)
        println(io, " • start_solution: ", r.start_solution)
        println(io, "")
        println(io, " • t: ", r.t)
        println(io, " • iterations: ", r.iterations)
        println(io, " • npredictions: $(r.npredictions)")
    end
end

TreeViews.hastreeview(::PathResult) = true
TreeViews.treelabel(io::IO, x::PathResult, ::MIME"application/prs.juno.inline") =
    print(io, "<span class=\"syntax--support syntax--type syntax--julia\">PathResult</span>")

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
    startsolution(pathresult)

Get the start solution of the solution ``x`` of the path.
"""
startsolution(r::PathResult) = r.start_solution

"""
    issuccess(pathresult)

Checks whether the path is successfull.
"""
LinearAlgebra.issuccess(r::PathResult) = r.returncode == :success

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
Base.isfinite(r::PathResult) = r.returncode == :success # we don't check isaffine to make other code easier

"""
    issingular(pathresult;tol=1e14)
    issingular(pathresult, tol)

Checks whether the path result is singular. This is true if
the winding number > 1 or if the condition number of the Jacobian
is larger than `tol`.
"""
issingular(r::PathResult;tol=1e14) = issingular(r, tol)
function issingular(r::PathResult, tol::Real)
    (r.windingnumber > 1 || r.condition_number > tol) && LinearAlgebra.issuccess(r)
end

"""
    isnonsingular(pathresult;tol=1e14)

Checks whether the path result is non-singular. This is true if
it is not singular.
"""
isnonsingular(r::PathResult;tol=1e14) = isnonsingular(r, tol)
isnonsingular(r::PathResult, tol::Real) = !issingular(r, tol) && LinearAlgebra.issuccess(r)


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
