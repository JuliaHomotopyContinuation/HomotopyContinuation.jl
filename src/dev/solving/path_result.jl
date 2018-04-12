import ..NewHomotopies

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
struct PathResult{T}
    returncode::Symbol
    solution::Vector{T}
    # isolated::Bool
    # singular::Bool

    residual::Float64
    newton_residual::Float64
    condition_number::Float64
    # windingnumber::Int
    angle_to_infinity::Float64
    # real::Bool

    start_solution::Vector{T}
    #
    iterations::Int
    # endgame_iterations::Int
    # npredictions::Int
    # predictions::Vector{Vector{T}}
end


function PathResult(::Problems.NullHomogenization,
    H::NewHomotopies.HomotopyWithCache,
    pathtracker_result::PathTracking.PathTrackerResult{<:ProjectiveVector}, x₁, t₀, v, J)
    realtol = 1e-6
    angle_to_infinity_tol = 1e-8
    # we need to make the solution affine + check for infinity
    x₀ = pathtracker_result.x.data
    returncode = pathtracker_result.returncode
    iters = pathtracker_result.iters

    angle_to_infinity = NaN

    NewHomotopies.evaluate_and_jacobian!(v, J, H, x₀, t₀)
    res = norm(v, Inf)
    newton_res = norm(J \ v, Inf)
    condition = cond(J)

    PathResult(returncode, solution, res, newton_res, condition, angle_to_infinity, x₁, iters)
end

function PathResult(::Problems.DefaultHomogenization,
    H::NewHomotopies.HomotopyWithCache,
    pathtracker_result::PathTracking.PathTrackerResult{<:ProjectiveVector}, x₁, t₀, v, J)
    realtol = 1e-6
    angle_to_infinity_tol = 1e-8
    # we need to make the solution affine + check for infinity
    x₀ = pathtracker_result.x.data
    returncode = pathtracker_result.returncode
    iters = pathtracker_result.iters

    hom_part = x₀[1]
    affine_part = x₀[2:end]

    angle_to_infinity = atan(abs(hom_part)/norm(affine_part))

    if angle_to_infinity < angle_to_infinity_tol
        returncode = :at_infinity
    end

    if returncode == :at_infinity
        res = Inf
        newton_res = Inf
        solution = x₀
        # isreal = false
        condition_number = 0.0
    else
        # We want to evaluate on the affine patch we are interested in
        scale!(x₀, inv(hom_part))
        NewHomotopies.evaluate_and_jacobian!(v, J, H, x₀, t₀)
        res = norm(v, Inf)
        newton_res = norm(J \ v, Inf)
        scale!(affine_part, inv(hom_part))
        solution = affine_part
        # isreal = maximum(z -> abs(imag(z)), z) < realtol

        J_affine = @view J[:,2:end]
        condition = cond(J_affine)
    end

    PathResult(returncode, solution, res, newton_res, condition, angle_to_infinity, x₁, iters)
end

function pathresults(solver::Solver, results)
    pathresults(solver.prob, results, solver.start_solutions, solver.t₀)
end

function pathresults(prob::Problems.AbstractProblem, trackedpath_results, start_solutions, t₀)
    r = trackedpath_results[1]
    H = NewHomotopies.HomotopyWithCache(prob.homotopy, r.x.data, r.t)
    v, J = NewHomotopies.evaluate_and_jacobian(H, r.x.data, r.t)
    pathresults(prob.homogenization_strategy, H, trackedpath_results, start_solutions, t₀, v, J)
end

 function pathresults(strategy::Problems.AbstractHomogenizationStrategy,
     H, trackedpath_results, start_solutions, t₀, v, J)
     map(trackedpath_results, start_solutions) do r, x₁
         PathResult(strategy, H, r, x₁, t₀, v, J)
     end
 end
