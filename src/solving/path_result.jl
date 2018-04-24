import ..Homotopies
using ..Utilities

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
    t::Float64

    residual::Float64
    condition_number::Float64
    multiplicity::Int

    start_solution::Vector{T}

    iterations::Int
    npredictions::Int
end


function PathResult(::Problems.NullHomogenization,
    H::Homotopies.HomotopyWithCache,
    x₁, t₀, x, t, returncode::Symbol, iters::Int, multiplicity::Int, npredictions::Int, v, J)

    Homotopies.evaluate_and_jacobian!(v, J, H, x₀, t₀)
    res = infinity_norm(v)
    condition = cond(J)

    PathResult(returncode, x, t, res, condition, multiplicity, x₁, iters, npredictions)
end

function PathResult(::Problems.DefaultHomogenization,
    H::Homotopies.HomotopyWithCache,
    x₁::Vector, t₀, x, t, returncode::Symbol, iters::Int, multiplicity, npredictions::Int, v, J)

    # We want to evaluate on the affine patch we are interested in
    hom_part = x[1]
    scale!(x, inv(x[1]))

    Homotopies.evaluate_and_jacobian!(v, J, H, x, t₀)
    res = infinity_norm(v)
    if res > 0.1
        returncode = :at_infinity
    end

    if returncode == :at_infinity
        solution = x
        condition = 0.0
    else
        solution = x[2:end]
        condition = cond(@view J[:,2:end])
    end

    PathResult(returncode, solution, t, res, condition, multiplicity, x₁, iters, npredictions)
end

function Base.show(io::IO, r::PathResult)
    println(io, typeof(r), ":")
    println(io, "* returncode: $(r.returncode)")
    println(io, "* solution: $(r.solution)")
    println(io, "* t: $(r.t)")
    println(io, "---------------------------------------------")
    println(io, "* iterations: $(r.iterations)")
    println(io, "* npredictions: $(r.npredictions)")
    println(io, "---------------------------------------------")
    println(io, "* residual: $(@sprintf "%.3e" r.residual)")
    println(io, "* log10 of the condition_number: $(@sprintf "%.3e" log10(r.condition_number))")
    println(io, "* multiplicity: $(r.multiplicity)")
end

function pathresults(solver::Solver, results)
    pathresults(solver.prob, results, solver.start_solutions, solver.t₀)
end


function pathresults(prob::Problems.AbstractProblem,
    results::Vector{<:Endgame.EndgamerResult},
    start_solutions, t₀)
    r = results[1]
    H = Homotopies.HomotopyWithCache(prob.homotopy, r.x, t₀)
    v, J = Homotopies.evaluate_and_jacobian(H, r.x, t₀)
    pathresults(prob.homogenization_strategy, H, results, start_solutions, t₀, v, J)
end

 function pathresults(strategy::Problems.AbstractHomogenizationStrategy, H,
     results::Vector{<:Endgame.EndgamerResult},
     start_solutions, t₀, v, J)
     map(results, start_solutions) do r, x₁
         PathResult(strategy, H, x₁, t₀, r.x.data, real(r.t), r.returncode, r.iters, r.windingnumber, r.npredictions, v, J)
     end
 end

function pathresults(prob::Problems.AbstractProblem,
    trackedpath_results::Vector{<:PathTracking.PathTrackerResult},
    start_solutions, t₀)
    r = trackedpath_results[1]
    H = Homotopies.HomotopyWithCache(prob.homotopy, r.x.data, r.t)
    v, J = Homotopies.evaluate_and_jacobian(H, r.x.data, r.t)
    pathresults(prob.homogenization_strategy, H, trackedpath_results, start_solutions, t₀, v, J)
end

function pathresults(strategy::Problems.AbstractHomogenizationStrategy, H,
     trackedpath_results::Vector{<:PathTracking.PathTrackerResult},
     start_solutions, t₀, v, J)
     map(trackedpath_results, start_solutions) do r, x₁
         PathResult(strategy, H, x₁, t₀, r.x.data, real(r.t), r.returncode, r.iters, 1, 0, v, J)
     end
end
