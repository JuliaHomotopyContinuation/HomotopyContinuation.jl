import MultivariatePolynomials
const MP = MultivariatePolynomials

import .Problems
import .Systems
import .NewHomotopies
import .Solving

export solve

# External
function solve(F::Vector{<:MP.AbstractPolynomial}; kwargs...)
    TDP = Problems.TotalDegreeProblem(F)
    start_solutions = Utilities.totaldegree_solutions(F)
    solve(TDP, start_solutions; kwargs...)
end

function solve(G::Vector{<:MP.AbstractPolynomial}, F::Vector{<:MP.AbstractPolynomial}, start_solutions; kwargs...)
    STP = Problems.StartTargetProblem(G, F)
    solve(STP, start_solutions; kwargs...)
end

# Internal
function solve(prob::Problems.AbstractDynamicProblem, start_solutions;
    system=Systems.SPSystem,
    homotopy=NewHomotopies.StraightLineHomotopy,
    kwargs...)
    P = Problems.ProjectiveStartTargetProblem(prob, system=system, homotopy=homotopy)
    solve(P, start_solutions; kwargs...)
end

function solve(prob::Problems.AbstractProblem, start_solutions, t₁=1.0, t₀=0.0; kwargs...)
    solver = Solving.Solver(prob, start_solutions, t₁, t₀)
    Solving.solve(solver)
end
