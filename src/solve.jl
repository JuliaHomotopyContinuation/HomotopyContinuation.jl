import MultivariatePolynomials
const MP = MultivariatePolynomials

import .Problems
import .Systems
import .Homotopies
import .Solving

export solve


"""
    solve(F::Vector{<:MP.AbstractPolynomial}; options...)

Solve the system `F` using a total degree homotopy.

    solve(G, F, start_solutions; options...)

Solve the system `F` by tracking the each solution of
`G` (as provided by `start_solutions`).
"""
function solve end

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
    seed=rand(100:10_000),
    system=Systems.SPSystem,
    homotopy=Homotopies.StraightLineHomotopy,
    kwargs...)
    println("Seed used: $(seed)")
    srand(seed)
    P = Problems.ProjectiveStartTargetProblem(prob, system=system, homotopy=homotopy)
    solve(P, start_solutions; kwargs...)
end

function solve(prob::Problems.AbstractProblem, start_solutions, t₁=1.0, t₀=0.0; threading=true, kwargs...)
    solver = Solving.Solver(prob, start_solutions, t₁, t₀; kwargs...)
    if threading && Threads.nthreads() > 1
        solvers = append!([solver], [deepcopy(solver) for _=2:Threads.nthreads()])
        Solving.solve(solvers, start_solutions)
    else
        Solving.solve(solver, start_solutions)
    end
end
