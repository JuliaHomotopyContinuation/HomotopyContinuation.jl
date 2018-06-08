import MultivariatePolynomials
const MP = MultivariatePolynomials

import .Problems
import .Systems
import .Homotopies
import .Solving

using .Utilities

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
function solve(F::Vector{<:MP.AbstractPolynomial}; seed=randseed(), kwargs...)
    srand(seed)
    F = filter(f -> !iszero(f), F)
    checkfinite_dimensional(F)
    TDP = Problems.TotalDegreeProblem(F)
    start_solutions = Utilities.totaldegree_solutions(F)
    solve(TDP, start_solutions, seed; kwargs...)
end

function solve(G::Vector{<:MP.AbstractPolynomial}, F::Vector{<:MP.AbstractPolynomial}, startsolutions; seed=randseed(), kwargs...)
    srand(seed)
    @assert length(G) == length(F)
    checkfinite_dimensional(F)
    STP = Problems.StartTargetProblem(G, F)

    solve(STP, promote_startsolutions(startsolutions), seed; kwargs...)
end

function checkfinite_dimensional(F::Vector{<:MP.AbstractPolynomial})
    N = MP.nvariables(F)
    n = length(F)
    # square system and each polynomial is non-zero
    if n ≥ N ||
       n == N - 1 && ishomogenous(F)
        return
    end
    throw(AssertionError("The input system will not result in a finite number of solutions."))
end

promote_startsolutions(xs::Vector{Vector{Complex128}}) = xs
function promote_startsolutions(xs::Vector{<:AbstractVector{<:Number}})
    PT = promote_type(typeof(xs[1][1]), Complex{Float64})
    map(s -> convert.(PT, s), xs)
end

randseed() = rand(1_000:1_000_000)

# Internal
function solve(prob::Problems.AbstractDynamicProblem, start_solutions, seed;
    system=Systems.FPSystem,
    homotopy=Homotopies.StraightLineHomotopy,
    kwargs...)

    P = Problems.ProjectiveStartTargetProblem(prob, system=system, homotopy=homotopy)
    solve(P, start_solutions, seed; kwargs...)
end

function solve(prob::Problems.AbstractProblem, start_solutions, seed; threading=true, kwargs...)
    solver = Solving.Solver(prob, start_solutions, 1.0, 0.0, seed; kwargs...)
    if threading && Threads.nthreads() > 1
        solvers = append!([solver], [deepcopy(solver) for _=2:Threads.nthreads()])
        Solving.solve(solvers, start_solutions)
    else
        Solving.solve(solver, start_solutions)
    end
end
