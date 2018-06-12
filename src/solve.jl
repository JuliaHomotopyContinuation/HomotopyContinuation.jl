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
function solve(F::Vector{<:MP.AbstractPolynomial}; seed=randseed(), kwargs...)
    srand(seed)
    TDP = Problems.TotalDegreeProblem(F)
    start_solutions = Utilities.totaldegree_solutions(F)
    solve(TDP, start_solutions, seed; kwargs...)
end

function solve(G::Vector{<:MP.AbstractPolynomial}, F::Vector{<:MP.AbstractPolynomial}, startsolutions; seed=randseed(), kwargs...)
    srand(seed)
    STP = Problems.StartTargetProblem(G, F)

    solve(STP, promote_startsolutions(startsolutions), seed; kwargs...)
end

function solve(F::Vector{<:MP.AbstractPolynomial}, p::Vector{<:MP.AbstractVariable}, a_1::Vector{<:Number}, a_2::Vector{<:Number}, startsolutions; seed=randseed(), homotopy=nothing, kwargs...)
    srand(seed)
    STP = Problems.ParameterProblem(F, p, a_1, a_2)

    @assert length(p) == length(a_1) "Number of parameters must match"
    @assert length(a_1) == length(a_2) "Start and target parameters must have the same length"
    solve(STP, promote_startsolutions(startsolutions), seed; homotopy = Homotopies.ParameterHomotopy, kwargs...)
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
