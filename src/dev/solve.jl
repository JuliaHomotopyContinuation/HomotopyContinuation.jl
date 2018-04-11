import .Problems
import .Systems
import .NewHomotopies
import .PathTracking
import MultivariatePolynomials
const MP = MultivariatePolynomials

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

function solve(prob::Problems.ProjectiveStartTargetProblem, start_solutions; kwargs...)
    xs = map(x -> Problems.embed(prob, x), start_solutions)
    @assert length(xs) > 1 "start solutions are empty"
    t₁ = 1.0
    t₀ = 0.0
    tracker = PathTracking.PathTracker(prob.homotopy, first(xs), t₁, t₀)

    PathTracking.track(tracker, xs, t₁, t₀)
end
