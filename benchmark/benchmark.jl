using HomotopyContinuation
import TypedPolynomials
using BenchmarkTools

function triangle()
    a = 5
    b = 4
    c = 3;

    TypedPolynomials.@polyvar sθ cθ
    f1 = cθ^2 + sθ^2 - (1.0 + 0im)
    f2 = (a*cθ - b)^2 + (1.0 + 0im) * (a*sθ)^2 - c^2
    F = PolySystem([f1, f2])

    # lets use a total degree homotopy
    G, start_solutions = totaldegree(F)

    # we will use a simple straight line homotopy
    H = StraightLineHomotopy(G, F)

    algorithm = PredictorCorrector.Spherical()
    #now we can solve
    results = solve(H, start_solutions, algorithm, report_progress=false)
    results = map(x -> x.solution, filter(x -> x.retcode == :Success, results))
    @assert (length(results) == 2)

    @benchmark solve($H, $start_solutions, $algorithm, report_progress=false)
end

println("\n\nRun polysystem benchmark:")
triangle_benchmark = triangle()
show(STDOUT, MIME"text/plain"(), triangle_benchmark)

println("\nDone.")
