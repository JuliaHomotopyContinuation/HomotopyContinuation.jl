using HomotopyContinuation
using DynamicPolynomials
using BenchmarkTools

function triangle()
    a = 5
    b = 4
    c = 3

    @polyvar sθ cθ
    f1 = cθ^2 + sθ^2 - (1.0 + 0im)
    f2 = (a*cθ - b)^2 + (1.0 + 0im) * (a*sθ)^2 - c^2

    # lets use a total degree straight line homotopy
    H, start_solutions = totaldegree(StraightLineHomotopy, [f1, f2])

    #now we can solve
    results = solve(H, start_solutions, SphericalPredictorCorrector())
    solutions = solution.(filter(issuccessfull, results))
    @assert (length(solutions) == 2)

    @benchmark solve($H, $start_solutions, $SphericalPredictorCorrector())
end

println("\n\nRun polysystem benchmark:")
triangle_benchmark = triangle()
show(STDOUT, MIME"text/plain"(), triangle_benchmark)

println("\nDone.")
