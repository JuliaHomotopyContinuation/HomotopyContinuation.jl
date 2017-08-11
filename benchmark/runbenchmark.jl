using HomotopyContinuation
using TypedPolynomials
using BenchmarkTools

function triangle()
    # lets fix side lengths
    a = 5
    b = 4
    c = 3;

    #create the variables
    @polyvar sθ cθ x0
    # now construct the system
    # In the future we will be able to avoid to explicitly write complex coefficients
    f1 = cθ^2 + sθ^2 - (1.0 + 0im)
    f2 = (a*cθ - b)^2 + (1.0 + 0im) * (a*sθ)^2 - c^2
    F = [f1, f2]

    # lets use a total degree homotopy
    G, start_solutions = total_degree(F)

    # we will use a simple straight line homotopy
    H = StraightLineHomotopy(G, F)

    #now we can solve
    results = solve(H, start_solutions, PredictorCorrector.Spherical(x0), report_progress=false)

    #lets ignore the solutions at infinity
    map(x -> x.solution, filter(x -> x.retcode == :Success, results))
end

println("--------------------")
println("Triangle:")
triangle_benchmark = @benchmark triangle() seconds=10
show(STDOUT, MIME"text/plain"(), triangle_benchmark)