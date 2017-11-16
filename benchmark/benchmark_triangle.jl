using HomotopyContinuation
using DynamicPolynomials
using BenchmarkTools

function triangle1()
    a = 5
    b = 4
    c = 3
    @polyvar sθ cθ
    f = [cθ^2 + sθ^2 - (1.0 + 0im), (a*cθ - b)^2 + (1.0 + 0im) * (a*sθ)^2 - c^2]
    H, s = totaldegree(GeodesicOnTheSphere, f)
    S = solve(H, s)
    @assert (length(isolated_solutions(S)) == 2)
    @assert (length(real_solutions(S)) == 2)
    @benchmark solve($H, $s)
end
println("\n\nRun benchmark 'triangle' with StraightLineHomotopy:")
triangle_benchmark = triangle1()
show(STDOUT, MIME"text/plain"(), triangle_benchmark)
println("\n.")

function triangle2()
    a = 5
    b = 4
    c = 3
    @polyvar sθ cθ
    f = [cθ^2 + sθ^2 - (1.0 + 0im), (a*cθ - b)^2 + (1.0 + 0im) * (a*sθ)^2 - c^2]
    H, s = totaldegree(StraightLineHomotopy, f)
    S = solve(H, s)
    @assert (length(isolated_solutions(S)) == 2)
    @assert (length(real_solutions(S)) == 2)
    @benchmark solve($H, $s)
end
println("\n\nRun benchmark 'triangle' with GeodesicOnTheSphere:")
triangle_benchmark = triangle2()
show(STDOUT, MIME"text/plain"(), triangle_benchmark)
println("\n.")
