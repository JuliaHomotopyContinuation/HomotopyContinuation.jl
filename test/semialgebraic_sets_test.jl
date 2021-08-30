@testset "SemialgebraicSets" begin
    solver = SemialgebraicSetsHCSolver(; compile = false)
    @polyvar x y
    V = SemialgebraicSets.@set x^2 == 1 && y^2 == 2 solver
    S = collect(V)
    S = sort!(map(s -> round.(s, digits = 2), S))
    @test S == [[-1.0, -1.41], [-1.0, 1.41], [1.0, -1.41], [1.0, 1.41]]

    V = SemialgebraicSets.@set x == 1 && x == 1 solver
    S = collect(V)
    S = sort!(map(s -> round.(s, digits = 2), S))
    @test S == [[1.0]]

    # Inspired from https://jump.dev/SumOfSquares.jl/v0.4.6/generated/Polynomial%20Optimization/
    V = SemialgebraicSets.algebraicset([
        -x - y + 1,
        -x*y + y^2 - y,
        -x^2 + y^2 - 2y + 1],
        solver,
    )
    S = collect(V)
    S = sort!(map(s -> round.(s, digits = 2), S))
    @test S == [[0.0, 1.0], [1.0, 0.0]]
end
