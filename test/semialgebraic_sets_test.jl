@testset "SemialgebraicSets" begin
    solver = SemialgebraicSetsHCSolver(; compile = false)
    @polyvar x y
    V = SemialgebraicSets.@set x^2 == 1 && y^2 == 2 solver
    S = collect(V)
    S = sort!(map(s -> round.(s, digits = 2), S))
    @test S == [[-1.0, -1.41], [-1.0, 1.41], [1.0, -1.41], [1.0, 1.41]]
end
