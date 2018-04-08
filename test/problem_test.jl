@testset "Problems" begin
    F = equations(katsura5())
    G = equations(cyclic6())

    P1 = Problems.StartTargetProblem(G, F)
    P2 = Problems.TotalDegreeProblem(F)

    PP1 = Problems.ProjectiveStartTargetProblem(P1)
    @test PP1 isa Problems.ProjectiveStartTargetProblem
    PP2 = Problems.ProjectiveStartTargetProblem(P2)
    @test PP2 isa Problems.ProjectiveStartTargetProblem

    @test Problems.homotopy(PP2) isa NewHomotopies.StraightLineHomotopy

    @test Problems.homogenization_strategy(PP2) == Problems.DefaultHomogenization()

    PolyImpl.@polyvar x y

    P = Problems.StartTargetProblem([x^2+y^2, x^4+y^4], [x^3,y^3])
end
