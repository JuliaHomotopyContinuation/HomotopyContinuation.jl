@testset "LinearSystem" begin
    A = rand(3, 3)
    b = rand(3)

    PolyImpl.@polyvar x y z
    F = A * [x, y, z] - b

    P1 = Problems.TotalDegreeProblem(F)
    P = Problems.ProjectiveStartTargetProblem(P1)

    sols = Utilities.totaldegree_solutions(F) |> collect

    start_sols = Problems.embed.(P, sols)

    iter = PathTrackers.pathtracker(P.homotopy, start_sols[1], 1.0, 0.0)
    PathTrackers.track!(iter)
    x = iter.state.x
    @test norm(x[2:end] / x[1] - A \ b) < 1e-8
end
