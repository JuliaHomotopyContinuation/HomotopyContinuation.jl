@testset "LinearSystem" begin
    A = rand(3, 3)
    b = rand(3)

    PolyImpl.@polyvar x y z
    F = A * [x, y, z] - b

    P1 = Problems.TotalDegreeProblem(F)
    P = Problems.ProjectiveStartTargetProblem(P1)

    sols = Utilities.totaldegree_solutions(F) |> collect

    start_sols = Problems.embed.(P, sols)

    s = start_sols[1]
    tracker = PathTracking.PathTracker(P.homotopy, s, 1.0, 0.0)
    PathTracking.track!(tracker)
    x = PathTracking.current_x(tracker)
    @test PathTracking.current_t(tracker) == 0.0
    @test norm(x[2:end] / x[1] - A \ b) < 1e-8

    PathTracking.track!(tracker, s, 1.0, 0.0)
    x = PathTracking.current_x(tracker)
    @test norm(x[2:end] / x[1] - A \ b) < 1e-8
end
