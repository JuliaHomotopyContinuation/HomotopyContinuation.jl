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
    # intilizae tracker with random start and target to see whether it is properly resetetted
    tracker = PathTracking.PathTracker(P.homotopy, s, rand(), rand())

    result = PathTracking.track(tracker, s, 1.0, 0.0)
    @test result.returncode == :success
    @test result.t == 0.0
    x = result.x
    @test norm(x[2:end] / x[1] - A \ b) < 1e-12

    x_inter = copy(s)
    retcode = PathTracking.track!(x_inter, tracker, s, 1.0, 0.1)
    @test retcode == :success
    x_final = zero(x_inter)
    retcode = PathTracking.track!(x_final, tracker, x_inter, 0.1, 0.0, false)
    @test retcode == :success
    tracker
    @test PathTracking.current_iter(tracker) < 3
    x = PathTracking.current_x(tracker)
    @test norm(x[2:end] / x[1] - A \ b) < 1e-12
end
