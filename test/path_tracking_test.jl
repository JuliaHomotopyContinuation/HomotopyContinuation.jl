
@testset "PathTracking" begin
    F = equations(katsura5())

    TDP = Problems.TotalDegreeProblem(F)
    P = Problems.ProjectiveStartTargetProblem(TDP)
    start_sols = Utilities.totaldegree_solutions(F) |> collect

    # test construction
    t1 = PathTracking.PathTracker(P, first(start_sols), 1.0, 0.1)

    @test t1 isa PathTracking.PathTracker
    @test length(t1.x) == 7

    t2 = PathTracking.PathTracker(P, Problems.embed(P, first(start_sols)), 1.0, 0.1)
    @test t2 isa PathTracking.PathTracker
    @test length(t2.x) == 7

    t3 = PathTracking.PathTracker(P, first(start_sols), 1.0, 0.1, predictor=Predictors.Euler())
    @test t3.predictor_corrector.predictor isa Predictors.Euler


    PathTracking.setup!(t1, Problems.embed(P, start_sols[2]), 1.0, 0.4)
    @test PathTracking.currstatus(t1) == :ok
    @test PathTracking.currt(t1) == 1.0

    PathTracking.setup!(t1, Problems.embed(P, start_sols[2]), 0.5, 0.4)
    @test PathTracking.currstatus(t1) == :invalid_startvalue
    @test PathTracking.currt(t1) == 0.5

    R = PathTracking.track(t1, Problems.embed(P, first(start_sols)), 1.0, 0.0)
    @test R isa PathTracking.PathTrackerResult
    @test R.returncode == :success
    @test R.res < 1e-7

    out = Problems.embed(P, first(start_sols))
    retcode = PathTracking.track!(out, t1, Problems.embed(P, first(start_sols)), 1.0, 0.0)
    @test retcode == :success
    @test out == R.x
end

@testset "LinearSystem" begin
    A = rand(3, 3)
    b = rand(3)

    @polyvar x y z
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
    retcode = PathTracking.track!(x_final, tracker, x_inter, 0.1, 0.0)
    @test retcode == :success
    tracker
    @test PathTracking.curriters(tracker) < 3
    x = PathTracking.currx(tracker)
    @test norm(x[2:end] / x[1] - A \ b) < 1e-12
end

@testset "fixedpatch" begin

    F = equations(katsura5())
    TDP = Problems.TotalDegreeProblem(F)
    P = Problems.ProjectiveStartTargetProblem(TDP)
    start_sols = Utilities.totaldegree_solutions(F) |> collect

    # test construction
    tracker = PathTracking.PathTracker(P, first(start_sols), 1.0, 0.1)
    r1 = PathTracking.track(tracker, Problems.embed(P, start_sols[1]), 1.0, 0.1)

    PathTracking.fixpatch!(tracker, true)
    PathTracking.track(tracker, r1.x, 0.1, 0.0)

    @test all(r1.x .≈ PathTracking.patch(tracker.cache))
end
