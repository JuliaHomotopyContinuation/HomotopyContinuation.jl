@testset "Projective Tracker" begin
    const PathTrackers = HomotopyContinuation.PathTrackers
    p1 = Problems.TotalDegreeProblem(equations(katsura7()))
    P = Problems.ProjectiveStartTargetProblem(p1)

    sols = Utilities.totaldegree_solutions(p1.system) |> collect

    tracker = PathTracking.PathTracker(P.homotopy, Problems.embed(P, sols[123]), 1.0, 0.0)
    @test tracker isa PathTracking.PathTracker{<:PathTrackers.Projective,
        <:PathTrackers.ProjectiveState, <:PathTrackers.ProjectiveCache}

    @test PathTracking.current_status(tracker) == :ok

    tracker = PathTracking.PathTracker(P.homotopy, collect(1:9), 1.0, 0.0)
    @test PathTracking.current_x(tracker) isa Vector{Complex{Float64}}
    @test PathTracking.current_status(tracker) == :invalid_startvalue

    @test_throws ErrorException PathTracking.PathTracker(P.homotopy, ones(Int, 8), 1.0, 0.0)
end
