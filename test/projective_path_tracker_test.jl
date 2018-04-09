@testset "Projective Tracker" begin
    p1 = Problems.TotalDegreeProblem(equations(katsura7()))
    P = Problems.ProjectiveStartTargetProblem(p1)

    sols = Utilities.totaldegree_solutions(p1.system) |> collect

    iter = PathTrackers.pathtracker(P.homotopy, Problems.embed(P, sols[123]), 1.0, 0.0)
    @test iter isa PathTrackers.TrackerIterator{<:PathTrackers.ProjectiveTracker,
        <:PathTrackers.ProjectiveState, <:PathTrackers.ProjectiveTrackerCache}

    @test iter.state.status == :default

    iter = PathTrackers.pathtracker(P.homotopy, collect(1:9), 1.0, 0.0)
    @test iter.state.x isa Vector{Complex{Float64}}
    @test iter.state.status == :invalid_startvalue

    @test_throws AssertionError PathTrackers.pathtracker(P.homotopy, ones(Int, 8), 1.0, 0.0)
end
