@testset "Projective Tracker" begin
    p1 = Problems.TotalDegreeProblem(equations(katsura7()))
    P = Problems.ProjectiveStartTargetProblem(p1)

    sols = Utilities.totaldegree_solutions(p1.system) |> collect

    tracker = PathTrackers.pathtracker(P)
    @test tracker isa PathTrackers.ProjectiveTracker
    iter = PathTrackers.iterator(tracker, Problems.embed(P, sols[123]), 1.0, 0.0)
    @test iter isa PathTrackers.TrackerIterator{<:PathTrackers.ProjectiveTracker,
        <:PathTrackers.ProjectiveState, <:PathTrackers.ProjectiveTrackerCache}
end
