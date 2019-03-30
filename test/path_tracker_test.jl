@testset "PathTracker" begin
    f = equations(griewank_osborne())

    tracker, startsolutions = pathtracker_startsolutions(f, seed=130793)
    S = collect(startsolutions)
    # This path has the special case that we track towards a non-singular solution
    # at infinity and the valuation doesn't stabilize fast enough
    result = track(tracker, S[3])
    @test result.return_code == PathTrackerStatus.at_infinity
end
