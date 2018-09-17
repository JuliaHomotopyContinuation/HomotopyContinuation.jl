@testset "pathtracker_startsolutions" begin
    @polyvar  x y

    tracker, S = pathtracker_startsolutions([x^2-y^2+4, x + y - 3])
    result = PathTracking.track(tracker, first(S), 1.0, 0.0)
    @test result.returncode == :success
end
