@testset "Pathtracker" begin
    H, s = randomhomotopy(StraightLineHomotopy, 4)

    tracker = Pathtracker(H, SphericalPredictorCorrector(), first(s))
    track!(tracker)
    result = PathtrackerResult(tracker)
    @test result.retcode == :success
end
