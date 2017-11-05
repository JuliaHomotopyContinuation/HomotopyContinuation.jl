@testset "CauchyEndgame" begin
     PolyImpl.@polyvar x
     H, solutions = totaldegree(StraightLineHomotopy{Complex128}, [(x - 2.0)^4])


     results = solve(H, solutions, SphericalPredictorCorrector(), CauchyEndgame())

     for result in results
          @test norm([2.0] - result.solution) ≈ 0.0 atol=1e-3
     end
end
