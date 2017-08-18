@testset "cauchyendgame" begin
     PolyImpl.@polyvar x
     F = PolySystem([(x - 2.0)^4])
     G, solutions = totaldegree(F)

     H = StraightLineHomotopy(G,F)

     results = solve(H, solutions)
     for result in results
          @test norm([2.0] - result.solution) ≈ 0.0 atol=1e-3
     end
end
