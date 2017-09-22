@testset "Spherical" begin

    @testset "simplest system" begin
        PolyImpl.@polyvar x
        H = StraightLineHomotopy([(x - 4.3im) * (x + (2.1 - 4im))], [(x - 2.0) * (x - (2.5+ 4.0im))])

        res = solve(H, [[-2.1 + 4.0im], [4.3im]], SphericalPredictorCorrector(), endgame_start=0.0)

        @test issuccessfull(res[1])
        @test issuccessfull(res[2])

        @test real(solution(res[1])) ≈ [2.0]
        @test real(solution(res[2])) ≈ [2.5]
        @test imag(solution(res[2])) ≈ [4.0]
    end

end
