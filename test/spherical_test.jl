@testset "Spherical" begin

    @testset "simplest system" begin
        PolyImpl.@polyvar x
        H = StraightLineHomotopy([(x - 4.3im) * (x + (2.1 - 4im))], [(x - 2.0) * (x - (2.5+ 4.0im))])

        res = solve(H, [[-2.1 + 4.0im], [4.3im]], Spherical(), report_progress=false)

        @test res[1].retcode == :Success
        @test res[2].retcode == :Success

        @test real(res[1].solution) ≈ [2.0]
        @test real(res[2].solution) ≈ [2.5]
        @test imag(res[2].solution) ≈ [4.0]
    end

end
