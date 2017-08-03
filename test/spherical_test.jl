using HomotopyContinuation

@testset "spherical" begin

    @testset "simplest system" begin
        x = MPoly.generator(Complex128, :x)
        f = (x - 2.0) * (x - (2.5+ 4.0im))
        g = (x - 4.3im) * (x + (2.1 - 4im))

        F = MPoly.system([f])
        G = MPoly.system([g])

        H = StraightLineHomotopy(G, F)

        res = solve(H, [[-2.1 + 4.0im], [4.3im]], PredictorCorrector.Spherical(), report_progress=false)

        @test res[1].retcode == :Success
        @test res[2].retcode == :Success

        @test real(res[1].solution) ≈ [2.0]
        @test real(res[2].solution) ≈ [2.5]
        @test imag(res[2].solution) ≈ [4.0]
    end

end