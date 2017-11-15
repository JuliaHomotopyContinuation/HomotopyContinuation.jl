@testset "SphericalPredictorCorrector" begin
    PolyImpl.@polyvar x
    H = StraightLineHomotopy(
        [(x - 4.3im) * (x + (2.1 - 4im))],
        [(x - 2.0) * (x - (2.5+ 4.0im))])

    res = solve(H, [[-2.1 + 4.0im], [4.3im]], SphericalPredictorCorrector(), endgame_start=0.0, apply_gammatrick=false)

    @test res[1].returncode == :isolated
    @test res[2].returncode == :isolated

    if real.(res[1].solution) ≈ [2.5]
        @test real.(res[1].solution) ≈ [2.5]
        @test imag.(res[1].solution) ≈ [4.0]
        @test real(res[2].solution) ≈ [2.0]
    else
        @test real.(res[2].solution) ≈ [2.5]
        @test imag.(res[2].solution) ≈ [4.0]
        @test real(res[1].solution) ≈ [2.0]
    end
end
