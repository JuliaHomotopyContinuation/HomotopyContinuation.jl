@testset "solve" begin
    @polyvar x
    @test nfinite(solve([x - 1])) == 1
    F = equations(katsura5())
    @test nfinite(solve(F, threading=false)) == 32
    @test nfinite(solve(F, system=Systems.SPSystem, threading=false)) == 32
    @test nfinite(solve(F, system=Systems.FPSystem, threading=false)) == 32

    @test nfinite(solve(F, homotopy=Homotopies.StraightLineHomotopy)) == 32
    result = solve(F, predictor=Predictors.Euler(), homotopy=Homotopies.StraightLineHomotopy)
    @test nresults(result) == 32
    @test nfinite(solve(F, tol=1e-5)) == 32

    result = solve(F)
    @test nfinite(result) == 32
    @test string(result) isa String

    @test nfinite(solve(F, patch=AffinePatches.RandomPatch())) == 32
    @test nfinite(solve(F, patch=AffinePatches.EmbeddingPatch())) ≤ 32
    @test nfinite(solve(F, patch=AffinePatches.OrthogonalPatch())) == 32


end

@testset "solve - no endgame" begin
    F = equations(katsura5())
    # no endgame
    @test nfinite(solve(F, endgame_start=0.0)) == 32

    @test nfinite(solve(F, endgame_start=0.0, threading=false)) == 32
end

@testset "Path Crossing" begin
    F = equations(katsura5())
    # this will have two crossed paths
    srand(120)
    @test nfinite(solve(F, tol=1e-1, threading=true)) == 32
    srand(120)
    @test nfinite(solve(F, tol=1e-1, threading=false)) == 32
end
