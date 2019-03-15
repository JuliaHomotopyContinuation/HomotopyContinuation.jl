function setup_square_corrector_test()
    F = SPSystem(equations(katsura(6)))
    x = rand(Complex{Float64}, 7)
    xnext = copy(x)
    t = rand()
    H = HomotopyWithCache(StraightLineHomotopy(F, F), x, t)
    H, x, xnext, t
end

@testset "Correctors" begin
    @testset "Newton" begin
        H, x, xnext, t = setup_square_corrector_test()

        corrector = NewtonCorrector()
        @test corrector isa NewtonCorrector
        corrector_cache = cache(corrector, H, x, t)
        @test corrector_cache isa HomotopyContinuation.NewtonCorrectorCache

        # check that this doesn't throw
        out = correct!(xnext, corrector, corrector_cache, H, x, t, euclidean_norm, 1e-7, 3)
        @test out isa CorrectorResult
    end
end
