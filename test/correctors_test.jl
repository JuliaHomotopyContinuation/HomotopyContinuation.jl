function setup_square_corrector_test()
    F = Systems.SPSystem(equations(katsura(6)))
    x = rand(Complex{Float64}, 7)
    xnext = copy(x)
    t = rand()
    H = Homotopies.HomotopyWithCache(Homotopies.StraightLineHomotopy(F, F), x, t)
    H, x, xnext, t
end

@testset "Correctors" begin
    @testset "Newton" begin
        H, x, xnext, t = setup_square_corrector_test()

        corrector = Correctors.Newton()
        @test corrector isa Correctors.Newton
        corrector_cache = Correctors.cache(corrector, H, x, t)
        @test corrector_cache isa Correctors.NewtonCache

        # check that this doesn't throw
        out = Correctors.correct!(xnext, corrector, corrector_cache, H, x, t, tol=1e-7, maxiters=3)
        @test out isa Correctors.Result
    end
end
