@testset "Correctors" begin
    F = Systems.SPSystem(equations(katsura6()))
    x = rand(Complex{Float64}, 7)
    xnext = copy(x)
    t = rand()
    H = Homotopies.HomotopyWithCache(Homotopies.StraightLineHomotopy(F, F), x, t)

    corrector = Correctors.Newton()
    @test corrector isa Correctors.Newton
    corrector_cache = Correctors.cache(corrector, H, x, t)
    @test corrector_cache isa Correctors.NewtonCache

    # check that this doesn't throw
    out = Correctors.correct!(xnext, corrector, corrector_cache, H, x, t, 1e-7, 3)
    @test out isa Correctors.Result
end
