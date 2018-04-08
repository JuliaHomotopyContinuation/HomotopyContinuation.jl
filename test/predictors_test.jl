@testset "Predictors" begin
    F = Systems.SPSystem(equations(katsura6()))
    x = rand(Complex{Float64}, 7)
    xnext = copy(x)
    t = rand()
    H = NewHomotopies.HomotopyWithCache(NewHomotopies.StraightLineHomotopy(F, F), x, t)

    predictor = Predictors.Euler()
    @test predictor isa Predictors.Euler
    predictor_cache = Predictors.cache(predictor, H, x, t)
    @test predictor_cache isa Predictors.EulerCache

    # check that this doesn't throw
    Predictors.predict!(xnext, predictor, predictor_cache, H, x, t, 0.05)
end
