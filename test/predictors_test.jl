function setup_prediction_test()
    F = Systems.SPSystem(equations(katsura(6)))
    x = rand(Complex{Float64}, 7)
    xnext = copy(x)
    t = rand()
    H = Homotopies.HomotopyWithCache(Homotopies.StraightLineHomotopy(F, F), x, t)
    H, x, xnext, t
end

@testset "Predictors.NullPredictor" begin
    H, x, xnext, t = setup_prediction_test()

    predictor = Predictors.NullPredictor()
    @test predictor isa Predictors.NullPredictor
    predictor_cache = Predictors.cache(predictor, H, x, t)
    @test predictor_cache isa Predictors.NullPredictorCache

    # check that this doesn't throw
    Predictors.predict!(xnext, predictor, predictor_cache, H, x, t, 0.05)
    @test xnext == x
end

@testset "Predictors.Euler" begin
    H, x, xnext, t = setup_prediction_test()

    predictor = Predictors.Euler()
    @test predictor isa Predictors.Euler
    predictor_cache = Predictors.cache(predictor, H, x, t)
    @test predictor_cache isa Predictors.EulerCache

    # check that this doesn't throw
    Predictors.predict!(xnext, predictor, predictor_cache, H, x, t, 0.05)
end

@testset "Predictors.RK4" begin
    H, x, xnext, t = setup_prediction_test()

    predictor = Predictors.RK4()
    @test predictor isa Predictors.RK4
    predictor_cache = Predictors.cache(predictor, H, x, t)
    @test predictor_cache isa Predictors.RK4Cache

    # check that this doesn't throw
    Predictors.predict!(xnext, predictor, predictor_cache, H, x, t, 0.05)
end
