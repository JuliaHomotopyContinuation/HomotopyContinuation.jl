import LinearAlgebra, ProjectiveVectors
function setup_prediction_test()
    F = Systems.SPSystem(equations(katsura(6)))
    x = ProjectiveVectors.PVector(rand(Complex{Float64}, 7))
    xnext = copy(x)
    t = rand()
    H = Homotopies.HomotopyWithCache(Homotopies.StraightLineHomotopy(F, F), x, t)
    J = Homotopies.jacobian(H, x, t)
    fac = Utilities.factorization(J)
    ẋ = fac \ -Homotopies.dt(H, x, t)
    H, x, xnext, t, ẋ, fac
end

function setup_overdetermined_prediction_test()
    @polyvar x y
    F = Systems.SPSystem([x^2+y, y^2-3x*y, y+x+3, y+2x-5])
    x = ProjectiveVectors.PVector(rand(Complex{Float64}, 2))
    xnext = copy(x)
    t = rand()
    H = Homotopies.HomotopyWithCache(Homotopies.StraightLineHomotopy(F, F), x, t)
    J = Homotopies.jacobian(H, x, t)
    fac = Utilities.factorization(J)
    ẋ = fac \ -Homotopies.dt(H, x, t)
    H, x, xnext, t, ẋ, fac
end

function test_predictor(PREDICTOR, CACHE)
    H, x, xnext, t, ẋ, fac = setup_prediction_test()

    predictor = PREDICTOR()
    @test predictor isa PREDICTOR
    predictor_cache = Predictors.cache(predictor, H, x, ẋ, t)
    Predictors.setup!(predictor_cache, H, x, ẋ, t, fac)
    @test predictor_cache isa CACHE
    # check that this doesn't throw
    Predictors.update!(predictor_cache, H, x, ẋ, t, fac)
    @test_nowarn Predictors.predict!(xnext, predictor, predictor_cache, H, x, t, 0.05, ẋ)


    H, x, xnext, t, ẋ, fac = setup_overdetermined_prediction_test()
    predictor = PREDICTOR()
    @test predictor isa PREDICTOR
    predictor_cache = Predictors.cache(predictor, H, x, ẋ, t)
    Predictors.setup!(predictor_cache, H, x, ẋ, t, fac)
    @test predictor_cache isa CACHE
    # check that this doesn't throw
    Predictors.update!(predictor_cache, H, x, ẋ, t, fac)
    @test_nowarn Predictors.predict!(xnext, predictor, predictor_cache, H, x, t, 0.05, ẋ)
end


@testset "Predictors" begin

    @testset "NullPredictor" begin
        test_predictor(Predictors.NullPredictor, Predictors.NullPredictorCache)
    end

    @testset "Euler" begin
        test_predictor(Predictors.Euler, Predictors.EulerCache)
    end

    @testset "Heun" begin
        test_predictor(Predictors.Heun, Predictors.HeunCache)
    end

    @testset "Midpoint" begin
        test_predictor(Predictors.Midpoint, Predictors.MidpointCache)
    end

    @testset "Ralston" begin
        test_predictor(Predictors.Ralston, Predictors.RalstonCache)
    end

    @testset "RK3" begin
        test_predictor(Predictors.RK3, Predictors.RK3Cache)
    end

    @testset "Pade21" begin
        test_predictor(Predictors.Pade21, Predictors.Pade21Cache)
    end

    @testset "RK4" begin
        test_predictor(Predictors.RK4, Predictors.RK4Cache)
    end
end
