import LinearAlgebra, ProjectiveVectors
function setup_prediction_test()
    F = SPSystem(equations(katsura(6)))
    x = ProjectiveVectors.PVector(rand(Complex{Float64}, 7))
    xnext = copy(x)
    t = rand()
    H = HomotopyWithCache(StraightLineHomotopy(F, F), x, t)
    J = jacobian(H, x, t)
    fac = HC.factorization(J)
    ẋ = fac \ -dt(H, x, t)
    H, x, xnext, t, ẋ, fac
end

function setup_overdetermined_prediction_test()
    @polyvar x y
    F = SPSystem([x^2+y, y^2-3x*y, y+x+3, y+2x-5])
    x = ProjectiveVectors.PVector(rand(Complex{Float64}, 2))
    xnext = copy(x)
    t = rand()
    H = HomotopyWithCache(StraightLineHomotopy(F, F), x, t)
    J = jacobian(H, x, t)
    fac = HC.factorization(J)
    ẋ = fac \ -dt(H, x, t)
    H, x, xnext, t, ẋ, fac
end

function test_predictor(PREDICTOR, CACHE)
    H, x, xnext, t, ẋ, fac = setup_prediction_test()

    predictor = PREDICTOR()
    @test predictor isa PREDICTOR
    predictor_cache = cache(predictor, H, x, ẋ, t)
    setup!(predictor_cache, H, x, ẋ, t, fac)
    @test predictor_cache isa CACHE
    # check that this doesn't throw
    update!(predictor_cache, H, x, ẋ, t, fac)
    @test_nowarn predict!(xnext, predictor, predictor_cache, H, x, t, 0.05, ẋ)


    H, x, xnext, t, ẋ, fac = setup_overdetermined_prediction_test()
    predictor = PREDICTOR()
    @test predictor isa PREDICTOR
    predictor_cache = cache(predictor, H, x, ẋ, t)
    setup!(predictor_cache, H, x, ẋ, t, fac)
    @test predictor_cache isa CACHE
    # check that this doesn't throw
    update!(predictor_cache, H, x, ẋ, t, fac)
    @test_nowarn predict!(xnext, predictor, predictor_cache, H, x, t, 0.05, ẋ)
end


@testset "Predictors" begin

    @testset "NullPredictor" begin
        test_predictor(NullPredictor, HomotopyContinuation.NullPredictorCache)
    end

    @testset "Euler" begin
        test_predictor(Euler, HomotopyContinuation.EulerCache)
    end

    @testset "Heun" begin
        test_predictor(Heun, HomotopyContinuation.HeunCache)
    end

    @testset "Midpoint" begin
        test_predictor(Midpoint, HomotopyContinuation.MidpointCache)
    end

    @testset "Ralston" begin
        test_predictor(Ralston, HomotopyContinuation.RalstonCache)
    end

    @testset "RK3" begin
        test_predictor(RK3, HomotopyContinuation.RK3Cache)
    end

    @testset "Pade21" begin
        test_predictor(Pade21, HomotopyContinuation.Pade21Cache)
    end

    @testset "RK4" begin
        test_predictor(RK4, HomotopyContinuation.RK4Cache)
    end
end
