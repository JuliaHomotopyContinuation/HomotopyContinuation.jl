function setup_square_corrector_test()
    F = SPSystem(equations(katsura(6)))
    x = rand(Complex{Float64}, 7)
    xnext = copy(x)
    t = rand()
    H = HomotopyWithCache(StraightLineHomotopy(F, F), x, t)
    Jac = HC.Jacobian(HC.jacobian(H, x, t))
    H, x, xnext, t, Jac
end


function setup_overdetermined_corrector_test()
    @polyvar x y
    F = SPSystem([x^2+y, y^2-3x*y, y+x+3, y+2x-5])
    x = rand(Complex{Float64}, 2)
    xnext = copy(x)
    t = rand()
    H = HomotopyWithCache(StraightLineHomotopy(F, F), x, t)
    Jac = HC.Jacobian(HC.jacobian(H, x, t))
    H, x, xnext, t, Jac
end


@testset "Correctors" begin
    @testset "Newton - square" begin
        H, x, xnext, t, Jac = setup_square_corrector_test()

        corrector = NewtonCorrector()
        @test corrector isa NewtonCorrector
        corrector_cache = cache(corrector, H, x, t)
        @test corrector_cache isa HomotopyContinuation.NewtonCorrectorCache

        # check that this doesn't throw
        out = correct!(xnext, corrector, corrector_cache, H, x, t, euclidean_norm, Jac, 1e-7, 3,
        HC.StepSizeModel())
        @test out isa CorrectorResult

        out = correct!(xnext, corrector, corrector_cache, H, x, t,
                    euclidean_norm, Jac, 1e-7, 3, HC.StepSizeModel(); update_jacobian_infos=true)
        @test out isa CorrectorResult
    end

    @testset "Newton - overdetermined" begin
        H, x, xnext, t, Jac = setup_overdetermined_corrector_test()

        corrector = NewtonCorrector()
        @test corrector isa NewtonCorrector
        corrector_cache = cache(corrector, H, x, t)
        @test corrector_cache isa HomotopyContinuation.NewtonCorrectorCache

        # check that this doesn't throw
        out = correct!(xnext, corrector, corrector_cache, H, x, t, euclidean_norm, Jac, 1e-7, 3, HC.StepSizeModel())
        @test out isa CorrectorResult

        out = correct!(xnext, corrector, corrector_cache, H, x, t,
                    euclidean_norm, Jac, 1e-7, 3, HC.StepSizeModel(); update_jacobian_infos=true)
        @test out isa CorrectorResult
    end
end
