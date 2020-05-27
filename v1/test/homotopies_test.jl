@testset "Homotopies" begin
    @testset "StraightLineHomotopy" begin
        F = SPSystem(equations(katsura(5)))
        G = SPSystem(equations(cyclic(6)))
        H = StraightLineHomotopy(F, G)
        @test H isa AbstractHomotopy
        @test size(H) == (6, 6)
        @test HC.gamma(H) isa ComplexF64
        HomotopyContinuation.homotopy_interface_test(H)
    end

    @testset "HomotopyWithCache" begin
        x = rand(Complex{Float64}, 6)
        t = rand()
        F = SPSystem(equations(katsura(5)))
        G = SPSystem(equations(cyclic(6)))
        H = HomotopyWithCache(StraightLineHomotopy(F, G), x, t)
        @test H isa AbstractHomotopy
        @test size(H) == (6, 6)

        m, n = size(H)
        u = zeros(Complex{Float64}, m)
        U = zeros(Complex{Float64}, m, n)

        evaluate!(u, H, x, t)
        @test evaluate(H, x, t) ≈ u

        dt!(u, H, x, t)
        @test dt(H, x, t) ≈ u

        jacobian!(U, H, x, t)
        @test jacobian(H, x, t) ≈ U

        evaluate_and_jacobian!(u, U, H, x, t)
        @test all(evaluate_and_jacobian(H, x, t) .≈ (u, U))
        @test all((u, U) .≈ (evaluate(H, x, t), jacobian(H, x, t)))

        jacobian_and_dt!(U, u, H, x, t)
        @test all(jacobian_and_dt(H, x, t) .≈ (U, u))
        @test all((U, u) .≈ (jacobian(H, x, t), dt(H, x, t)))
    end

    @testset "ParameterHomotopy" begin
        F = equations(katsura(5))
        vars = variables(F)
        H = ParameterHomotopy(F, vars[5:6], p₁ = rand(2), p₀ = rand(2))
        @test H isa AbstractHomotopy
        @test size(H) == (6, 4)
        @test size(H, 1) == 6
        @test size(H, 2) == 4
        @test length(H) == 6
        @test HC.nvariables(H) == 4

        HomotopyContinuation.homotopy_interface_test(H)

        H = ParameterHomotopy(
            F,
            vars[5:6],
            p₁ = rand(2),
            p₀ = rand(2),
            γ₁ = randn(ComplexF64),
            γ₀ = randn(ComplexF64),
        )
        @test H isa ParameterHomotopy
        HomotopyContinuation.homotopy_interface_test(H)
    end

    @testset "FixedPointHomotopy" begin
        F = SPSystem(equations(katsura(5)))
        H = FixedPointHomotopy(F, rand(ComplexF64, 6))
        @test H isa AbstractHomotopy
        @test size(H) == (6, 6)

        HomotopyContinuation.homotopy_interface_test(H)
    end

    @testset "PatchedHomotopy" begin
        F = SPSystem(equations(katsura(5)))
        G = SPSystem(equations(cyclic(6)))
        x = ProjectiveVectors.embed(rand(Complex{Float64}, 5))
        H = PatchedHomotopy(StraightLineHomotopy(F, G), OrthogonalPatch(), x)
        @test H isa AbstractHomotopy
        @test size(H) == (7, 6)

        HomotopyContinuation.homotopy_interface_test(H, x)
    end

    @testset "Coefficient Homotopy" begin
        E = [[2 1 0; 0 0 0], [1 0; 1 0]]
        start = [[1.0 + 0im, -3.0, 2.0], [2.0 + 0im, -2.0]]
        target = [[2.0 + 0im, -2.0, 5.0], [3.0 + 0im, -1.0]]
        H = CoefficientHomotopy(E, start, target)
        x = randn(ComplexF64, 2)
        HomotopyContinuation.homotopy_interface_test(H, x)
    end

    @testset "ConstantHomotopy" begin
        F = SPSystem(equations(katsura(5)))
        H = ConstantHomotopy(F)
        @test H isa AbstractHomotopy
        @test size(H) == (6, 6)

        HomotopyContinuation.homotopy_interface_test(H)
    end
end
