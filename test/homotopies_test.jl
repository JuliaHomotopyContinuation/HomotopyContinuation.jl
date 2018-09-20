@testset "Homotopies" begin
    @testset "StraightLineHomotopy" begin
        F = Systems.SPSystem(equations(katsura(5)))
        G = Systems.SPSystem(equations(cyclic(6)))
        H = Homotopies.StraightLineHomotopy(F, G)
        @test H isa Homotopies.AbstractHomotopy
        @test size(H) == (6, 6)
        @test Homotopies.gamma(H) isa ComplexF64
        InterfaceTest.homotopy(H)
    end

    @testset "HomotopyWithCache" begin
        x = rand(Complex{Float64}, 6)
        t = rand()
        F = Systems.SPSystem(equations(katsura(5)))
        G = Systems.SPSystem(equations(cyclic(6)))
        H = Homotopies.HomotopyWithCache(Homotopies.StraightLineHomotopy(F, G), x, t)
        @test H isa Homotopies.AbstractHomotopy
        @test size(H) == (6, 6)

        m, n = size(H)
        u = zeros(Complex{Float64}, m)
        U = zeros(Complex{Float64}, m, n)

        Homotopies.evaluate!(u, H, x, t)
        @test Homotopies.evaluate(H, x, t) ≈ u

        Homotopies.dt!(u, H, x, t)
        @test Homotopies.dt(H, x, t) ≈ u

        Homotopies.jacobian!(U, H, x, t)
        @test Homotopies.jacobian(H, x, t) ≈ U

        Homotopies.evaluate_and_jacobian!(u, U, H, x, t)
        @test all(Homotopies.evaluate_and_jacobian(H, x, t) .≈ (u, U))
        @test all((u, U) .≈ (Homotopies.evaluate(H, x, t), Homotopies.jacobian(H, x, t)))

        Homotopies.jacobian_and_dt!(U, u, H, x, t)
        @test all(Homotopies.jacobian_and_dt(H, x, t) .≈ (U, u))
        @test all((U, u) .≈ (Homotopies.jacobian(H, x, t), Homotopies.dt(H, x, t)))
    end

    @testset "ParameterHomotopy" begin
        F = equations(katsura(5))
        vars = variables(F)
        H = Homotopies.ParameterHomotopy(F, vars[1:4], vars[5:6], rand(2), rand(2))
        @test H isa Homotopies.AbstractHomotopy
        @test size(H) == (6, 4)
        # @test Homotopies.gamma(H) isa ComplexF64
        # @test Homotopies.γ(H) isa ComplexF64

        InterfaceTest.homotopy(H)
    end

    @testset "FixedPointHomotopy" begin
        F = Systems.SPSystem(equations(katsura(5)))
        H = Homotopies.FixedPointHomotopy(F, rand(ComplexF64, 6))
        @test H isa Homotopies.AbstractHomotopy
        @test size(H) == (6, 6)

        InterfaceTest.homotopy(H)
    end

    @testset "PatchedHomotopy" begin
        F = Systems.SPSystem(equations(katsura(5)))
        G = Systems.SPSystem(equations(cyclic(6)))
        x = ProjectiveVectors.PVector(rand(Complex{Float64}, 6), 1)
        H = Homotopies.PatchedHomotopy(Homotopies.StraightLineHomotopy(F, G),
            AffinePatches.OrthogonalPatch(),
            x)
        @test H isa Homotopies.AbstractHomotopy
        @test size(H) == (7, 6)

        InterfaceTest.homotopy(H, x)
    end
end
