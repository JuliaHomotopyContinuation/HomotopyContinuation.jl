@testset "Composition" begin
    @testset "total degree" begin
        @polyvar a b c x y z u v
        e = [u + 1, v - 2]
        f = [a * b - 2, a * c - 1]
        g = [x + y, y + 3, x + 2]
        res = solve(e ∘ f ∘ g)
        @test nnonsingular(res) == 2
        @test nnonsingular(solve(e ∘ f ∘ g, system_scaling = nothing)) == 2
        @test nnonsingular(solve(e ∘ f ∘ g, system_scaling = :equations_and_variables)) == 2

        res = solve(e ∘ f ∘ g, system = SPSystem)
        @test nnonsingular(res) == 2
    end

    @testset "polyhedral" begin
        @polyvar x y u v
        result = solve(
            [2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y] ∘ [u, v];
            seed = 32241,
            start_system = :polyhedral,
        )
        @test nsolutions(result) == 6
    end

    @testset "parameter homotopy" begin
        @polyvar p q a b c x y z u v

        f = [a * b - 2, a * c - 1]
        g = [x + y, y + 3, x + 2]
        res = solve(f ∘ g; system = FPSystem, threading = false)
        @test nnonsingular(res) == 2
        # parameters at the end
        f2 = [a * b - q, a * c - p]
        g = [x + y, y + 3, x + 2]
        r = solve(
            f2 ∘ g,
            solutions(res);
            parameters = [p, q],
            p₁ = [1, 2],
            p₀ = [2, 3],
            threading = false,
        )
        @test nnonsingular(r) == 2

        # parameters at the beginning
        f = [a * b - 2, a * c - 1]
        g2 = [x + y, y + u, x + v]
        r = solve(
            f ∘ g2,
            solutions(res);
            parameters = [u, v],
            p₁ = [3, 2],
            p₀ = [-2, 3],
            threading = false,
        )
        @test nnonsingular(r) == 2

        # parameter in the middle
        e = [u + 1, v - 2]
        res2 = solve(e ∘ f ∘ g; system = SPSystem)
        f2 = [a * b - q, a * c - p]
        r = solve(
            e ∘ f2 ∘ g,
            solutions(res2);
            parameters = [p, q],
            p₁ = [1, 2],
            p₀ = [2, 3],
            threading = false,
        )
        @test nnonsingular(r) == 2
    end
end
