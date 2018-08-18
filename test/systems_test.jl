@testset "Systems" begin

    @testset "SPSystem" begin
        fs = equations(katsura(5))

        F = Systems.SPSystem(fs)

        @test size(F) == (6, 6)
        @test length(F) == 6
        x = rand(Complex{Float64}, 6)
        cache = Systems.cache(F, x)
        @test cache isa Systems.NullCache
        u = zeros(Complex{Float64}, 6)
        Systems.evaluate!(u, F, x, cache)
        @test Systems.evaluate(F, x, cache) ≈ u
        @test Systems.evaluate(F, x) ≈ u

        U = zeros(Complex{Float64}, 6, 6)
        Systems.jacobian!(U, F, x, cache)
        @test Systems.jacobian(F, x, cache) ≈ U
        @test Systems.jacobian(F, x) ≈ U

        Systems.evaluate_and_jacobian!(u, U, F, x, cache)
        @test all(Systems.evaluate_and_jacobian(F, x, cache) .≈ (u, U))
    end

    @testset "FixedHomotopy" begin
        f = Systems.SPSystem(equations(katsura(5)))
        g = Systems.SPSystem(equations(cyclic(6)))
        H = Homotopies.StraightLineHomotopy(f, g)

        F = Systems.FixedHomotopy(H, 0.23)
        x = rand(Complex{Float64}, 6)
        u = zeros(Complex{Float64}, 6)
        U = zeros(Complex{Float64}, 6, 6)

        cache = Systems.cache(F, x)
        @test cache isa Systems.FixedHomotopyCache
        u = zeros(Complex{Float64}, 6)
        Systems.evaluate!(u, F, x, cache)
        @test Systems.evaluate(F, x, cache) ≈ u

        U = zeros(Complex{Float64}, 6, 6)
        Systems.jacobian!(U, F, x, cache)
        @test Systems.jacobian(F, x, cache) ≈ U

        Systems.evaluate_and_jacobian!(u, U, F, x, cache)
        @test all(Systems.evaluate_and_jacobian(F, x, cache) .≈ (u, U))
    end

    @testset "FPSystem" begin
        F = Systems.FPSystem(equations(katsura(5)))
        x = rand(Complex{Float64}, 6)
        u = zeros(Complex{Float64}, 6)
        U = zeros(Complex{Float64}, 6, 6)

        cache = Systems.cache(F, x)
        @test cache isa Systems.FPSystemCache
        u = zeros(Complex{Float64}, 6)
        Systems.evaluate!(u, F, x, cache)
        @test Systems.evaluate(F, x, cache) ≈ u

        U = zeros(Complex{Float64}, 6, 6)
        Systems.jacobian!(U, F, x, cache)
        @test Systems.jacobian(F, x, cache) ≈ U

        Systems.evaluate_and_jacobian!(u, U, F, x, cache)
        @test all(Systems.evaluate_and_jacobian(F, x, cache) .≈ (u, U))
    end

    @testset "TotalDegreeSystem" begin
        @polyvar x y z

        F = Systems.SPSystem([x^4-z^4, y^3-z^3])
        w = rand(ComplexF64, 3)

        G = Systems.TotalDegreeSystem([x^4-z^4, y^3-z^3], [x, y, z], z)
        @test Systems.evaluate(F, w) ≈ Systems.evaluate(G, w)
        @test Systems.jacobian(F, w) ≈ Systems.jacobian(G, w)
        u, U = Systems.evaluate_and_jacobian(F, w)
        @test u ≈ Systems.evaluate(F, w)
        @test U ≈ Systems.jacobian(G, w)

        @polyvar x y z t
        F = Systems.SPSystem([x^4-z^4, y^3-z^3, t-z])
        G = Systems.TotalDegreeSystem([x^4-z^4, y^3-z^3, t-z], [x, y, z, t], z)
        w = rand(ComplexF64, 4)
        @test Systems.evaluate(F, w) ≈ Systems.evaluate(G, w)
        @test Systems.jacobian(F, w) ≈ Systems.jacobian(G, w)
        u, U = Systems.evaluate_and_jacobian(G, w)
        @test u ≈ Systems.evaluate(F, w)
        @test U ≈ Systems.jacobian(F, w)
    end
end
