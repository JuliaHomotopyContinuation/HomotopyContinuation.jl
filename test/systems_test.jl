@testset "Systems" begin

    @testset "SPSystem" begin
        fs = equations(katsura(5))

        F = SPSystem(fs)

        @test size(F) == (6, 6)
        @test length(F) == 6
        x = rand(Complex{Float64}, 6)
        system_cache = cache(F, x)
        @test system_cache isa HomotopyContinuation.SystemNullCache
        u = zeros(Complex{Float64}, 6)
        evaluate!(u, F, x, system_cache)
        @test evaluate(F, x, system_cache) ≈ u
        @test evaluate(F, x) ≈ u

        U = zeros(Complex{Float64}, 6, 6)
        jacobian!(U, F, x, system_cache)
        @test jacobian(F, x, system_cache) ≈ U
        @test jacobian(F, x) ≈ U

        evaluate_and_jacobian!(u, U, F, x, system_cache)
        @test all(evaluate_and_jacobian(F, x, system_cache) .≈ (u, U))
    end

    @testset "FixedHomotopy" begin
        f = SPSystem(equations(katsura(5)))
        g = SPSystem(equations(cyclic(6)))
        H = StraightLineHomotopy(f, g)

        F = FixedHomotopy(H, 0.23)
        x = rand(Complex{Float64}, 6)
        u = zeros(Complex{Float64}, 6)
        U = zeros(Complex{Float64}, 6, 6)

        system_cache = cache(F, x)
        @test system_cache isa HomotopyContinuation.FixedHomotopyCache
        u = zeros(Complex{Float64}, 6)
        evaluate!(u, F, x, system_cache)
        @test evaluate(F, x, system_cache) ≈ u

        U = zeros(Complex{Float64}, 6, 6)
        jacobian!(U, F, x, system_cache)
        @test jacobian(F, x, system_cache) ≈ U

        evaluate_and_jacobian!(u, U, F, x, system_cache)
        @test all(evaluate_and_jacobian(F, x, system_cache) .≈ (u, U))
    end

    @testset "FPSystem" begin
        F = FPSystem(equations(katsura(5)))
        x = rand(Complex{Float64}, 6)
        u = zeros(Complex{Float64}, 6)
        U = zeros(Complex{Float64}, 6, 6)

        system_cache = cache(F, x)
        @test system_cache isa HomotopyContinuation.FPSystemCache
        u = zeros(Complex{Float64}, 6)
        evaluate!(u, F, x, system_cache)
        @test evaluate(F, x, system_cache) ≈ u

        U = zeros(Complex{Float64}, 6, 6)
        jacobian!(U, F, x, system_cache)
        @test jacobian(F, x, system_cache) ≈ U

        evaluate_and_jacobian!(u, U, F, x, system_cache)
        @test all(evaluate_and_jacobian(F, x, system_cache) .≈ (u, U))
    end

    @testset "TotalDegreeSystem" begin
        @polyvar x y z

        F = SPSystem([x^4-z^4, y^3-z^3])
        w = rand(ComplexF64, 3)

        G = TotalDegreeSystem([x^4-z^4, y^3-z^3])
        @test evaluate(F, w) ≈ evaluate(G, w)
        @test jacobian(F, w) ≈ jacobian(G, w)
        u, U = evaluate_and_jacobian(F, w)
        @test u ≈ evaluate(F, w)
        @test U ≈ jacobian(G, w)

        # @polyvar x y z t
        # F = SPSystem([x^4-z^4, y^3-z^3, t-z])
        # G = TotalDegreeSystem([x^4-z^4, y^3-z^3, t-z], [x, y, z, t], z)
        # w = rand(ComplexF64, 4)
        # @test evaluate(F, w) ≈ evaluate(G, w)
        # @test jacobian(F, w) ≈ jacobian(G, w)
        # u, U = evaluate_and_jacobian(G, w)
        # @test u ≈ evaluate(F, w)
        # @test U ≈ jacobian(F, w)
    end

    @testset "FixedParameterSystem" begin
        @polyvar x y a b
        f = SPSystem([x^2+y-a, b*x+y*a^2]; parameters=[a, b])
        F = FixedParameterSystem(f, [2.0, 3.0])

        @test size(F) == (2, 2)
        @test length(F) == 2
        x = rand(Complex{Float64}, 2)
        system_cache = cache(F, x)
        @test system_cache isa HomotopyContinuation.FixedParameterSystemCache
        u = zeros(Complex{Float64}, 2)
        evaluate!(u, F, x, system_cache)
        @test evaluate(F, x, system_cache) ≈ u
        @test evaluate(F, x) ≈ u

        U = zeros(Complex{Float64}, 2, 2)
        jacobian!(U, F, x, system_cache)
        @test jacobian(F, x, system_cache) ≈ U
        @test jacobian(F, x) ≈ U

        evaluate_and_jacobian!(u, U, F, x, system_cache)
        @test all(evaluate_and_jacobian(F, x, system_cache) .≈ (u, U))

        set_parameters!(F, (4.0, 2.3))
        @test F.p == [4.0, 2.3]
    end
end
