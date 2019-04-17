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
        w = ProjectiveVectors.PVector(rand(ComplexF64, 3))

        G = TotalDegreeSystem([x^4-z^4, y^3-z^3])
        @test size(G) == (2, 3)
        @test evaluate(F, w) ≈ evaluate(G, w)
        @test jacobian(F, w) ≈ jacobian(G, w)
        u, U = evaluate_and_jacobian(F, w)
        @test u ≈ evaluate(F, w)
        @test U ≈ jacobian(G, w)

        @polyvar x y z

        F_affine = SPSystem([x^4-1, y^3-1])
        w = rand(ComplexF64, 3)

        G_affine = TotalDegreeSystem([x^4-1, y^3-1]; affine=true)
        @test size(G_affine) == (2, 2)
        @test evaluate(F_affine, w) ≈ evaluate(G_affine, w)
        @test jacobian(F_affine, w) ≈ jacobian(G_affine, w)

        u, U = evaluate_and_jacobian(F_affine, w)
        @test u ≈ evaluate(F_affine, w)
        @test U ≈ jacobian(G_affine, w)
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

    @testset "SquaredUpSystem" begin
        # affine
        @polyvar x y z

        F = SPSystem([x^4-1, y^3-1, z^2+x, x+y+z-1, x^2+z^2-3])
        w = rand(ComplexF64, 3)
        A = randn(ComplexF64, 3, 2)
        S = SquaredUpSystem(F, A)
        system_cache = cache(S, w)
        @test system_cache isa AbstractSystemCache

        @test size(S) == (3, 3)
        @test evaluate(S, w, system_cache) ≈ [LinearAlgebra.I A] * evaluate(F, w) atol=1e-12
        @test jacobian(S, w, system_cache) ≈ [LinearAlgebra.I A] * jacobian(F, w) atol=1e-12
        u, U = evaluate_and_jacobian(S, w, system_cache)
        @test u ≈ evaluate(S, w)
        @test U ≈ jacobian(S, w)


        # projective
        @polyvar x y z w

        A = randn(ComplexF64, 3, 2)
        F = [x^4-1, y^3-1, z^2+x, x^2+z^2-3, x+y+z-1]
        G = SPSystem(homogenize([LinearAlgebra.I A] * F, w))
        F = SPSystem(homogenize(F, w))
        v = ProjectiveVectors.PVector(rand(ComplexF64, 4))

        S = SquaredUpSystem(F, A, [4, 3, 2, 2, 1])
        system_cache = cache(S, v)
        @test system_cache isa AbstractSystemCache

        @test size(S) == (3, 4)
        @test evaluate(S, v, system_cache) ≈ evaluate(G, v) atol=1e-12
        @test jacobian(S, v, system_cache) ≈ jacobian(G, v) atol=1e-12
        u, U = evaluate_and_jacobian(S, v, system_cache)
        @test u ≈ evaluate(S, v)
        @test U ≈ jacobian(S, v)

        @test HC.check_homogeneous_degrees(S) == [4, 3, 2]
    end
end
