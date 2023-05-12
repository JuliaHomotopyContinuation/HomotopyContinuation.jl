function test_homotopy_evaluate(homotopy, symbolic_homotopy)
    m, n = size(homotopy)
    u = zeros(ComplexF64, m)
    x = randn(ComplexF64, n)
    t = randn(ComplexF64)

    evaluate!(u, homotopy, x, t)
    @test u ≈ symbolic_homotopy(x, t) rtol = 1e-12
end

function test_homotopy_jacobian(homotopy, symbolic_homotopy)
    m, n = size(homotopy)
    u = zeros(ComplexF64, m)
    U = zeros(ComplexF64, m, n)
    x = randn(ComplexF64, n)
    t = randn(ComplexF64)

    evaluate_and_jacobian!(u, U, homotopy, x, t)
    H_x = differentiate(symbolic_homotopy.expressions, symbolic_homotopy.variables)

    @test u ≈ symbolic_homotopy(x, t) rtol = 1e-12
    @test U ≈ evaluate(H_x, [symbolic_homotopy.variables; symbolic_homotopy.t] => [x; t]) rtol =
        1e-12
end

function test_homotopy_taylor(homotopy, symbolic_homotopy)
    m, n = size(homotopy)
    v = zeros(ComplexF64, m)
    t = randn(ComplexF64)
    ṫ = rand()
    X = TaylorVector{4}(randn(ComplexF64, 4, n))

    for incr in [false, true]
        for K = 1:4
            @testset "Taylor, K = $K, incr = $incr" begin
                k = K
                x = TaylorVector{K}(X)
                @var λ
                tx = [sum(xi .* λ .^ (0:length(xi)-1)) for xi in eachcol(x.data[1:K, :])]
                true_taylor_value =
                    transpose.([
                        (differentiate(
                            Expression.(symbolic_homotopy(tx, t + ṫ * λ)),
                            λ,
                            j,
                        )).(λ => 0) / factorial(j) for j = 0:K
                    ])
                v .= 0

                w = TaylorVector{K + 1}(vcat(true_taylor_value...))
                v .= 0

                taylor!(v, Val(k), homotopy, x, (t, ṫ))
                @test v ≈ last(vectors(w)) rtol = 1e-12

                u = TaylorVector{K + 1}(zeros(ComplexF64, K + 1, m))
                taylor!(u, Val(k), homotopy, x, (t, ṫ))

                j = 0
                for (ui, wi) in zip(vectors(u), vectors(w))
                    j += 1
                    @test ui ≈ wi rtol = 1e-12
                end

            end
        end
    end

end

function test_homotopy(homotopy, symbolic_homotopy)

    test_homotopy_evaluate(homotopy, symbolic_homotopy)
    test_homotopy_jacobian(homotopy, symbolic_homotopy)
    test_homotopy_taylor(homotopy, symbolic_homotopy)

end


@testset "Homotopies" begin
    @testset "ParameterHomotopy" begin
        @var x y a b c

        f = [(2 * x^2 + b^2 * y^3 + 2 * a * x * y)^3, (a + c)^4 * x + y^2]
        F = System(f, [x, y], [a, b, c])

        p = [5.2, -1.3, 9.3]
        q = [2.6, 3.3, 2.3]

        @var tvar
        h = Homotopy(
            f([x, y] => [x, y], [a, b, c] => tvar * p + (1 - tvar) * q),
            [x, y],
            tvar,
        )

        H = ParameterHomotopy(F, p, q; compile = false)

        test_homotopy(H, h)
    end
    @testset "StraightLineHomotopy" begin
        @var x y
        a, b, c, = rand(3)
        f = [(2 * x^2 + b^2 * y^3 + 2 * a * x * y)^3, (a + c)^4 * x + y^2]
        F = System(f, [x, y])
        g = [(2 * y^2 + b^2 * x^3 + 2 * a * x * y)^2, (a - c)^3 * y + x^2]
        G = System(g, [x, y])
        H = StraightLineHomotopy(F, G; compile = false)

        @var t
        h = Homotopy(t .* F([x, y]) + (1 - t) .* G([x, y]), [x, y], t)

        test_homotopy(H, h)
    end

    @testset "ExtrinsicLinearSubspaceHomotopy" begin
        L₁ = rand_subspace(4; codim = 2)
        L₂ = rand_subspace(4; codim = 2)
        @var x[1:4]
        f1 = rand_poly(x, 2)
        f2 = rand_poly(x, 2)
        F = System([f1, f2], x)
        H = ExtrinsicLinearSubspaceHomotopy(F, L₁, L₂; compile = false)
        @var t
        h = Homotopy([f1; f2; t .* L₁(x) + (1 - t) .* L₂(x)], x, t)
        test_homotopy(H, h)
    end

    @testset "SystemHomotopyComposition" begin
        L₁ = rand_subspace(4; codim = 3)
        L₂ = rand_subspace(4; codim = 3)
        @var x[1:4]
        f1 = rand_poly(x, 2)
        f2 = rand_poly(x, 2)
        F = System([f1, f2], x)
        H = ExtrinsicLinearSubspaceHomotopy(F, L₁, L₂; compile = false)
        @var y[1:5]
        g1 = rand_poly(y, 3)
        g2 = rand_poly(y, 2)
        G = System([g1, g2], y)
        K = SystemHomotopyComposition(fixed(G), H)

        @var t
        h = [f1; f2; t .* L₁(x) + (1 - t) .* L₂(x)]
        k = Homotopy(G(h), x, t)
        test_homotopy(K, k)
    end


    @testset "SystemHomotopyComposition" begin
        L₁ = rand_subspace(5; codim = 2)
        L₂ = rand_subspace(5; codim = 2)
        @var x[1:5]
        f1 = rand_poly(x, 3)
        f2 = rand_poly(x, 2)
        F = System([f1, f2], x)
        H = IntrinsicLinearSubspaceHomotopy(F, L₁, L₂)

        @var y[1:3]
        @var t

        h = Homotopy(F(t .* intrinsic(L₁)(y) + (1 - t) .* intrinsic(L₂)(y)), y, t)
        test_homotopy(H, h)
    end


    @testset "LinearSubspaceGeodesicHomotopy" begin
        @var x[1:4]
        f1 = rand_poly(x, 2)
        f2 = rand_poly(x, 2)
        F = System([f1, f2], x)
        A = rand_subspace(4; codim = 2)
        B = rand_subspace(4, codim = 2)
        H = LinearSubspaceGeodesicHomotopy(fixed(F; compile = false), A, B)
        @unpack Q, Q_cos, Θ = H.geodesic

        @var v[1:3] t
        s = sin.(t .* Θ)
        c = cos.(t .* Θ)
        γ = Q_cos .* c' .+ Q .* s'

        h = Homotopy([F((γ*v)[1:end-1]); (γ*v)[end] - 1], v, t)

        x0 = randn(ComplexF64, 4)
        t0 = ComplexF64(rand())
        v0 = zeros(ComplexF64, 3)
        HomotopyContinuation.set_solution!(v0, H, x0, t0)

        out0 = zeros(ComplexF64, 3)
        evaluate!(out0, H, v0, t0)

        out1 = zeros(ComplexF64, 3)
        evaluate!(out1, InterpretedHomotopy(h), v0, t0)

        @test out0 ≈ out1 rtol = 1e-12


        test_homotopy(H, h)
    end


    @testset "MixedHomotopy" begin
        @var x y t
        a, b, c, = rand(3)
        f = [(2 * x^2 + b^2 * y^3 + 2 * a * x * y)^3, (a + c)^4 * x + y^2]
        F = System(f, [x, y])
        g = [(2 * y^2 + b^2 * x^3 + 2 * a * x * y)^2, (a - c)^3 * y + x^2]
        G = System(g, [x, y])
        h = Homotopy(t .* F([x, y]) + (1 - t) .* G([x, y]), [x, y], t)
        H = MixedHomotopy(h)
        test_homotopy(H, h)

        @var x t
        H = Homotopy([(1 - 2 * t) + (-2 + 2 * t) * x + x^2], [x], t)
        test_homotopy(MixedHomotopy(H), H)
    end
end
