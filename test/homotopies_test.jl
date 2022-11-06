function test_homotopy_evaluate(homotopy, symbolic_homotopy)

    u = zeros(ComplexF64, length(homotopy))
    x = randn(ComplexF64, nvariables(homotopy))
    t = randn(ComplexF64)

    evaluate!(u, homotopy, x, t)
    @test u ≈ symbolic_homotopy(x, t) rtol = 1e-12
end

function test_homotopy_jacobian(homotopy, symbolic_homotopy)

    u = zeros(ComplexF64, length(homotopy))
    U = zeros(ComplexF64, length(homotopy), nvariables(homotopy))
    x = randn(ComplexF64, nvariables(homotopy))
    t = randn(ComplexF64)

    evaluate_and_jacobian!(u, U, homotopy, x, t)
    H_x = differentiate(symbolic_homotopy.expressions, symbolic_homotopy.variables)

    @test u ≈ symbolic_homotopy(x, t) rtol = 1e-12
    @test U ≈ evaluate(H_x, [symbolic_homotopy.variables; symbolic_homotopy.t] => [x; t]) rtol =
        1e-12
end



function test_homotopy_taylor(homotopy, symbolic_homotopy)
    for K = 1:3
        h = homotopy
        H = symbolic_homotopy

        u = zeros(ComplexF64, length(h))
        v = zeros(ComplexF64, length(h))
        w = zeros(ComplexF64, length(h))
        x = TaylorVector{K}(randn(ComplexF64, K, nvariables(h)))
        t = randn(ComplexF64)

        taylor!(u, Val(K), H, x, t)
        @var λ
        tx = [sum(xi .* λ .^ (0:length(xi)-1)) for xi in eachcol(x.data)]
        for k = 0:K
            true_value =
                (differentiate(Expression.(h(tx, t + λ)), λ, k)).(λ => 0) / factorial(k)
            @test u ≈ true_value rtol = 1e-12
            if k > 0
                v .= 0
                taylor!(v, Val(k), H, x, t)
                @test v ≈ true_value rtol = 1e-12
            end
        end
    end
end



@testset "Homotopies" begin
    @testset "ParameterHomotopy" begin
        using HomotopyContinuation
        @var x y a b c

        f = [(2 * x^2 + b^2 * y^3 + 2 * a * x * y)^3, (a + c)^4 * x + y^2]
        F = System(f, [x, y], [a, b, c])

        p = [5.2, -1.3, 9.3]
        q = [2.6, 3.3, 2.3]

        @var tvar
        f([x, y] => [x, y], [a, b, c] => tvar * p + (1 - tvar) * q)
        h = Homotopy(
            f([x, y] => [x, y], [a, b, c] => tvar * p + (1 - tvar) * q),
            [x, y],
            tvar,
        )

        H = HC.ParameterHomotopy(F, p, q)

        v = [0.192, 2.21]
        t = 0.232

        u = zeros(ComplexF64, 2)
        U = zeros(ComplexF64, 2, 2)

        HC.evaluate!(u, H, v, t)
        @test u ≈ f([x, y] => [x, y], [a, b, c] => t * p + (1 - t) * q)

        HC.taylor!(u, Val(1), H, TaylorVector{1}(Matrix(v')), t)
        @test u ≈ let
            @var s sp[1:3] sq[1:3]
            pt = s .* sp .+ (1 .- s) .* sq
            differentiate(subs(f, [a, b, c] => pt), s)(
                [x, y] => v,
                sp => p,
                sq => q,
                s => t,
            )
        end

        HC.evaluate_and_jacobian!(u, U, H, v, t)
        @test U ≈ differentiate(f, [x, y])([x, y] => v, [a, b, c] => t * p + (1 - t) * q)
    end

    @testset "AffineSubspaceHomotopy" begin
        @var x[1:4]
        f1 = rand_poly(x, 2)
        f2 = rand_poly(x, 2)
        F = System([f1, f2], x)
        A = rand_subspace(4; codim = 2)
        B = rand_subspace(4, codim = 2)
        H = IntrinsicSubspaceHomotopy(F, A, B)
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
    end
end
