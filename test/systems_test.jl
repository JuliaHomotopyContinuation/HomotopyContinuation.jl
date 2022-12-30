function test_system_evaluate(system, symbolic_system)
    m, n = size(system)
    k = nparameters(system)
    u = zeros(ComplexF64, m)
    x = randn(ComplexF64, n)
    p = randn(ComplexF64, k)

    evaluate!(u, system, x, p)
    @test u ≈ symbolic_system(x, p) rtol = 1e-12
end

function test_system_jacobian(system, symbolic_system)
    m, n = size(system)
    k = nparameters(system)
    u = zeros(ComplexF64, m)
    U = zeros(ComplexF64, m, n)
    x = randn(ComplexF64, n)
    p = randn(ComplexF64, k)

    evaluate_and_jacobian!(u, U, system, x, p)
    H_x = differentiate(symbolic_system.expressions, symbolic_system.variables)

    @test u ≈ symbolic_system(x, p) rtol = 1e-12
    @test U ≈
          evaluate(H_x, [symbolic_system.variables; symbolic_system.parameters] => [x; p]) rtol =
        1e-12
end

function test_system_taylor(system, symbolic_system)
    m, n = size(system)
    r = nparameters(system)
    v = zeros(ComplexF64, m)

    for K = 1:4
        k = K
        x = TaylorVector{K}(randn(ComplexF64, K, n))
        p = TaylorVector{K}(randn(ComplexF64, K, r))

        @var λ
        tx = [sum(xi .* λ .^ (0:length(xi)-1)) for xi in eachcol(x.data)]
        tp = [sum(p_i .* λ .^ (0:length(p_i)-1)) for p_i in eachcol(p.data)]
        true_value =
            (differentiate(Expression.(symbolic_system(tx, tp)), λ, k)).(λ => 0) /
            factorial(k)
        v .= 0
        taylor!(v, Val(k), system, x, p)
        @test v ≈ true_value rtol = 1e-12

    end
end

function test_system(system, symbolic_system)
    test_system_evaluate(system, symbolic_system)
    test_system_jacobian(system, symbolic_system)
    test_system_taylor(system, symbolic_system)
end


@testset "System" begin
    @testset "CompositionSystem" begin
        @var x y a b
        f = System([y^2 + 2x + 3, x - 1])
        g = System([(x^4 + y * a)^3 - 3y, x - 1 * b^2], parameters = [a, b])
        comp = compose(g, f; compile = false)

        test_system(comp, System(comp))
    end
    @testset "CompositionSystem - Compose" begin
        @var x y a b

        f = System([y^2 + 2x + 3, x - 1])
        g = System([x + y * a, x - 1 * b], parameters = [a, b])
        @test parameters(g ∘ f) == [a, b]
        @test nparameters(g ∘ f) == 2

        @test parameters(g ∘ g) == [a, b]
        @test nparameters(g ∘ g) == 2

        @test parameters(f ∘ g) == [a, b]
        @test nparameters(f ∘ g) == 2

        @test parameters(f ∘ f) == []
        @test nparameters(f ∘ f) == 0

        @test parameters(f ∘ g ∘ f) == [a, b]
        @test nparameters(f ∘ g ∘ f) == 2
    end

    @testset "MixedSystem" begin
        @var x y a b
        g = System([(x^4 + y * a)^3 - 3y, x - 1 * b^2], parameters = [a, b])

        test_system(MixedSystem(g), g)
    end

    @testset "AffineChartSystem" begin
        @var x y a b
        g = System([x^2 + 3 * x * y + y^2 * b], parameters = [a, b])
        G = on_affine_chart(g)
        test_system(G, System(G))
    end

    @testset "FixedParameterSystem" begin
        @var x y a b
        g = System([x^2 + 3 * x * y + y^2 * b], parameters = [a, b])
        G = FixedParameterSystem(g, [2.3, 4.1])
        test_system(G, System(G))
    end

    @testset "RandomizedSystem" begin
        @var x y a b
        g = System(
            [
                x^2 + 3 * x * y + y^2 * b,
                (x + 3y - 2)^2 / 2,
                3 * a * b * (2x - y + y^3 - 4x^2),
            ],
            parameters = [a, b],
        )
        G = RandomizedSystem(g, 2)
        test_system(G, System(G))
    end

    @testset "SlicedSystem" begin
        @var x y a b
        g = System(
            [
                x^2 + 3 * x * y + y^2 * b,
                (x + 3y - 2)^2 / 2,
                3 * a * b * (2x - y + y^3 - 4x^2),
            ],
            parameters = [a, b],
        )
        L = rand_subspace(2; dim = 1)
        G = SlicedSystem(g, L)
        test_system(G, System(G))
    end
end
