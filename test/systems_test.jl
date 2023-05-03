function test_system_evaluate(system, symbolic_system)
    m, n = size(system)
    k = nparameters(system)
    u = zeros(ComplexF64, m)
    ū = zeros(ComplexDF64, m)
    x = randn(ComplexF64, n)
    x̄ = ComplexDF64.(x)
    p = randn(ComplexF64, k)

    evaluate!(u, system, x, p)
    @test u ≈ symbolic_system(x, p) rtol = 1e-12

    evaluate!(ū, system, x̄, p)
    @test ū ≈ symbolic_system(x, p) rtol = 1e-12
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
        true_taylor_value = [
            transpose(
                (differentiate(Expression.(symbolic_system(tx, tp)), λ, j)).(λ => 0) /
                factorial(j),
            ) for j = 0:K
        ]

        w = TaylorVector{K + 1}(vcat(true_taylor_value...))
        v .= 0

        taylor!(v, Val(k), system, x, p)
        @test v ≈ last(vectors(w)) rtol = 1e-12

        u = TaylorVector{K + 1}(zeros(ComplexF64, K + 1, m))
        taylor!(u, Val(k), system, x, p)

        for (ui, wi) in zip(vectors(u), vectors(w))
            @test ui ≈ wi rtol = 1e-12
        end
    end
end

function test_system(system, symbolic_system)
    test_system_evaluate(system, symbolic_system)
    test_system_jacobian(system, symbolic_system)
    test_system_taylor(system, symbolic_system)
end


@testset "System" begin
    @testset "InterpretedSystem" begin
        @var x y a b
        f = System([y^2 + 2x + 3, (x - 1)^3])

        test_system(InterpretedSystem(f), f)
    end

    @testset "CompiledSystem" begin
        @var x y a b
        f = System([y^2 + 2x + 3, (x - 1)^3])

        test_system(CompiledSystem(f), f)
    end


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

    @testset "LinearSystem" begin
        @var x[1:8]
        A = randn(ComplexF64, 4, 8)
        b = randn(ComplexF64, 4)
        f = System(A * x - b)
        L = LinearSystem(A, b; variables = x)

        test_system(L, f)
    end


    @testset "StackedSystem - 1" begin
        @var x y a b
        f = [y^2 + 2x + 3, x - 1]
        F = System(f)
        g = [(x^4 + y * a)^3 - 3y, x - 1 * b^2]
        G = System(g, parameters = [a, b])

        test_system(StackedSystem(F, G), System([f; g]; parameters = [a, b]))
    end

    @testset "StackedSystem - 2" begin
        @var x y
        f = [y^2 + 2x + 3, x - 1]
        F = System(f)
        A = randn(ComplexF64, 3, 2)
        b = randn(ComplexF64, 3)
        g = A * [x, y] - b
        G = LinearSystem(A, b; variables = [x, y])

        test_system(StackedSystem(F, G), System([f; g]))
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
