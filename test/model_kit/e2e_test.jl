include("../test_systems.jl")

@testset "interpreter: $name" for (name, system) in TEST_SYSTEM_COLLECTION
    @testset "interpreter symbolic" begin
        # skip bacillus since it has fp constants
        if name !== "bacillus"
            F = ModelKit.interpreter(Expression, system)

            u = Vector{Any}(randn(ComplexF64, length(system)))
            U = Matrix{Any}(randn(ComplexF64, length(system), nvariables(system)))
            x = variables(system)
            p = parameters(system)

            ModelKit.execute!(u, F, x, p)
            @test sum(expand.(u - system.expressions)) == 0
        end
    end

    @testset "$mode [ComplexF64]" for mode in [InterpretedSystem, CompiledSystem]
        F = mode(system)
        jf = System(
            vec(jacobian(system));
            variables = variables(system),
            parameters = parameters(system),
        )

        u = randn(ComplexF64, length(system))
        U = randn(ComplexF64, length(system), nvariables(system))
        x = randn(ComplexF64, nvariables(system))
        p = randn(ComplexF64, nparameters(system))

        J(x, p) = begin
            M = jf(x, p)
            map(M) do m
                if m isa Expression
                    ModelKit.to_number(expand(m))
                else
                    m
                end
            end
        end
        evaluate!(u, F, x, p)
        @test u ≈ system(x, p) rtol = 1e-12

        u .= 0
        evaluate_and_jacobian!(u, U, F, x, p)
        @test u ≈ system(x, p) rtol = 1e-12
        @test vec(U) ≈ J(x, p) rtol = 1e-12

        U .= 0
        jacobian!(U, F, x, p)
        @test vec(U) ≈ J(x, p) rtol = 1e-12
    end

    @testset "$mode evaluate [ComplexDF64]" for mode in [InterpretedSystem, CompiledSystem]
        F = mode(system)

        u = zeros(ComplexDF64, length(system))
        x = ComplexDF64.(randn(ComplexF64, nvariables(system)))
        p = randn(ComplexF64, nparameters(system))

        evaluate!(u, F, x, p)
        @test u ≈ system(x, p) rtol = 1e-12
    end

    @testset "$mode taylor order $K [ComplexF64]" for mode in
                                                      [InterpretedSystem, CompiledSystem],
        K = 1:3

        F = mode(system)

        u = TaylorVector{K + 1}(ComplexF64, length(system))
        v = zeros(ComplexF64, length(system))
        x = TaylorVector{K}(randn(ComplexF64, K, nvariables(system)))
        p = randn(ComplexF64, nparameters(system))

        taylor!(u, Val(K), F, x, p)
        @var λ
        tx = [sum(xi .* λ .^ (0:length(xi)-1)) for xi in eachcol(x.data)]
        for k = 0:K
            true_value =
                (differentiate(Expression.(system(tx, p)), λ, k)).(λ => 0) / factorial(k)
            @test vectors(u)[k+1] ≈ true_value rtol = 1e-12
            if k > 0
                v .= 0
                taylor!(v, Val(k), F, x, p)
                @test v ≈ true_value rtol = 1e-12
            end
        end

    end

    @testset "interpreted system [acb]" begin
        F = InterpretedSystem(system)
        jf = System(
            vec(jacobian(system));
            variables = variables(system),
            parameters = parameters(system),
        )

        u = Arblib.AcbRefVector(length(system))
        U = Arblib.AcbRefMatrix(length(system), nvariables(system))
        x = randn(ComplexF64, nvariables(system))
        p = randn(ComplexF64, nparameters(system))
        acb_x = Arblib.AcbRefVector(x)
        acb_p = isempty(p) ? Arblib.AcbRefVector(0) : Arblib.AcbRefVector(p)

        J(x, p) = begin
            M = jf(x, p)
            map(M) do m
                if m isa Expression
                    ModelKit.to_number(expand(m))
                else
                    m
                end
            end
        end
        evaluate!(u, F, acb_x, acb_p)
        @test ComplexF64.(u) ≈ system(x, p) rtol = 1e-12

        ModelKit.zero!(u)
        evaluate_and_jacobian!(u, U, F, acb_x, acb_p)
        @test ComplexF64.(u) ≈ system(x, p) rtol = 1e-12
        @test ComplexF64.(vec(U)) ≈ J(x, p) rtol = 1e-12

        ModelKit.zero!(U)
        jacobian!(U, F, acb_x, acb_p)
        @test ComplexF64.(vec(U)) ≈ J(x, p) rtol = 1e-12
    end

    @testset "$mode [ComplexF64]" for mode in [InterpretedHomotopy, CompiledHomotopy]
        @var __t__
        h = Homotopy(
            (1 - __t__) .* system.expressions + __t__^2 .* system.expressions,
            variables(system),
            __t__,
            parameters(system),
        )
        H = mode(h)
        jf = System(
            vec(jacobian(system));
            variables = variables(system),
            parameters = parameters(system),
        )

        u = randn(ComplexF64, length(system))
        U = randn(ComplexF64, length(system), nvariables(system))
        x = randn(ComplexF64, nvariables(system))
        t = randn(ComplexF64)
        p = randn(ComplexF64, nparameters(system))

        J(x, p) = begin
            M = jf(x, p)
            j = map(M) do m
                if m isa Expression
                    ModelKit.to_number(expand(m))
                else
                    m
                end
            end
            (1 - t) .* j + t^2 .* j
        end
        evaluate!(u, H, x, t, p)
        @test u ≈ (1 - t) .* system(x, p) + t^2 .* system(x, p) rtol = 1e-12

        u .= 0
        evaluate_and_jacobian!(u, U, H, x, t, p)
        @test u ≈ (1 - t) .* system(x, p) + t^2 .* system(x, p) rtol = 1e-12
        @test vec(U) ≈ J(x, p) rtol = 1e-12

        U .= 0
        jacobian!(U, H, x, t, p)
        @test vec(U) ≈ J(x, p) rtol = 1e-12
    end

    @testset "$mode evaluate [ComplexDF64]" for mode in
                                                [InterpretedHomotopy, CompiledHomotopy]
        @var __t__
        h = Homotopy(
            (1 - __t__) .* system.expressions + __t__^2 .* system.expressions,
            variables(system),
            __t__,
            parameters(system),
        )
        H = mode(h)


        u = zeros(ComplexDF64, length(system))
        x = ComplexDF64.(randn(ComplexF64, nvariables(system)))
        t = randn(ComplexF64)
        p = randn(ComplexF64, nparameters(system))
        evaluate!(u, H, x, t, p)
        @test u ≈ (1 - t) .* system(x, p) + t^2 .* system(x, p) rtol = 1e-12
    end

    @testset "$mode taylor order $K [ComplexF64]" for mode in [
            InterpretedHomotopy,
            CompiledHomotopy,
        ],
        K = 1:3

        @var __t__
        h = Homotopy(
            (1 - __t__) .* system.expressions + __t__^2 .* system.expressions,
            variables(system),
            __t__,
            parameters(system),
        )

        F = mode(h)

        u = TaylorVector{K + 1}(ComplexF64, length(h))
        v = zeros(ComplexF64, length(h))
        x = TaylorVector{K}(randn(ComplexF64, K, nvariables(h)))
        t = randn(ComplexF64)
        p = randn(ComplexF64, nparameters(h))

        taylor!(u, Val(K), F, x, t, p)
        @var λ
        tx = [sum(xi .* λ .^ (0:length(xi)-1)) for xi in eachcol(x.data)]
        for k = 0:K
            true_value =
                (differentiate(Expression.(h(tx, t + λ, p)), λ, k)).(λ => 0) / factorial(k)
            @test vectors(u)[k+1] ≈ true_value rtol = 1e-12
            if k > 0
                v .= 0
                taylor!(v, Val(k), F, x, t, p)
                @test v ≈ true_value rtol = 1e-12
            end
        end

    end
end
