@testset "ModelKit" begin
    @testset "System (no parameters)" begin
        ModelKit.@var x y z
        function F_test(x)
            let (x, y, z) = (x[1], x[2], x[3])
                [
                 (x^2 + y + z + 2)^2 - 3,
                 4 * x^2 * z^2 * y + 4z - 6x * y * z^7,
                 (-z) / (x + y),
                ]
            end
        end

        f = ModelKit.System(F_test([x, y, z]), [x, y, z])

        @test sprint(show, f) == """
        System
         variables: x, y, z

         (x ^ 2 + y + z + 2) ^ 2 - 3
         (4 * x ^ 2 * z ^ 2 * y + 4z) - 6 * x * y * z ^ 7
         -z / (x + y)"""

        F = ModelKit.compile(f)

        @test size(F) == (3, 3)
        @test length(F) == 3
        @test size(F, 1) == 3
        @test size(F, 2) == 3
        u = [0.0, 0.0, 0]
        U = zeros(3, 3)
        v = [0.4814, -0.433, -0.82709]

        fdm = FD.central_fdm(5, 1)
        û = F_test(v)
        @test evaluate(F, v) ≈ û
        @test F(v) ≈ û
        @test ModelKit.interpreted(F)(v) ≈ û
        @test evaluate(ModelKit.interpreted(F), v) ≈ û
        @test evaluate!(u, F, v) ≈ û
        @test u ≈ û
        Û = FD.jacobian(fdm, F_test, v)[1]
        @test jacobian(F, v) ≈ Û
        @test jacobian(ModelKit.interpreted(F), v) ≈ Û
        @test jacobian!(U, F, v) ≈ Û
        @test U ≈ Û rtol = 1e-8
        u .= 0
        U .= 0
        @test evaluate_and_jacobian!(u, U, F, v) === nothing
        @test u ≈ û
        @test U ≈ Û
    end

    @testset "System (parameters)" begin
        ModelKit.@var x y z a b c
        function F_test(x, p)
            let (x, y, z, b, c, a) = (x[1], x[2], x[3], p[1], p[2], p[3])
                [
                 (x^2 + y + z * a^2 + 2)^2 - 3,
                 4 * x^2 * z^2 * y + b + 4z - 6x * y * z^2,
                 a * c * (-z) / (x + y),
                ]
            end
        end

        f = ModelKit.System(
            [
             (x^2 + y + z * a^2 + 2)^2 - 3,
             4 * x^2 * z^2 * y + b + 4z - 6x * y * z^2,
             a * c * (-z) / (x + y),
            ],
            [x, y, z],
            [b, c, a],
        )
        @test sprint(show, f) == """
        System
         variables: x, y, z
         parameters: b, c, a

         (x ^ 2 + y + z * a ^ 2 + 2) ^ 2 - 3
         (4 * x ^ 2 * z ^ 2 * y + b + 4z) - 6 * x * y * z ^ 2
         (a * c * -z) / (x + y)"""

        F = ModelKit.compile(f)
        @test startswith(
            sprint(show, F),
            "Compiled{TSystem{3,3,3,#",
        )

        @test size(F) == (3, 3)
        @test length(F) == 3
        @test size(F, 1) == 3
        @test size(F, 2) == 3
        u = [0.0, 0.0, 0]
        U = zeros(3, 3)
        v = [0.4814, -0.433, -0.82709]
        p = [0.1312, 2.1232, 1.2411]

        fdm = FD.central_fdm(5, 1)
        û = F_test(v, p)
        @test evaluate(F, v, p) ≈ û
        @test F(v, p) ≈ û
        @test evaluate(ModelKit.interpreted(F), v, p) ≈ û
        @test ModelKit.interpreted(F)(v, p) ≈ û
        @test evaluate!(u, F, v, p) ≈ û
        @test u ≈ û
        Û = FD.jacobian(fdm, v -> F_test(v, p), v)[1]
        @test jacobian(F, v, p) ≈ Û
        @test jacobian(ModelKit.interpreted(F), v, p) ≈ Û
        @test jacobian!(U, F, v, p) ≈ Û
        @test U ≈ Û rtol = 1e-8
        u .= 0
        U .= 0
        @test evaluate_and_jacobian!(u, U, F, v, p) === nothing
        @test u ≈ û
        @test U ≈ Û
        v, V = evaluate_and_jacobian(F, v, p)
        @test v ≈ û
        @test V ≈ Û
    end

    @testset "Homotopy" begin
        ModelKit.@var x y z t

        function H_test(x, s)
            let (x, y, z, t) = (x[1], x[2], x[3], s)
                [x^2 + y + z + 2t, 4 * x^2 * z^2 * y + 4z - 6x * y * z^2]
            end
        end
        h = ModelKit.Homotopy(H_test([x, y, z], t), [x, y, z], t)

        @test sprint(show, h) == """
        Homotopy in t
         variables: x, y, z

         x ^ 2 + y + z + 2t
         (4 * x ^ 2 * z ^ 2 * y + 4z) - 6 * x * y * z ^ 2"""

        @test h == ModelKit.Homotopy(H_test([x, y, z], t), [x, y, z], t)

        H = ModelKit.compile(h)
        @test typeof(H) == typeof(ModelKit.compile(h))
        @test size(H) == (2, 3)
        @test size(H, 1) == 2
        @test size(H, 2) == 3
        @test length(H) == 2

        u = [0.0, 0.0]
        U = zeros(2, 3)
        v = [0.4814, -0.433, -0.82709]
        s = 0.359

        fdm = FD.central_fdm(5, 1)
        û = H_test(v, s)
        @test evaluate(H, v, s) ≈ û
        @test H(v, s) ≈ û
        @test evaluate(ModelKit.interpreted(H), v, s) ≈ û
        @test ModelKit.interpreted(H)(v, s) ≈ û
        @test evaluate!(u, H, v, s) ≈ û
        @test u ≈ û
        Û = FD.jacobian(fdm, v -> H_test(v, s), v)[1]
        @test jacobian(H, v, s) ≈ Û
        @test jacobian(ModelKit.interpreted(H), v, s) ≈ Û
        @test jacobian!(U, H, v, s) ≈ Û
        @test U ≈ Û rtol = 1e-8
        u .= 0
        U .= 0
        @test evaluate_and_jacobian!(u, U, H, v, s) === nothing
        @test u ≈ û
        @test U ≈ Û

        û = FD.grad(fdm, s -> H_test(v, s), s)[1]
        @test dt(H, v, s) ≈ û
        @test dt(ModelKit.interpreted(H), v, s) ≈ û
        @test dt!(u, H, v, s) ≈ û
        u .= 0
        U .= 0
        @test jacobian_and_dt!(u, U, H, v, s) === nothing
        @test u ≈ û
        @test U ≈ Û
    end

    @testset "Homotopy (parameters)" begin
        ModelKit.@var x y z t a b c

        function H_test(x, s, p)
            let (x, y, z, t, a, b, c) = (x[1], x[2], x[3], s, p[1], p[2], p[3])
                [
                 x^2 + y + z + a^2 - b + 2t,
                 4 * x^2 * z^2 * a * b * y + c * 4z - 6x * y * z^2,
                ]
            end
        end

        h = ModelKit.Homotopy(H_test([x, y, z], t, [a, b, c]), [x, y, z], t, [a, b, c])
        @test sprint(show, h) == """
        Homotopy in t
         variables: x, y, z
         parameters: a, b, c

         ((x ^ 2 + y + z + a ^ 2) - b) + 2t
         (4 * x ^ 2 * z ^ 2 * a * b * y + c * 4 * z) - 6 * x * y * z ^ 2"""

        H = ModelKit.compile(h)
        @test startswith(
            sprint(show, H),
            "Compiled{THomotopy{2,3,3,#",
        )

        @test size(H) == (2, 3)
        @test size(H, 1) == 2
        @test size(H, 2) == 3
        @test length(H) == 2

        u = [0.0, 0.0]
        U = zeros(2, 3)
        v = [0.4814, -0.433, -0.82709]
        s = 0.359
        p = [-2.123, 0.412, 1.21321]

        fdm = FD.central_fdm(5, 1)
        û = H_test(v, s, p)
        @test evaluate(H, v, s, p) ≈ û
        @test H(v, s, p) ≈ û
        @test evaluate(ModelKit.interpreted(H), v, s, p) ≈ û
        @test ModelKit.interpreted(H)(v, s, p) ≈ û
        @test evaluate!(u, H, v, s, p) ≈ û
        @test u ≈ û
        Û = FD.jacobian(fdm, v -> H_test(v, s, p), v)[1]
        @test jacobian(H, v, s, p) ≈ Û
        @test jacobian(ModelKit.interpreted(H), v, s, p) ≈ Û
        @test jacobian!(U, H, v, s, p) ≈ Û
        @test U ≈ Û rtol = 1e-8
        u .= 0
        U .= 0
        @test evaluate_and_jacobian!(u, U, H, v, s, p) === nothing
        @test u ≈ û
        @test U ≈ Û

        û = FD.grad(fdm, s -> H_test(v, s, p), s)[1]
        @test dt(H, v, s, p) ≈ û
        @test dt(ModelKit.interpreted(H), v, s, p) ≈ û
        @test dt!(u, H, v, s, p) ≈ û
        u .= 0
        U .= 0
        @test jacobian_and_dt!(u, U, H, v, s, p) === nothing
        @test u ≈ û
        @test U ≈ Û

        v, V = jacobian_and_dt(H, v, s, p)
        @test v ≈ û
        @test V ≈ Û
    end

    @testset "Modeling" begin
        @testset "Bottleneck" begin
            ModelKit.@var x y z
            f = [(0.3 * x^2 + 0.5z + 0.3x + 1.2 * y^2 - 1.1)^2 +
                 (0.7 * (y - 0.5x)^2 + y + 1.2 * z^2 - 1)^2 - 0.3]

            I = let
                x = ModelKit.variables(f)
                n, m = length(x), length(f)
                ModelKit.@unique_var y[1:n] v[1:m] w[1:m]
                J = [ModelKit.differentiate(fᵢ, xᵢ) for fᵢ in f, xᵢ in x]
                f′ = [ModelKit.subs(fᵢ, x => y) for fᵢ in f]
                J′ = [ModelKit.subs(gᵢ, x => y) for gᵢ in J]
                Nx = (x - y) - J' * v
                Ny = (x - y) - J′' * w
                ModelKit.System([f; f′; Nx; Ny], [x; y; v; w])
            end
            @test I isa ModelKit.System
            @test size(I) == (8, 8)
        end

        @testset "Steiner" begin
            ModelKit.@var x[1:2] a[1:5] c[1:6] y[1:2, 1:5]

            #tangential conics
            f = sum([a; 1] .* ModelKit.monomials(x, 2))
            ∇ = ModelKit.differentiate(f, x)
            #5 conics
            g = sum(c .* ModelKit.monomials(x, 2))
            ∇_2 = ModelKit.differentiate(g, x)
            #the general system
            #f_a_0 is tangent to g_b₀ at x₀
            function Incidence(f, a₀, g, b₀, x₀)
                fᵢ = f(x => x₀, a => a₀)
                ∇ᵢ = [∇ᵢ(x => x₀, a => a₀) for ∇ᵢ in ∇]
                Cᵢ = g(x => x₀, c => b₀)
                ∇_Cᵢ = [∇ⱼ(x => x₀, c => b₀) for ∇ⱼ in ∇_2]

                [fᵢ; Cᵢ; det([∇ᵢ ∇_Cᵢ])]
            end
            ModelKit.@var v[1:6, 1:5]
            I = vcat(map(i -> Incidence(f, a, g, v[:, i], y[:, i]), 1:5)...)
            F = ModelKit.System(I, [a; vec(y)], vec(v))
            @test size(F) == (15, 15)
        end

        @testset "Reach plane curve" begin
            ModelKit.@var x y
            f = (x^3 - x * y^2 + y + 1)^2 * (x^2 + y^2 - 1) + y^2 - 5
            ∇ = ModelKit.differentiate(f, [x; y]) # the gradient
            H = ModelKit.differentiate(∇, [x; y]) # the Hessian

            g = ∇ ⋅ ∇
            v = [-∇[2]; ∇[1]]
            h = v' * H * v
            dg = ModelKit.differentiate(g, [x; y])
            dh = ModelKit.differentiate(h, [x; y])

            ∇σ = g .* dh - ((3 / 2) * h) .* dg

            F = ModelKit.System([v ⋅ ∇σ; f], [x, y])
            @test size(F) == (2, 2)
        end
    end
end
