@testset "ModelKit - Symbolic" begin
    @testset "SymEngine" begin
        @test Expression(MathConstants.catalan) isa Expression
        @test Expression(MathConstants.e) isa Expression
        @test Expression(MathConstants.pi) isa Expression
        @test Expression(MathConstants.γ) isa Expression
        @test Expression(Float16(2.1)) isa Expression
        @test Expression(Float32(2.1)) isa Expression
        @test Expression(big(2.1)) isa Expression
        @test Expression(Int32(2)) isa Expression
        @test Expression(Int16(2)) isa Expression
        @test Expression(Int128(2)) isa Expression
        @test Expression(2 // 3) isa Expression
        @test convert(BigFloat, Expression(big(2.1))) == big(2.1)
        @test ModelKit.is_number(Expression(2))
        @test ModelKit.to_number(Expression(big(2.1))) isa BigFloat
        @test ModelKit.to_number(Expression(2 // 3)) isa Rational{Int32}
        @test ModelKit.to_number(Expression(2) // Expression(3)) isa Rational{Int32}
        @test ModelKit.to_number(Expression(big(2)^129)) isa BigInt

        @test BigFloat(Expression(big(2.1))) == big(2.1)
        @test Float64(Expression(2)) == 2.0
        @test Rational{Int32}(Expression(2 // 3)) == 2 // 3
        @test UInt(Expression(2)) == UInt(2)

        @test complex(Expression(1), Expression(2)) == Expression(1) + im * 2
        @test complex(Expression(1), 2) == Expression(1) + im * 2
        @test complex(1, Expression(2)) == Expression(1) + im * 2
        @test +Expression(2) == 2
        @test Expression(HC.DoubleDouble.DoubleF64(2)) == big(2.0)
    end

    @testset "Variables" begin
        @var a b x[1:2] y[1:2, 1:3]

        @test a isa Variable
        @test b isa Variable
        @test x isa Vector{Variable}
        @test length(x) == 2
        @test y isa Matrix{Variable}
        @test size(y) == (2, 3)
        @test join(x, " ") == "x₁ x₂"

        @test reshape(sort(vec(y)), 2, 3) == y
        @var x[9:11]
        @test sort(x) == x

        @var c, d
        @test c isa Variable
        @test d isa Variable
        @test variables(c) == [c]
        @test c === c
        @test c !== copy(c)

        let
            c2 = @unique_var c
            @test c != c2
        end

        f = a + b
        @test variables(f) == [a, b]
        @test nvariables(f) == 2
        g = x[1]
        @test variables([f, g]) == [a, b, x[1]]
        @test nvariables([f, g]) == 3
    end

    @testset "Subs" begin
        @var x y w z u
        f = x^2 * (x + y * w)

        @test subs(f, x => z) == z^2 * (z + w * y)
        @test subs([f], x => z) == [z^2 * (z + w * y)]
        @test subs(f, [x, y] => [z^2, z + 2]) == z^4 * (w * (2 + z) + z^2)
        @test subs(f, [x, y] => [z^2, z + 2], w => u) == z^4 * (u * (2 + z) + z^2)
        @test subs(f, x => z^2, y => 3, w => u) == z^4 * (3 * u + z^2)
        @test subs(f, Dict(x => z^2, y => 3, w => u)) == z^4 * (3 * u + z^2)
    end

    @testset "Evaluation" begin
        @var x y w
        f = x^2 * (x + y * w)
        @test evaluate([f, f], [x, y, w] => [2, 3, -5]) == [-52, -52]
        @test [f, 2f]([x, y, w] => [2, 3, -5]) == [-52, -104]
    end

    @testset "Linear Algebra" begin
        @var x[1:2, 1:2]
        @test det(x[1:1, 1:1]) == x[1, 1]
        @test det(x) == -x[2, 1] * x[1, 2] + x[2, 2] * x[1, 1]
        @test det(x') == -x[2, 1] * x[1, 2] + x[2, 2] * x[1, 1]
        @test det(transpose(x)) == -x[2, 1] * x[1, 2] + x[2, 2] * x[1, 1]
    end

    @testset "Differentiation" begin
        @var x y

        f = x^2 + y^2
        g = x^3 + y^3

        @test differentiate(f, x) == 2x
        @test differentiate(f, [x, y]) == [2x, 2y]

        @test differentiate([f, g], x) == [2x, 3 * x^2]
        @test differentiate([f, g], [x, y]) == [2x 2y; 3*x^2 3*y^2]
    end

    @testset "Expand" begin
        @var x y
        @test expand((x + y)^2) == 2 * x * y + x^2 + y^2
    end

    @testset "to_dict" begin
        @var x y a
        @test ModelKit.to_dict(x^2 + y * x * a + y * x, [x, y]) ==
              Dict([2, 0] => Expression(1), [1, 1] => a + 1)
        @test ModelKit.degree(x * y^5 + y^2, [x, y]) == 6
    end

    @testset "Modeling" begin
        @testset "Bottleneck" begin
            @var x y z
            f = [
                (0.3 * x^2 + 0.5z + 0.3x + 1.2 * y^2 - 1.1)^2 +
                (0.7 * (y - 0.5x)^2 + y + 1.2 * z^2 - 1)^2 - 0.3,
            ]

            I = let x = [x, y, z]
                n, m = length(x), length(f)
                @unique_var y[1:n] v[1:m] w[1:m]
                J = differentiate(f, x)
                f′ = f(x => y)
                J′ = J(x => y)
                Nx = (x - y) - J' * v
                Ny = (x - y) - J′' * w
                System([f; f′; Nx; Ny], [x; y; v; w])
            end
            @test I isa System
            @test size(I) == (8, 8)
        end

        @testset "Steiner" begin
            @var x[1:2] a[1:5] c[1:6] y[1:2, 1:5]

            #tangential conics
            f = sum([a; 1] .* monomials([x; 1], 2; affine = false))
            ∇ = differentiate(f, x)
            #5 conics
            g = sum(c .* monomials([x; 1], 2; affine = false))
            ∇_2 = differentiate(g, x)
            #the general system
            #f_a_0 is tangent to g_b₀ at x₀
            function Incidence(f, a₀, g, b₀, x₀)
                fᵢ = f(x => x₀, a => a₀)
                ∇ᵢ = ∇(x => x₀, a => a₀)
                Cᵢ = g(x => x₀, c => b₀)
                ∇_Cᵢ = ∇_2(x => x₀, c => b₀)

                [fᵢ; Cᵢ; det([∇ᵢ ∇_Cᵢ])]
            end
            @var v[1:6, 1:5]
            I = vcat(map(i -> Incidence(f, a, g, v[:, i], y[:, i]), 1:5)...)
            F = System(I, [a; vec(y)], vec(v))
            @test size(F) == (15, 15)
        end

        @testset "Reach plane curve" begin
            @var x y
            f = (x^3 - x * y^2 + y + 1)^2 * (x^2 + y^2 - 1) + y^2 - 5
            ∇ = differentiate(f, [x; y]) # the gradient
            H = differentiate(∇, [x; y]) # the Hessian

            g = ∇ ⋅ ∇
            v = [-∇[2]; ∇[1]]
            h = v' * H * v
            dg = differentiate(g, [x; y])
            dh = differentiate(h, [x; y])

            ∇σ = g .* dh - ((3 / 2) * h) .* dg

            F = System([v ⋅ ∇σ; f], [x, y])
            @test size(F) == (2, 2)
        end
    end

    @testset "Horner" begin
        # use fano scheme as test
        n = 4
        s = 1
        @var x[1:n]
        ν = monomials(x, 5; affine = false)
        N = length(ν)
        q₀ = randn(ComplexF64, N - 1)
        @var q[1:N-1]
        F = sum(q[i] * ν[i] for i = 1:N-1) + 1

        @var a[1:n-1] b[1:n-1] t
        L = [a; 1] .* t + [b; 0]

        G = subs(F, x => L)

        FcapL = collect(values(ModelKit.to_dict(expand(G), [t])))

        @test all(FcapL) do f
            g = horner(f)
            expand(g) == f
        end
        @test all(FcapL) do f
            g = horner(f, [a; b])
            expand(g) == f
        end
    end

    @testset "Rand / dense poly" begin
        @var x y
        f, c = dense_poly([x, y], 3)
        @test length(c) == 10
        g = rand_poly(Float64, [x, y], 3)
        @test subs(f, c => coefficients(g, [x, y])) == g
        _, coeffs = exponents_coefficients(g, [x, y]; expanded = true)
        @test subs(f, c => coeffs) == g

        @var x[1:3]
        f, c = dense_poly(x, 3, coeff_name = :c)
        g = x[1]^3 + x[2]^3 + x[3]^3 - 1
        gc = coeffs_as_dense_poly(g, x, 3)
        @test subs(f, c => gc) == g
    end

    @testset "System" begin
        @var x y a b
        f = [(x + y)^3 + x^2 + x + 5y + 3a, 2 * x^2 + b]
        F = System(f, [x, y], [b, a])
        @test nvariables(F) == 2
        @test nparameters(F) == 2
        @test variables(F) == [x, y]
        @test parameters(F) == [b, a]
        show_F = """
        System of length 2
         2 variables: x, y
         2 parameters: b, a

         3*a + x + 5*y + x^2 + (x + y)^3
         b + 2*x^2"""
        @test sprint(show, F) == show_F
        @test degrees(F) == [3, 2]

        T = CompiledSystem(F; optimizations = false)
        F2 = System(T)
        @test F == F2
        @test sprint(show, T) == "Compiled: $show_F"
        @test size(T) == size(F) == (2, 2)

        T2 = CompiledSystem(F; optimizations = true)
        F3 = System(T2)
        @test expand.(F.expressions) == expand.(F3.expressions)

        F4 = System(f, parameters = [b, a])
        @test variables(F4) == [x, y]

        @var x y a b
        f = Any[(x+y)^3+x^2+x+5y+3a, 2*x^2+b]
        F = System(f, [x, y], [b, a])
        @test F isa System
    end

    @testset "System variables groups + homogeneous" begin
        @var x y z
        f = System([
            (x^2 + y^2 + x * y - 3 * z^2) * (x + 3z),
            (x^2 + y^2 + x * y - 3 * z^2) * (y - x + 2z),
            2x + 5y - 3z,
        ])
        @test is_homogeneous(f)
        multi_degrees(f) == [1 1; 2 0]

        @var x y z v w
        g = System([x * y - 2v * w, x^2 - 4 * v^2], variable_groups = [[x, v], [y, w]])
        @test is_homogeneous(g)
        multi_degrees(g) == [1 1; 2 0]
    end

    @testset "Homotopy" begin
        @var x y z t

        h = [x^2 + y + z + 2t, 4 * x^2 * z^2 * y + 4z - 6x * y * z^2]
        H = Homotopy(h, [x, y, z], t)

        @test nvariables(H) == 3
        @test nparameters(H) == 0
        @test variables(H) == [x, y, z]
        @test parameters(H) == []
        show_H = """
        Homotopy in t of length 2
         3 variables: x, y, z

         2*t + y + z + x^2
         4*z - 6*x*y*z^2 + 4*x^2*y*z^2"""

        @test sprint(show, H) == show_H

        H_any = Homotopy(convert(Vector{Any}, h), [x, y, z], t)
        @test H_any isa Homotopy

        T = CompiledHomotopy(H)
        H2 = ModelKit.interpret(T)
        @test H == H2
        @test sprint(show, T) == "Compiled: $show_H"
        @test size(T) == size(H) == (2, 3)
    end

    @testset "Optimizations for rational system" begin
        @var y[1:6] q[1:8]
        y₁, y₂, y₃, y₄, y₅, y₆ = y
        q₁, q₂, q₃, q₄, q₅, q₆, q₇, q₈ = q

        f =
            q₁ / y₁ - q₂ / (-y₁ + y₂) +
            q₅ * y₄ / (y₁ * y₄ - y₃ * y₂) +
            q₈ * y₆ / (y₁ * y₆ - y₂ * y₅)
        F = ModelKit.System([f], y, q)
        I = ModelKit.interpreter(ModelKit.optimize(F))
        u = Expression[0]
        ModelKit.execute!(u, I, y, q, ModelKit.InterpreterCache(Expression, I))

        @test expand(u[1] - f) == 0
    end
end
