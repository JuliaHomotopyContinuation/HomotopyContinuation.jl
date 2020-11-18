@testset "ModelKit - SLP" begin

    @testset "CompiledHomotopy/InterpretedHomotopy" begin
        @var x y a b c
        f = [(2 * x^2 + b^2 * y^3 + 2 * a * x * y)^3, (a + c)^4 * x + y^2]

        @var s sp[1:3] sq[1:3]
        g = subs(f, [a, b, c] => s .* sp .+ (1 .- s) .* sq)

        h = Homotopy(g, [x, y], s, [sp; sq])

        p = [5.2, -1.3, 9.3]
        q = [2.6, 3.3, 2.3]

        for H in [InterpretedHomotopy(h), CompiledHomotopy(h)]
            @test size(H) == (2, 2)

            v = [0.192, 2.21]
            t = 0.232

            u = zeros(ComplexF64, 2)
            U = zeros(ComplexF64, 2, 2)

            evaluate!(u, H, v, t, [p; q])
            @test u ≈ f([x, y] => v, [a, b, c] => t * p + (1 - t) * q)

            taylor!(u, Val(1), H, TaylorVector{1}(Matrix(v')), t, [p; q])
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

            evaluate_and_jacobian!(u, U, H, v, t, [p; q])
            @test U ≈
                  differentiate(f, [x, y])([x, y] => v, [a, b, c] => t * p + (1 - t) * q)
        end
    end

    @testset "Codegen helpers" begin
        @test ModelKit.sqr(3 + 2im) == (3 + 2im)^2
        @test ModelKit.sqr(3) == 3^2
    end

    @testset "Homotopy codegen (Katsura(3))" begin
        n = 3
        @var x[0:n] ẋ[0:n] ẍ[0:n] x3[0:n] t γ
        K = [
            (
                sum(x[abs(l)+1] * x[abs(m - l)+1] for l = -n:n if abs(m - l) <= n) - x[m+1] for m = 0:n-1
            )...,
            x[1] + 2 * sum(x[i+1] for i = 1:n) - 1,
        ]

        h = γ .* t .* [x[1:n] .^ 2 .- 1; x[n+1] - 1] + (1 - t) .* K
        H = ModelKit.Homotopy(h, x, t, [γ])
        TH = CompiledHomotopy(H)

        tx3 = TaylorVector{4}(zeros(Expression, 4, 4))
        y, y1, y2, y3 = vectors(tx3)
        y .= x
        y1 .= ẋ
        y2 .= ẍ
        y3 .= x3

        dt1 =
            ModelKit.taylor!(zeros(Expression, 4), Val(1), TH, TaylorVector{1}(tx3), t, [γ])
        dt2 =
            ModelKit.taylor!(zeros(Expression, 4), Val(2), TH, TaylorVector{2}(tx3), t, [γ])
        dt3 =
            ModelKit.taylor!(zeros(Expression, 4), Val(3), TH, TaylorVector{3}(tx3), t, [γ])
        dt4 =
            ModelKit.taylor!(zeros(Expression, 4), Val(4), TH, TaylorVector{4}(tx3), t, [γ])

        @var λ
        Hd1 = subs(H.expressions, t => t + λ)
        true_dt1 = subs(differentiate(Hd1, λ, 1), λ => 0)

        Hd2 = subs(H.expressions, x => x .+ λ .* ẋ, t => t + λ)
        true_dt2 = subs(differentiate(Hd2, λ, 2), λ => 0) / 2

        Hd3 = subs(H.expressions, x => x .+ λ .* ẋ .+ λ^2 .* ẍ, t => t + λ)
        true_dt3 = subs(differentiate(Hd3, λ, 3), λ => 0) / 6

        Hd4 = subs(H.expressions, x => x .+ λ .* ẋ .+ λ^2 .* ẍ .+ λ^3 .* x3, t => t + λ)
        true_dt4 = subs(differentiate(Hd4, λ, 4), λ => 0) / 24

        @test expand.(TH(x, t, [γ])) == expand.(h)
        @test expand.(ModelKit.jacobian(TH, x, t, [γ])) == expand.(differentiate(h, x))
        @test expand.(dt1) == expand.(true_dt1)
        @test expand.(dt2) == expand.(true_dt2)
        @test expand.(dt3) == expand.(true_dt3)
        @test expand.(dt4) == expand.(true_dt4)
    end

    @testset "taylor! - system" begin
        @var x[1:2] ẋ[1:2] ẍ[1:2] x3[1:2] p[1:2] ṗ[1:2]
        f = [(x[1] + x[2])^3 + x[1]^2 + x[1] + 5 * x[2] + 3 * p[1], 2 * x[1]^2 + p[2]]
        F = System(f, x, p)
        TF = CompiledSystem(F)

        tx = TaylorVector{4}([x'; ẋ'; ẍ'; x3'])
        tv = TaylorVector{5}(Expression, 2)
        v, v1, v2, v3, v4 = vectors(tv)
        @var λ

        ModelKit.taylor!(TaylorVector{3}(tv), Val(2), TF, TaylorVector{2}(tx), p)
        Fd2 = subs(F.expressions, x => x .+ λ .* ẋ)
        true_v1 = subs(differentiate(Fd2, λ, 1), λ => 0)
        true_v2 = subs(differentiate(Fd2, λ, 2), λ => 0) / 2
        @test expand.(v) == expand.(f)
        @test expand.(v1) == expand.(true_v1)
        @test expand.(v2) == expand.(true_v2)

        ModelKit.taylor!(TaylorVector{4}(tv), Val(3), TF, TaylorVector{3}(tx), p)
        Fd3 = subs(F.expressions, x => x .+ λ .* ẋ .+ λ^2 .* ẍ)
        true_v1 = subs(differentiate(Fd3, λ, 1), λ => 0)
        true_v2 = subs(differentiate(Fd3, λ, 2), λ => 0) / 2
        true_v3 = subs(differentiate(Fd3, λ, 3), λ => 0) / 6
        @test expand.(v) == expand.(f)
        @test expand.(v1) == expand.(true_v1)
        sub = variables(v2) => rand(1:10_000, 6)
        @test v2(sub) == true_v2(sub)
        @test v3(sub) == true_v3(sub)
    end

    @testset "Monomials" begin
        @var x y
        F = System([x^2 - 1, y])
        u = zeros(ComplexF64, 2)
        U = zeros(ComplexF64, 2, 2)
        w = randn(ComplexF64, 2)
        evaluate_and_jacobian!(u, U, CompiledSystem(F), w)
        @test sum(abs, U) > 0
        U .= 0
        evaluate_and_jacobian!(u, U, InterpretedSystem(F), w)
        @test sum(abs, U) > 0
    end

    @testset "Evaluation of Arb" begin
        @var x y
        # define the polynomials
        f₁ = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        f₂ = x^2 + 2x * y^2 - 2 * y^2 - 1 / 2
        F = System([f₁, f₂])
        for I in [
            Arblib.AcbMatrix(randn(ComplexF64, 2, 1)),
            Arblib.AcbRefVector(randn(ComplexF64, 2)),
        ]
            @test all(!iszero, F(I))
            @test F(real.(I)) isa Arblib.ArbVector
            @test all(!iszero, System(F)(real.(I)))

            u = Arblib.AcbVector(2)
            evaluate!(u, InterpretedSystem(System(F)), I)
            @test all(!iszero, u)
            v = Arblib.ArbVector(2)
            evaluate!(v, InterpretedSystem(System(F)), real.(I))
            @test all(!iszero, v)
            U = Arblib.AcbMatrix(2, 2)
            jacobian!(U, InterpretedSystem(System(F)), I)
            @test all(!iszero, U)
            V = Arblib.ArbMatrix(2, 2)
            jacobian!(V, InterpretedSystem(System(F)), real.(I))
            @test all(!iszero, V)
        end
    end
end
