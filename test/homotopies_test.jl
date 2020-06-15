@testset "Homotopies" begin
    @testset "ParameterHomotopy" begin
        @var x y a b c

        f = [(2 * x^2 + b^2 * y^3 + 2 * a * x * y)^3, (a + c)^4 * x + y^2]
        F = System(f, [x, y], [a, b, c])

        p = [5.2, -1.3, 9.3]
        q = [2.6, 3.3, 2.3]

        H = HC.ParameterHomotopy(F, p, q)

        v = [0.192, 2.21]
        t = 0.232

        u = zeros(ComplexF64, 2)
        U = zeros(ComplexF64, 2, 2)

        HC.evaluate!(u, H, v, t)
        @test u ≈ f([x, y] => v, [a, b, c] => t * p + (1 - t) * q)

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
end
