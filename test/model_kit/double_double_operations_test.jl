@testset "Fused double-double polynomial operations" begin
    D = HC.DoubleDouble.DoubleF64
    C = HC.DoubleDouble.ComplexDF64
    rng = MersenneTwister(0xf053)
    function sample_dd()
        hi = ldexp(randn(rng), rand(rng, -200:200))
        D(hi, 0.25 * randn(rng) * eps(hi))
    end
    setprecision(BigFloat, 256) do
        for _ = 1:200
            x, y, z, w = ntuple(_ -> C(sample_dd(), sample_dd()), 4)
            for (a, b, c, d) in ((x, y, z, w), (x, y, -x, y), (x, y, -x * y, one(C)))
                A, B, Cb, Db = Complex{BigFloat}.((a, b, c, d))
                checks = (
                    (ModelKit.op_mul(a, b), A * B, abs(A) * abs(B)),
                    (ModelKit.op_muladd(a, b, c), A * B + Cb, abs(A) * abs(B) + abs(Cb)),
                    (ModelKit.op_mulsub(a, b, c), A * B - Cb, abs(A) * abs(B) + abs(Cb)),
                    (ModelKit.op_submul(a, b, c), Cb - A * B, abs(A) * abs(B) + abs(Cb)),
                    (
                        ModelKit.op_mulmuladd(a, b, c, d),
                        A * B + Cb * Db,
                        abs(A) * abs(B) + abs(Cb) * abs(Db),
                    ),
                    (
                        ModelKit.op_mulmulsub(a, b, c, d),
                        A * B - Cb * Db,
                        abs(A) * abs(B) + abs(Cb) * abs(Db),
                    ),
                    (ModelKit.op_sqr(a), A^2, abs(A)^2),
                    (ModelKit.op_cb(a), A^3, abs(A)^3),
                    (ModelKit.op_mul3(a, b, c), A * B * Cb, abs(A) * abs(B) * abs(Cb)),
                    (
                        ModelKit.op_mul4(a, b, c, d),
                        A * B * Cb * Db,
                        abs(A) * abs(B) * abs(Cb) * abs(Db),
                    ),
                )
                for (actual, expected, scale) in checks
                    @test abs(Complex{BigFloat}(actual) - expected) <=
                          64 * big(eps(Float64))^2 * scale
                end
            end
        end
    end
    # Low limbs survive cancellation even when the high product terms cancel.
    δ = D(2.0^-60)
    x, y = C(1 + δ, 1), C(1, 1)
    @test real(ModelKit.op_mul(x, y)) == δ
    @test ModelKit.op_mul(zero(C), one(C)) == zero(C)
end
