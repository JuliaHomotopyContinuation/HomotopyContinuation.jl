using HomotopyContinuation.DoubleDouble
const DD = HomotopyContinuation.DoubleDouble

@testset "DoubleF64" begin
    @test DoubleF64(Float32(2.0)) == DoubleF64(2.0)
    @test DoubleF64(Float16(2.0)) == DoubleF64(2.0)
    @test DoubleF64(2) == DoubleF64(2.0)
    @test DoubleF64(big(π)) isa DoubleF64
    @test DoubleF64(π) == DoubleF64(big(π))

    @test DD.hi(DoubleF64(2.0)) == 2.0
    @test DD.lo(DoubleF64(2.0)) == 0.0

    @test isbitstype(DoubleF64)
    @test isbits(DoubleF64(3))

    @test zero(DoubleF64(2.0)) == zero(DoubleF64)

    @test convert(Integer, DoubleF64(2.0)) isa Int64
    @test 32 ≤ length(string(DoubleF64(rand()))) ≤ 40

    @test convert(DoubleF64, 2.0) isa DoubleF64
    @test convert(DoubleF64, Float32(2.0)) isa DoubleF64
    @test convert(DoubleF64, Float16(2.0)) isa DoubleF64
    @test convert(DoubleF64, 2) isa DoubleF64
    @test convert(DoubleF64, π) isa DoubleF64

    @test promote_type(DoubleF64, Float64) == DoubleF64
    @test promote_type(DoubleF64, Float32) == DoubleF64
    @test promote_type(DoubleF64, Float16) == DoubleF64

    @test DD.wide_sub(rand(), rand()) isa DoubleF64
    @test DD.wide_add(rand(), rand()) isa DoubleF64
    @test DD.wide_div(rand(), rand()) isa DoubleF64
    @test DD.wide_mul(rand(), rand()) isa DoubleF64

    x = DoubleF64(2.321)
    @test x isa DoubleF64
    @test x / 3.2 isa DoubleF64

    @test x + 2 isa DoubleF64
    @test 2 + x isa DoubleF64

    @test im * DoubleF64(rand()) isa Complex{DoubleF64}
    @test im * DoubleF64(rand()) isa ComplexDF64
    @test x < 2 * x
    @test x ≤ x
    @test convert(Float64, x) ≤ x
    @test x ≤ 2 * convert(Float64, x)
    @test DoubleF64(2.0) == 2.0
    @test 2.0 == DoubleF64(2.0)

    @test isnan(DoubleF64(NaN))
    @test isinf(DoubleF64(Inf))

    @test floor(DoubleF64(3.2)) == 3.0
    @test floor(Int, DoubleF64(3.2)) == 3
    @test ceil(DoubleF64(3.2)) == 4.0
    @test ceil(Int, DoubleF64(3.2)) == 4
    @test trunc(DoubleF64(3.2)) == 3.0
    @test trunc(Int, DoubleF64(3.2)) == 3
    @test isinteger(DoubleF64(3.2)) == false
    @test isinteger(DoubleF64(3.0)) == true

    @test DoubleF64(0.0)^0.14 == 1.0

    for k = 1:10
        x = DoubleF64(rand()) * 20 - 10
        y = DoubleF64(rand()) * 20 - 10
        @test x * y ≈ big(x) * big(y) atol = 1e-29
        @test x + y ≈ big(x) + big(y) atol = 1e-30
        @test x - y ≈ big(x) - big(y) atol = 1e-30
        @test x / y ≈ big(x) / big(y) atol = 1e-26
        u = rand()
        @test x / u ≈ big(x) / big(u) atol = 1e-26
        @test x^5 ≈ BigFloat(x)^5 atol = 1e-25

        @test exp(x) ≈ exp(big(x)) atol = 1e-26
        @test log(abs(x)) ≈ log(big(abs(x))) atol = 1e-29
        @test sin(x) ≈ sin(big(x)) atol = 1e-30
        @test cos(x) ≈ cos(big(x)) atol = 1e-30
        @test atan(x) ≈ atan(big(x)) atol = 1e-30
        @test atan(y, x) ≈ atan(big(y), big(x)) atol = 1e-30
        @test tan(x) ≈ tan(big(x)) atol = 1e-21
        sin_x, cos_x = DD.sincos(x)
        @test sin_x ≈ sin(big(x)) atol = 1e-30
        @test cos_x ≈ cos(big(x)) atol = 1e-30

        z = DoubleF64(rand() * 2 - 1)
        @test asin(z) ≈ asin(big(z)) atol = 1e-30
        @test acos(z) ≈ acos(big(z)) atol = 1e-30

        n, p, d = Base.decompose(x)
        @test x ≈ BigFloat(n) * BigFloat(2)^p / BigFloat(d)
    end
end
