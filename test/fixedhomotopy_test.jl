@testset "FixedHomotopy" begin
    PolyImpl.@polyvar x y a
    H = StraightLineHomotopy([y^2-x], [x^2+3.0y])

    H1 = FixedHomotopy(H, 1.0)
    H0 = FixedHomotopy(H, 0.0)
    HH = StraightLineHomotopy(H1, H0)

    @test typeof(H1)<:FixedHomotopy{Float64, typeof(H), Float64}

    @test evaluate(H1, [2.0, 1.0]) == evaluate(H, [2.0, 1.0], 1.0)
    @test H1([2.0, 1.0]) == H([2.0, 1.0], 1.0)
    @test H0([2.0, 1.0]) == H([2.0, 1.0], 0.0)
    @test HH([2.0, 1.0], 1.0) == H([2.0, 1.0], 1.0)
    @test HH([2.0, 1.0], 0.0) == H([2.0, 1.0], 0.0)

    J_H = differentiate(H)
    J_H1 = differentiate(H1)
    J_H0 = differentiate(H0)

    @test J_H0([1.0, 1.0]) == J_H([1.0, 1.0], 0.0)
    @test J_H1([1.0, 1.0]) == J_H([1.0, 1.0], 1.0)

    @test degrees(H1) == [2]

    H = FixedHomotopy(StraightLineHomotopy([y^2-x, 3x + y], [x^2+3y, 2y-3x]), 0.34)

    @test homogenized(H) == false
    K = homogenize(H)
    @test ishomogenous(H) == false
    @test homogenized(K)
    @test ishomogenous(K)
    @test nvariables(K) == 3
end
