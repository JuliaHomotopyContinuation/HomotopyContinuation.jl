@testset "StraightLineHomotopy" begin
    TP.@polyvar x y a
    H = StraightLineHomotopy([y^2-x], [x^2+3.0y])
    @test typeof(H)<:StraightLineHomotopy{Float64}
    @test_throws ErrorException StraightLineHomotopy([x^2+3.1y], [a])
    @test_throws ErrorException StraightLineHomotopy([x^2+3y, y^2-x], [x^2+3y])

    @test evaluate(H, [2.0, 1.0], 1.0) == [-1.0]
    @test H([2.0, 1.0], 1.0) == [-1.0]
    @test H([2.0, 1.0], 0.0) == [7.0]

    J_H = differentiate(H)

    @test J_H([1.0, 1.0], 0.0) == [2 3]
    @test J_H([1.0, 1.0], 1.0) == [-1 2]

    ∂H∂t = dt(H)

    ## time derivate is independent of t
    @test ∂H∂t([1.0,1.0], 0.23) == [-4]

    @test degrees(H) == [2]
    @test startsystem(H) == PolySystem([y^2-x])
    @test targetsystem(H) == PolySystem([x^2+3.0y])

    H = StraightLineHomotopy([y^2-x, 3x + y], [x^2+3y, 2y-3x])

    @test homogenized(H) == false
    K = homogenize(H)
    @test ishomogenous(H) == false
    @test homogenized(K)
    @test ishomogenous(K)
    @test nvariables(K) == 3
end
