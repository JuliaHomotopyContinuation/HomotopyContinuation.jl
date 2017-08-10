@testset "straight_line" begin

    @TP.polyvar x y a
    f = [x^2+3.0*y]
    g = [y^2-1.0*x]

    H = StraightLineHomotopy(g, f)
    @test typeof(H)<:StraightLineHomotopy{Float64}
    @test_throws MethodError StraightLineHomotopy(g, [TP.polynomial(a)])
    @test_throws ErrorException StraightLineHomotopy([x^2+3y, y^2-x], [x^2+3y])

    @test evaluate(H, [2.0, 1.0], 1.0) == [-1.0]
    @test evaluate(H, [2.0, 1.0], 0.0) == [7.0]

    J_H = jacobian(H)

    @test J_H([1.0, 1.0], 0.0) == [2 3]
    @test J_H([1.0, 1.0], 1.0) == [-1 2]

    ∂H∂t = dt(H)

    ## time derivate is independent of t
    @test ∂H∂t([1.0,1.0], 0.23) == [-4]

    @test degrees(H) == [2]
    @test startsystem(H) == g
    @test targetsystem(H) == f

    @testset "homogenize" begin
        @TP.polyvar x y z
        f = [x^2+3y, 2y-3x]
        g = [y^2-x, 3x + y]

        H = StraightLineHomotopy(g, f)

        K = homogenize(H, z)
        
        @test nvars(K) == 3
    end
end