using HomotopyContinuation
using Base.Test

@testset "straight_line" begin

    @testset "interface" begin
        x, y = MPoly.generators(Float64, :x, :y)
        f = MPoly.system(x^2+3y)
        g = MPoly.system(y^2-x)

        H = StraightLineHomotopy(g, f)
        @test typeof(H)<:StraightLineHomotopy{Float64}
        a = MPoly.generators(Complex128, :a)
        @test_throws MethodError StraightLineHomotopy(g, MPoly.system(a))
        @test_throws ErrorException StraightLineHomotopy(MPoly.system(x^2+3y, y^2-x), MPoly.system(x^2+3y))

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
    end

    @testset "homogenize" begin
        x, y = MPoly.generators(Float64, :x, :y)
        f = MPoly.system(x^2+3y, 2y-3x)
        g = MPoly.system(y^2-x, 3x + y)

        H = StraightLineHomotopy(g, f)

        K = homogenize(H)
        
        @test nvars(K) == 3
    end
end