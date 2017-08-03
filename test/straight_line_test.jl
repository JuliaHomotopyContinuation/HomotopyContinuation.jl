using HomotopyContinuation
using Base.Test

@testset "straight_line" begin

    @testset "interface" begin
        x, y = MPoly.generators(Float64, :x, :y)
        f = MPoly.system(x^2+3y)
        g = MPoly.system(y^2-x)

        H = StraightLineHomotopy(g, f)
        @test typeof(H)<:StraightLineHomotopy{Float64}

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


    #     H = StraightLineHomotopy(g, f)
    #     @test typeof(H)<:StraightLineHomotopy{Float64}

    #     @test evaluate(H, [2.0, 1.0], 1.0) == [-1.0]
    #     @test evaluate(H, [2.0, 1.0], 0.0) == [7.0]

    #     J_H = jacobian(H)

    #     @test J_H([1.0, 1.0], 0.0) == [2 3]
    #     @test J_H([1.0, 1.0], 1.0) == [-1 2]

    #     ∂H∂t = dt(H)

    #     ## time derivate is independent of t
    #     @test ∂H∂t([1.0,1.0], 0.23) == [-4]

    #     @test degrees(H) == [2]
    #     @test startsystem(H) == g
    #     @test targetsystem(H) == f
    end

    # @testset "eval" begin
    #     x, y = MPoly.generators(Float64, :x, :y)
    #     f = MPoly.system(3x^2*y^2+4x^3*y+2x*y^2)
    #     @test evaluate(f, [3.0, 5.0]) == [1365.0]
    # end

    # @testset "jacobian" begin
    #     x, y = MPoly.generators(Float64, :x, :y)
    #     f = polysystem(3x^2*y^2+4x^3*y+2x*y^2)

    #     df = jacobian(f)
    #     @test df([2.0,2.0]) == [152.0 96.0]
    # end

    # @testset "is_homogenous" begin
    #     x, y = MPoly.generators(Float64, :x, :y)

    #     @test is_homogenous(MPoly.system(3x^2*y+4x^3+2x*y^2)) == true
    #     @test is_homogenous(MPoly.system(3x^2*y+4x^3+2x*y)) == false
    # end

    # @testset "degrees, nequations, nvars" begin
    #     x, y = MPoly.generators(Float64, :x, :y)
    #     F = MPoly.system(x^4, 2x+2y^2,3x*y^2)

    #     @test degrees(F) == [4,2,3]
    #     @test nequations(F) == 3
    #     @test nvars(F) == 2
    # end

    # @testset "weyl dot" begin
    #     x,y = MPoly.generators(Complex128, :x,:y)

    #     f = 3.0x^2 + 2.0x*y - 1.0*y^2
    #     g = (-2.5+2im)x^2 + -3.0x*y + 4.0*y^2

    #     F = MPoly.system(f)
    #     G = MPoly.system(g)

    #     @test weyl_dot(F,G) == 3.0 * conj(-2.5+2im) + 2.0 * (-3.0) / 2 + (-1.0) * 4.0
    #     @test weyl_dot(F,F) == 9.0 + 4.0  / 2 + 1.0

    #     @test weyl_dot(F,F) ≈ weyl_norm(F)^2
    # end
end