using HomotopyContinuation
using Base.Test

@testset "straight_line" begin

    @testset "interface" begin
        x, y = MPoly.generators(Complex128, :x, :y)
        f = MPoly.system(x^2+3y)
        g = MPoly.system(y^2-x)

        H = GammaTrickHomotopy(g, f, 1.0+0im)
        @test typeof(H)<:GammaTrickHomotopy{Complex128}

        @test evaluate(H, [2.0+0im, 1.0], 1.0) == [-1.0+0im]
        @test evaluate(H, [2.0+0im, 1.0], 0.0) == [7.0+0im]

        J_H = jacobian(H)

        @test J_H([1.0+0im, 1.0], 0.0) == [2 3+0im]
        @test J_H([1.0+0im, 1.0], 1.0) == [-1 2+0im]

        ∂H∂t = dt(H)

        ## time derivate is independent of t
        @test ∂H∂t([1.0+0im,1.0], 0.23) == [-4+0im]

        @test degrees(H) == [2]
        @test startsystem(H) == g
        @test targetsystem(H) == f
    end

    @testset "homogenize" begin
        x, y = MPoly.generators(Complex128, :x, :y)
        f = MPoly.system(x^2+3y, 2y-3x)
        g = MPoly.system(y^2-x, 3x + y)

        H = GammaTrickHomotopy(g, f)

        K = homogenize(H)
        
        @test nvars(K) == 3
        @test H.γ == K.γ
    end

    @testset "constructor" begin
        x, y = MPoly.generators(Complex128, :x, :y)
        f = MPoly.system(x^2+3y)
        g = MPoly.system(y^2-x)

        H = GammaTrickHomotopy(g, f)
        @test norm(H.γ) ≈ 1.0 atol=1e-8

        H = GammaTrickHomotopy(g, f, 1.3+4im)
        @test H.γ ≈ 1.3+4im

    end

end