using HomotopyContinuation
using Base.Test

@testset "MPoly" begin
    @testset "Primitives" begin

        f = MPoly.zeropoly(Float64, :x)
        @test f.vars == (:x,);

        g = MPoly.zeropoly(Float64, :x,:y)
        @test g.vars == (:x,:y);

        h = MPoly.zero(g)
        @test h.vars == (:x, :y);

        A = MPoly.zeros(h, 3, 4)
        @test size(A) == (3,4)

        f1 = MPoly.onepoly(Float64, :a)
        @test f1[0] == 1.0

        g1 = MPoly.onepoly(Float64, :x,:y)
        @test g1[0,0] == 1.0

        h1 = MPoly.one(g1)
        @test h1[0,0] == 1.0
        @test h1.vars == (:x, :y);
    end

    @testset "generators" begin
        x, y, z = MPoly.generators(Float64, :x, :y, :z)
        @test x.vars == (:x, :y, :z)
        @test x[(1,0,0)] == 1.0
        @test y[(0,1,0)] == 1.0
        @test z[(0,0,1)] == 1.0

        (a,) = MPoly.generators(Int, :a)
        @test a.vars == (:a,)
        @test a[1] == 1
    end

    @testset "embed" begin
        f = MPoly.zeropoly(Float64, :a,:y)
        f[0,2] = 3.0
        f[3,1] = 2.0

        g = MPoly.embed(f, [:a,:b,:y])
        @test g[0,0,2] == 3.0
        @test g[3,0,1] == 2.0
    end

    @testset "homogenize" begin
        f = MPoly.zeropoly(Float64, :a,:y)
        f[0,2] = 3.0
        f[3,1] = 2.0

        g = MPoly.homogenize(f, :b)
        @test g[2,0,2] == 3.0
        @test g[0,3,1] == 2.0
    end

    @testset "addition" begin
        x, y = MPoly.generators(Float64, :x, :y)
        @test (x + x)[1, 0] == 2.0
        @test (x + y)[1, 0] == 1.0
        @test (x + y)[0, 1] == 1.0

        y_ = MPoly.generator(Float64, :y)

        @test (x + y_ + y)[0, 1] == 2.0
        @test (x + y_)[0, 1] == 1.0

        @test (x + 3.2)[0, 0] == 3.2

    end

    @testset "operators + print" begin
        x, y = MPoly.generators(Float64, :x, :y)

        f = 3*x^2*y + 2x - 4
        @test string(f) == "3.0x^2y+2.0x-4.0"
        g = (-x^2+y) * y - y^2 + 1

        @test string(g) == "-x^2y+0.0y^2+1.0"
    end


    @testset "differentiated" begin
        x, y = MPoly.generators(Float64, :x, :y)

        f = 3*x^2*y + 2x - 4

        df = MPoly.gradient(f)

        @test string(df[1]) == "6.0xy+2.0"
        @test string(df[2]) == "3.0x^2"
    end
end;
