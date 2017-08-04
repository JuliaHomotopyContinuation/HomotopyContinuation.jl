using HomotopyContinuation
using Base.Test

@testset "polysystem" begin

    @testset "Constructor from Poly" begin
        x, y = MPoly.generators(Float64, :x, :y)
        f = MPoly.system(x^2 + 3 * x * y^2 + 3)

        @test typeof(f)<:PolySystem{Float64}
    end

    @testset "eval" begin
        x, y = MPoly.generators(Float64, :x, :y)
        f = MPoly.system(3x^2 * y^2 + 4x^3 * y + 2x * y^2)
        @test evaluate(f, [3.0, 5.0]) == [1365.0]
    end

    @testset "jacobian" begin
        x, y = MPoly.generators(Float64, :x, :y)
        f = MPoly.system(3x^2 * y^2 + 4x^3 * y + 2x * y^2)

        df = jacobian(f)
        @test df([2.0, 2.0]) == [152.0 96.0]
    end

    @testset "is_homogenous" begin
        x, y = MPoly.generators(Float64, :x, :y)

        @test is_homogenous(MPoly.system(3x^2 * y + 4x^3 + 2x * y^2)) == true
        @test is_homogenous(MPoly.system(3x^2 * y + 4x^3 + 2x * y)) == false
    end

    @testset "degrees, nequations, nvars" begin
        x, y = MPoly.generators(Float64, :x, :y)
        F = MPoly.system(x^4, 2x + 2y^2, 3x * y^2)

        @test degrees(F) == [4, 2, 3]
        @test nequations(F) == 3
        @test nvars(F) == 2
    end

    @testset "weyl dot" begin
        x, y = MPoly.generators(Complex128, :x, :y)

        f = 3.0x^2 + 2.0x * y - 1.0 * y^2
        g = (-2.5 + 2im)x^2 + -3.0x * y + 4.0 * y^2

        F = MPoly.system(f)
        G = MPoly.system(g)

        @test weyl_dot(F, G) == 3.0 * conj(-2.5 + 2im) + 2.0 * (-3.0) / 2 + (-1.0) * 4.0
        @test weyl_dot(F, F) == 9.0 + 4.0  / 2 + 1.0

        @test weyl_dot(F, F) ≈ weyl_norm(F)^2
    end

    @testset "total_degree" begin
        x, y = MPoly.generators(Complex128, :x, :y)

        F = MPoly.system(3.0x^2 + 2.0x * y - 1.0 * y^2,  (-2.5 + 2im)x^2 + -3.0x * y + 4.0 * y^3)

        G, solutions = total_degree(F)

        @test length(solutions) == 6
        @test degrees(G) == [2, 3]

        for sol in solutions
            @test norm(evaluate(G, sol)) ≈ 0.0 atol = 1e-12
        end

        @test any(x -> abs(1 - x) > 1e-8, norm.(solutions)) == true

        G, solutions = total_degree(F, unit_roots=true)
        @test all(x -> norm(1 - abs2.(x)) < 1e-8, solutions)
    end
end