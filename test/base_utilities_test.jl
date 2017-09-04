@testset "utilities" begin
    @testset "affine/projective" begin
        @test projective([1, 2, 3]) == [1, 1, 2, 3]
        @test affine(projective([1, 2, 3])) == [1, 2, 3]
    end


    @testset "totaldegree" begin
        @PolyImpl.polyvar x y

        F = PolySystem([
            3.0x^2 + 2.0x * y - 1.0 * y^2,
            (-2.5 + 2im)x^2 + -3.0x * y + 4.0 * y^3])

        G, solutions = totaldegree(F)

        @test length(solutions) == 6
        @test degrees(G) == [2, 3]

        for sol in solutions
            @test norm(evaluate(G, sol)) ≈ 0.0 atol = 1e-12
        end

        @test any(x -> abs(1 - x) > 1e-8, norm.(solutions)) == true

        G, solutions = totaldegree(F, unit_roots=true)
        @test all(x -> norm(1 - abs2.(x)) < 1e-8, solutions)

        @test_throws ErrorException totaldegree(PolySystem([x-1, x-2]))
    end

    @testset "randomsystem" begin
        F = randomsystem(Complex64, 2, 3, mindegree=2, maxdegree=5)

        @test nvariables(F) == 3
        @test length(F) == 2
        @test F isa FixedPolySystem.PolySystem{Complex64}

        @test all(d -> 2 ≤ d ≤ 5, degrees(F))

        @test randomsystem(2, 3) isa FixedPolySystem.PolySystem{Complex128}

        G = randomsystem([2, 3, 4], [:x, :y])
        @test nvariables(G) == 2
        @test length(G) == 3
        @test degrees(G) == [2, 3, 4]
    end
end
