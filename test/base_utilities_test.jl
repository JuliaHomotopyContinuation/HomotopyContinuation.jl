@testset "utilities" begin
    @testset "affine/projective" begin
        @test projective([1, 2, 3]) == [1, 1, 2, 3]
        @test affine(projective([1, 2, 3])) == [1, 2, 3]
    end


    @testset "total_degree" begin
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
    end
end
