@testset "polysystem" begin
    @TP.polyvar x y z
    f = [3.0x^2 * y^2 + 4x^3 * y + 2x * y^2, x^2 + y^3]
    @test evaluate(f, [3.0, 5.0]) == [1365.0, 134.0]
    df = jacobian(f)
    @test df([2.0, 2.0]) == [152.0 96.0;4.0 12.0]


    @test is_homogenous(f) == false
    @test is_homogenous([x^2+y^2, x+y]) == true
    
    @test is_homogenous(homogenize(f, z)) == true

    F = [x^4, 2x + 2y^2, 3x * y^2]
    @test degrees(F) == [4, 2, 3]
    @test nequations(F) == 3
    @test nvars(F) == 2


    f = 3.0x^2 + 2.0x * y - (1.0 + 0im) * y^2
    F = [f, f]
    g = (-2.5 + 2im)x^2 + -3.0x * y + 4.0 * y^2
    G = [g, g]
    @test weyl_dot(F, G) == 2 * weyl_dot(f,g)
    @test weyl_dot(F, F) == 2 * weyl_dot(f,f)
    @test weyl_dot(F, F) ≈ weyl_norm(F)^2

    F =[3.0x^2 + 2.0x * y - 1.0 * y^2, (-2.5 + 2im)x^2 + -3.0x * y + 4.0 * y^3]
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