@testset "Problems" begin
    F = equations(katsura(5))
    G = equations(cyclic(6))

    P1 = Input.StartTarget(G, F, [rand(ComplexF64, 6), rand(ComplexF64, 6)])
    (PP1, start1) = Problems.problem_startsolutions(P1)
    @test PP1 isa Problems.Projective
    @test length(start1) == 2

    P2 = Input.TotalDegree(G)
    (PP2, start2) = Problems.problem_startsolutions(P2)
    @test PP2 isa Problems.Projective
    @test length(start2) == 720

    @test Problems.homotopy(PP2) isa Homotopies.AbstractHomotopy
    @test Problems.homvars(PP2) == (7,)

    @polyvar x y z
    F = [x^2+y^2+2z^2, x+y+3z]
    P, _ = Problems.problem_startsolutions(Input.TotalDegree(F), homvar=y)
    @test Problems.homvars(P) == (2,)

    P, _ = Problems.problem_startsolutions(Input.TotalDegree(F))
    @test Problems.homvars(P) == nothing

    @test_throws ErrorException Problems.problem_startsolutions(Input.StartTarget(
        [x^2+y^2+z^2, x^4+y^4+z^3], [x^3+z^3,y^3-z^3], [rand(ComplexF64, 3)]))

    P, _ = Problems.problem_startsolutions(Input.StartTarget(
        [x^2+y^2+z^2, x^4+y^4+z^4], [x^3+z^3,y^3-z^3], [rand(ComplexF64, 3)]))
    @test Problems.homvars(P) == nothing

    P, _ = Problems.problem_startsolutions(Input.StartTarget(
        [x^2+y^2+z^2, x^4+y^4+z^4], [x^3+z^3,y^3-z^3], [rand(ComplexF64, 3)]), homvar=z)
    @test Problems.homvars(P) == (3,)


    F = Systems.FPSystem(Utilities.homogenize(equations(cyclic(6))))
    degrees = Input.check_homogenous_degrees(F)
    P, startvals = Problems.problem_startsolutions(Input.TotalDegree(F, degrees))
    @test Problems.homvars(P) == nothing
    @test length(startvals) == 720

    F = Systems.FPSystem(Utilities.homogenize(equations(cyclic(6))))
    degrees = Input.check_homogenous_degrees(F)
    P, startvals = Problems.problem_startsolutions(Input.TotalDegree(F, degrees), homvar=5)
    @test Problems.homvars(P) == (5,)
    @test length(startvals) == 720
end
