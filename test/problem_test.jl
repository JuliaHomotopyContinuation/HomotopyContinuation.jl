@testset "Problems" begin
    F = equations(katsura(5))
    G = equations(cyclic(6))

    P1 = Input.StartTarget(G, F, [rand(ComplexF64, 6), rand(ComplexF64, 6)])
    (PP1, start1) = problem_startsolutions(P1)
    @test PP1 isa ProjectiveProblem
    @test length(start1) == 2

    P2 = Input.TotalDegree(G)
    (PP2, start2) = problem_startsolutions(P2)
    @test PP2 isa ProjectiveProblem
    @test length(start2) == 720

    @test homotopy(PP2) isa Homotopies.AbstractHomotopy
    @test homvars(PP2) == (7,)

    @polyvar x y z
    F = [x^2+y^2+2z^2, x+y+3z]
    P, _ = problem_startsolutions(Input.TotalDegree(F), homvar=y)
    @test homvars(P) == (2,)

    P, _ = problem_startsolutions(Input.TotalDegree(F))
    @test homvars(P) == nothing

    @test_throws ErrorException problem_startsolutions(Input.StartTarget(
        [x^2+y^2+z^2, x^4+y^4+z^3], [x^3+z^3,y^3-z^3], [rand(ComplexF64, 3)]))

    P, _ = problem_startsolutions(Input.StartTarget(
        [x^2+y^2+z^2, x^4+y^4+z^4], [x^3+z^3,y^3-z^3], [rand(ComplexF64, 3)]))
    @test homvars(P) == nothing

    P, _ = problem_startsolutions(Input.StartTarget(
        [x^2+y^2+z^2, x^4+y^4+z^4], [x^3+z^3,y^3-z^3], [rand(ComplexF64, 3)]), homvar=z)
    @test homvars(P) == (3,)


    F = Systems.FPSystem(Utilities.homogenize(equations(cyclic(6))))
    degrees = Input.check_homogenous_degrees(F)
    P, startvals = problem_startsolutions(Input.TotalDegree(F, degrees))
    @test homvars(P) == nothing
    @test length(startvals) == 720

    F = Systems.FPSystem(Utilities.homogenize(equations(cyclic(6))))
    degrees = Input.check_homogenous_degrees(F)
    P, startvals = problem_startsolutions(Input.TotalDegree(F, degrees), homvar=5)
    @test homvars(P) == (5,)
    @test length(startvals) == 720
end
