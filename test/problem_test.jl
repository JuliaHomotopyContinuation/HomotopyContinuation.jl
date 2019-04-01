@testset "Problems" begin
    F = equations(katsura(5))
    G = equations(cyclic(6))

    P1 = StartTargetInput(G, F, [rand(ComplexF64, 6), rand(ComplexF64, 6)])
    (PP1, start1) = problem_startsolutions(P1)
    @test PP1 isa ProjectiveProblem
    @test length(start1) == 2

    P2 = TotalDegreeInput(G)
    (PP2, start2) = problem_startsolutions(P2)
    @test PP2 isa ProjectiveProblem
    @test length(start2) == 720

    @test homotopy(PP2) isa AbstractHomotopy
    @test homvars(PP2) == (7,)

    @polyvar x y z
    F = [x^2+y^2+2z^2, x+y+3z]
    P, _ = problem_startsolutions(TotalDegreeInput(F), homvar=y)
    @test homvars(P) == (2,)

    P, _ = problem_startsolutions(TotalDegreeInput(F))
    @test homvars(P) === nothing

    @test_throws ErrorException problem_startsolutions(StartTargetInput(
        [x^2+y^2+z^2, x^4+y^4+z^3], [x^3+z^3,y^3-z^3], [rand(ComplexF64, 3)]))

    P, _ = problem_startsolutions(StartTargetInput(
        [x^2+y^2+z^2, x^4+y^4+z^4], [x^3+z^3,y^3-z^3], [rand(ComplexF64, 3)]))
    @test homvars(P) === nothing

    _, starts = problem_startsolutions(HC.input(
        [x^2+y^2+z^2, x^4+y^4+z^4], [x^3+z^3,y^3-z^3], rand(ComplexF64, 3)))
    @test length(starts) == 1

    P, _ = problem_startsolutions(StartTargetInput(
        [x^2+y^2+z^2, x^4+y^4+z^4], [x^3+z^3,y^3-z^3], [rand(ComplexF64, 3)]), homvar=z)
    @test homvars(P) === (3,)


    F = FPSystem(homogenize(equations(cyclic(6))))
    P, startvals = problem_startsolutions(TotalDegreeInput(F))
    @test homvars(P) == nothing
    @test length(startvals) == 720

    F = FPSystem(homogenize(equations(cyclic(6))))
    P, startvals = problem_startsolutions(TotalDegreeInput(F), homvar=5)
    @test homvars(P) == (5,)
    @test length(startvals) == 720
end
