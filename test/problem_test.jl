@testset "Problems" begin
    F = equations(katsura(5))
    G = equations(cyclic(6))

    P1 = StartTargetInput(G, F)
    (PP1, start1) = problem_startsolutions(P1, [rand(ComplexF64, 6), rand(ComplexF64, 6)])
    @test PP1 isa Problem{AffineTracking}
    @test length(start1) == 2

    P2 = TargetSystemInput(G)
    (PP2, start2) = problem_startsolutions(P2)
    @test PP2 isa Problem{AffineTracking}
    @test length(start2) == 720

    @test homotopy(PP2) isa AbstractHomotopy
    @test homvars(PP2) == nothing

    @polyvar x y z
    F = [x^2+y^2+2z^2, x+y+3z]
    P, _ = problem_startsolutions(TargetSystemInput(F), nothing, homvar=y)
    @test homvars(P) == (2,)

    P, _ = problem_startsolutions(TargetSystemInput(F), nothing)
    @test homvars(P) === nothing

    @test_throws ArgumentError problem_startsolutions(StartTargetInput(
        [x^2+y^2+z^2, x^4+y^4+z^3], [x^3+z^3,y^3-z^3]), [rand(ComplexF64, 3)])

    P, _ = problem_startsolutions(StartTargetInput(
        [x^2+y^2+z^2, x^4+y^4+z^4], [x^3+z^3,y^3-z^3]), [rand(ComplexF64, 3)])
    @test P isa Problem{ProjectiveTracking}
    @test homvars(P) === nothing

    P, _ = problem_startsolutions(StartTargetInput(
        [x^2+y^2+z^2, x^4+y^4+z^2], [x^3+z^3,y^3-z^2]), [rand(ComplexF64, 3)])
    @test P isa Problem{AffineTracking}

    # not supported right now
    @test_throws MethodError problem_startsolutions(
        StartTargetInput(SPSystem([x^2+y^2+z^2, x^4+y^4+z^2]), SPSystem([x^3+z^3,y^3-z^2])),
        [rand(ComplexF64, 3)])

    P, starts = problem_startsolutions(HC.input_startsolutions(
        [x^2+y^2+1, x^4+y^4+1], [x^3+1,y^3-1], rand(ComplexF64, 3))...)
    @test P isa Problem{AffineTracking}
    @test length(starts) == 1

    P, _ = problem_startsolutions(StartTargetInput(
        [x^2+y^2+z^2, x^4+y^4+z^4], [x^3+z^3,y^3-z^3]), [rand(ComplexF64, 3)], homvar=z)
    @test homvars(P) === (3,)


    F = FPSystem(homogenize(equations(cyclic(6))))
    P, startvals = problem_startsolutions(TargetSystemInput(F))
    @test P isa Problem{ProjectiveTracking}
    @test homvars(P) == nothing

    F = homogenize(equations(cyclic(6)))
    P, _ = problem_startsolutions(TargetSystemInput(F))
    @test P isa Problem{ProjectiveTracking}

    F = FPSystem(homogenize(equations(cyclic(6))))
    P, startvals = problem_startsolutions(TargetSystemInput(F), homvar=5)
    @test P isa Problem{ProjectiveTracking}
    @test homvars(P) == (5,)
    @test length(startvals) == 720
end
