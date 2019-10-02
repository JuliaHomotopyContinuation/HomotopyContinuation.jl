@testset "Problems" begin
    @testset "General Problem" begin
        F = equations(katsura(5))
        G = equations(cyclic(6))

        P1 = StartTargetInput(G, F)
        (PP1, start1) = problem_startsolutions(
            P1,
            [rand(ComplexF64, 6), rand(ComplexF64, 6)],
        )
        @test PP1 isa Problem{AffineTracking}
        @test length(start1) == 2

        P2 = TargetSystemInput(G)
        (PP2, start2) = problem_startsolutions(P2)
        @test PP2 isa Problem{AffineTracking}
        @test length(start2) == 720

        @test homotopy(PP2) isa AbstractHomotopy
        @test homvars(PP2) == nothing

        @polyvar x y z
        F = [x^2 + y^2 + 2 * z^2, x + y + 3z]
        P, _ = problem_startsolutions(TargetSystemInput(F), nothing, homvar = y)
        @test homvars(P) == (2,)

        P, _ = problem_startsolutions(TargetSystemInput(F), nothing)
        @test homvars(P) === nothing

        @test_throws ArgumentError problem_startsolutions(
            StartTargetInput([x^2 + y^2 + z^2, x^4 + y^4 + z^3], [x^3 + z^3, y^3 - z^3]),
            [rand(ComplexF64, 3)],
        )

        P, _ = problem_startsolutions(
            StartTargetInput([x^2 + y^2 + z^2, x^4 + y^4 + z^4], [x^3 + z^3, y^3 - z^3]),
            [rand(ComplexF64, 3)],
        )
        @test P isa Problem{ProjectiveTracking}
        @test homvars(P) === nothing

        P, _ = problem_startsolutions(
            StartTargetInput([x^2 + y^2 + z^2, x^4 + y^4 + z^2], [x^3 + z^3, y^3 - z^2]),
            [rand(ComplexF64, 3)],
        )
        @test P isa Problem{AffineTracking}

        # not supported right now
        @test_throws MethodError problem_startsolutions(
            StartTargetInput(
                SPSystem([x^2 + y^2 + z^2, x^4 + y^4 + z^2]),
                SPSystem([x^3 + z^3, y^3 - z^2]),
            ),
            [rand(ComplexF64, 3)],
        )

        P, starts = problem_startsolutions(HC.input_startsolutions(
            [x^2 + y^2 + 1, x^4 + y^4 + 1],
            [x^3 + 1, y^3 - 1],
            rand(ComplexF64, 3),
        )...)
        @test P isa Problem{AffineTracking}
        @test length(starts) == 1

        P, _ = problem_startsolutions(
            StartTargetInput([x^2 + y^2 + z^2, x^4 + y^4 + z^4], [x^3 + z^3, y^3 - z^3]),
            [rand(ComplexF64, 3)],
            homvar = z,
        )
        @test homvars(P) === (3,)


        F = FPSystem(homogenize(equations(cyclic(6))))
        P, startvals = problem_startsolutions(TargetSystemInput(F))
        @test P isa Problem{ProjectiveTracking}
        @test homvars(P) == nothing

        F = homogenize(equations(cyclic(6)))
        P, _ = problem_startsolutions(TargetSystemInput(F))
        @test P isa Problem{ProjectiveTracking}

        F = FPSystem(homogenize(equations(cyclic(6))))
        P, startvals = problem_startsolutions(TargetSystemInput(F), homvar = 5)
        @test P isa Problem{ProjectiveTracking}
        @test homvars(P) == (5,)
        @test length(startvals) == 720
    end

    @testset "PolyhedralProblem" begin
        @polyvar x y z
        F = [x^2 + y^2 + 2 * z^2, x + y + 3z]
        P, starts = problem_startsolutions(F; start_system = :polyhedral)
        @test P isa HC.PolyhedralProblem
        @test starts isa HC.PolyhedralStartSolutionsIterator

        @polyvar u v w
        F = [x^2 + y^2 + 2 * z^2, x + y + 3z] ∘ [u, v, w]
        P, starts = problem_startsolutions(F; start_system = :polyhedral)
        @test P isa HC.PolyhedralProblem
        @test starts isa HC.PolyhedralStartSolutionsIterator
    end

    @testset "multi-homogenous" begin
        @polyvar z[1:6, 1:3]
        F = let
            p = [1, 1, 0]
            α = randn(5)
            a = randn(9)
            ## define the system of polynomials
            f = [z[i, :] ⋅ z[i, :] for i = 2:5]
            g = [z[i, :] ⋅ z[i+1, :] for i = 1:5]
            h = sum(a[i] .* (z[i, :] × z[i+1, :]) for i = 1:3) +
                sum(a[i+4] .* z[i, :] for i = 2:5)
            F′ = [f .- 1; g .- cos.(α); h .- p]
            ## assign values to z₁ and z₆
            [subs(f, z[1, :] => [1, 0, 0], z[6, :] => [1, 0, 0]) for f in F′]
        end
        group1 = [[z[2, :]; z[4, :]], [z[3, :]; z[5, :]]]
        group2 = [z[2, :], z[4, :], z[3, :], z[5, :]]

        @test bezout_number(F) == 1024
        @test bezout_number(F; variable_groups = group1) == 320
        @test bezout_number(F; variable_groups = group2) == 576

        _, starts = problem_startsolutions(F; variable_groups = group1)
        @test length(starts) == 320
        @test length(collect(starts)) == 320

        _, starts = problem_startsolutions(F; variable_groups = group2)
        @test length(starts) == 576
        @test length(collect(starts)) == 576
    end

    @testset "Overdetermined" begin
        @polyvar x y z
        prob, starts = problem_startsolutions([
            x - 2,
            y^2 + 3 * z^2,
            z^3 + x^3,
            z + x^2 + 3,
        ])
        @test prob isa HC.OverdeterminedProblem
        @test starts isa HC.TotalDegreeSolutionIterator
        @test starts.degrees == [3, 2, 2]

        # overdetermined, homogenous
        prob, starts = problem_startsolutions([x - 2y, y^2 + 3 * x^2])
        @test prob isa HC.OverdeterminedProblem
        @test starts isa HC.TotalDegreeSolutionIterator
        @test starts.degrees == [2]

        prob, starts = problem_startsolutions([x - 2y, y^2 + 3 * x^2, x^3 + y^3])
        @test prob isa HC.OverdeterminedProblem
        @test starts isa HC.TotalDegreeSolutionIterator
        @test starts.degrees == [3]

        # Overdetermined, polyhedral
        prob, starts = problem_startsolutions(
            [x - 2, y^2 + 3 * z^2, z^3 + x^3, z + x^2 + 3];
            start_system = :polyhedral,
        )
        @test prob isa HC.OverdeterminedProblem
        @test prob.problem isa HC.PolyhedralProblem
        @test starts isa HC.PolyhedralStartSolutionsIterator

        # polyhedral needs to expand Composition
        @polyvar u v w
        prob, starts = problem_startsolutions(
            [x - 2, y^2 + 3 * z^2, z^3 + x^3, z + x^2 + 3] ∘ [u, v, w];
            start_system = :polyhedral,
        )
        @test prob isa HC.OverdeterminedProblem
        @test prob.problem isa HC.PolyhedralProblem
        @test starts isa HC.PolyhedralStartSolutionsIterator

        # Abstract Systems
        prob, starts = problem_startsolutions(FPSystem([
            x - 2z,
            y^2 + 3 * z^2,
            z^3 + x^3,
            z + x,
        ]))
        @test prob isa HC.OverdeterminedProblem{HC.ProjectiveTracking}
        @test starts isa HC.TotalDegreeSolutionIterator
        @test starts.degrees == [1, 2, 3]
    end

    @testset "Invalid input" begin
        @polyvar x y z
        @test_throws ArgumentError problem_startsolutions([x - 2y + 2, 0])
        @test_throws ArgumentError problem_startsolutions(
            [x - 2, y - 2],
            [x - 2, y - 2, y + 2],
            [[2, -3]],
        )
        @test_throws ArgumentError problem_startsolutions(
            [x - 2z, y^2 + 3z, z^3 + x^3],
            homvar = z,
        )

        # constant term
        @test_throws ArgumentError problem_startsolutions([subs(x + 2, x => 2), y^2 + 3x])
        @polyvar u v
        @test_throws ArgumentError problem_startsolutions([subs(x + 2, x => 2), y^2 + 3x] ∘
                                                          [u, v])
        @test_throws ArgumentError problem_startsolutions([x - 2z, y^2 + 3z], homvar = z)
        # test numerical homogeneous check fails
        @polyvar x y z
        G = FPSystem([x - 2z, y^2 + 3z])
        @test_throws ArgumentError problem_startsolutions(G, homvar = 3)
    end
end
