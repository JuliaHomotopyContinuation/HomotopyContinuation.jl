@testset "Monodromy" begin
    function toric_ed(A)
        d, n = size(A)
        @var t[1:d] y[1:n] u[1:n]

        φ = map(j -> prod(i -> t[i]^A[i, j], 1:d), 1:n)
        Dφ = [differentiate(φ[i], t[j]) for i = 1:n, j = 1:d]

        System([φ + y - u; Dφ' * y], parameters = u)
    end

    @testset "monodromy_solve" begin
        F = toric_ed([3 2 1 0; 0 1 2 3])

        result = monodromy_solve(
            F,
            target_solutions_count = 21,
            max_loops_no_progress = 20,
            threading = false,
        )
        @test is_success(result)
        @test length(solutions(result)) == 21
        @test isempty(multiplicities(solutions(result)))
        @test isempty(sprint(show, result)) == false

        # test seed
        result2 = monodromy_solve(
            F,
            target_solutions_count = 21,
            max_loops_no_progress = 20,
            threading = false,
            seed = result.seed,
        )
        @test result2.statistics.tracked_loops[] == result.statistics.tracked_loops[]

        result = monodromy_solve(
            F,
            target_solutions_count = 21,
            max_loops_no_progress = 20,
            threading = true,
        )
        @test is_success(result)
        @test length(solutions(result)) == 21

        # test that timeout works
        result_timeout = monodromy_solve(F, target_solutions_count = 21, timeout = 1e-12)
        @test length(solutions(result_timeout)) < 21

        # test input of length > 1
        x₀, p₀ = find_start_pair(F)
        result = monodromy_solve(
            F,
            [x₀ for _ = 1:30],
            p₀,
            target_solutions_count = 21,
            max_loops_no_progress = 20,
        )
        @test length(solutions(result)) == 21

        # test with false input
        @test_logs (:warn, "None of the provided solutions is a valid start solution.") (
            result = monodromy_solve(F, [rand(6) for _ = 1:10], p₀)
        )
        @test result.returncode == :invalid_startvalue
        result =
            monodromy_solve(F.expressions, [x₀, rand(6)], p₀, parameters = F.parameters)
        @test length(solutions(result)) == 21

        # different distance function
        result = monodromy_solve(F, x₀, p₀, distance = (x, y) -> 0.0)
        @test length(solutions(result)) == 1

        # Test stop heuristic using too high target_solutions_count
        result = monodromy_solve(F, x₀, p₀, target_solutions_count = 25)
        @test is_heuristic_stop(result)
        # Test stop heuristic with no target solutions count
        result = monodromy_solve(F, x₀, p₀)
        @test is_heuristic_stop(result)

        roots_of_unity(s) = begin
            t = cis(π * 2 / 3)
            t² = t * t
            (vcat(t * s[1], t * s[2], s[3:end]), vcat(t² * s[1], t² * s[2], s[3:end]))
        end
        result = monodromy_solve(
            F,
            x₀,
            p₀,
            target_solutions_count = 21,
            max_loops_no_progress = 100,
            equivalence_classes = false,
            group_action = roots_of_unity,
        )
        @test length(solutions(result)) == 21

        # equivalence classes
        result = monodromy_solve(
            F,
            x₀,
            p₀,
            equivalence_classes = true,
            target_solutions_count = 7,
            max_loops_no_progress = 20,
            group_actions = roots_of_unity,
        )
        @test nresults(result) == 7
        # Test that equivalence classes are on by default if we supply a group action
        result = monodromy_solve(
            F,
            x₀,
            p₀,
            group_action = roots_of_unity,
            max_loops_no_progress = 20,
        )
        @test nsolutions(result) == 7

        # AbstractSystem as input
        result = monodromy_solve(
            ModelKitSystem(F),
            group_action = roots_of_unity,
            target_solutions_count = 7,
            max_loops_no_progress = 20,
        )
        @test nsolutions(result) == 7

        # test reuse strategies
        result = monodromy_solve(
            F,
            x₀,
            p₀,
            group_action = roots_of_unity,
            target_solutions_count = 7,
            reuse_loops = :all,
            max_loops_no_progress = 200,
        )
        @test nresults(result) == 7

        result = monodromy_solve(
            F,
            x₀,
            p₀,
            group_action = roots_of_unity,
            target_solutions_count = 7,
            reuse_loops = :random,
            max_loops_no_progress = 200,
        )
        @test nresults(result) == 7

        result = monodromy_solve(
            F,
            x₀,
            p₀,
            group_action = roots_of_unity,
            target_solutions_count = 7,
            reuse_loops = :none,
            max_loops_no_progress = 200,
        )
        @test nresults(result) == 7
    end

    @testset "Method of moments" begin
        f = moments3()
        S₃ = SymmetricGroup(3)
        relabeling(v) = map(p -> [v[1:3][p]..., v[4:6][p]..., v[7:9][p]...], S₃)

        R = monodromy_solve(
            f;
            group_action = relabeling,
            show_progress = false,
            max_loops_no_progress = 1,
        )
        @test nsolutions(R) ≤ 225

        for threading in [true, false]
            R = monodromy_solve(
                f;
                group_action = relabeling,
                show_progress = false,
                max_loops_no_progress = 20,
                target_solutions_count = 225,
                threading = false,
            )
            @test nsolutions(R) == 225
        end
    end

    @testset "Projective + Group Actions" begin
        @var a[1:2] x[1:2] s[1:2] z m[0:5]

        f0 = a[1] + a[2]
        f1 = a[1] * x[1] + a[2] * x[2]
        f2 = a[1] * (x[1]^2 + s[1] * z) + a[2] * (x[2]^2 + s[2] * z)
        f3 = a[1] * (x[1]^3 + 3 * s[1] * x[1] * z) + a[2] * (x[2]^3 + 3 * s[2] * x[2] * z)
        f4 =
            a[1] * (x[1]^4 + 6 * s[1] * x[1]^2 * z + 3 * s[1]^2 * z^2) +
            a[2] * (x[2]^4 + 6 * s[2] * x[2]^2 * z + 3 * s[2]^2 * z^2)
        f5 =
            a[1] * (x[1]^5 + 10 * s[1] * x[1]^3 * z + 15 * x[1] * s[1]^2 * z^2) +
            a[2] * (x[2]^5 + 10 * s[2] * x[2]^3 * z + 15 * x[2] * s[2]^2 * z^2)

        M2 = System(
            [f0, f1, f2, f3, f4, f5] - m .* z .^ (1:6),
            variables = [a; x; s; z],
            parameters = m,
        )

        relabeling = let S₂ = SymmetricGroup(2)
            v -> map(p -> [v[1:2][p]..., v[3:4][p]..., v[5:6][p]..., v[7]], S₂)
        end

        R = monodromy_solve(
            M2;
            group_action = relabeling,
            show_progress = false,
            threading = false,
            max_loops_no_progress = 10,
        )

        @test nsolutions(R) == 9
    end

    @testset "permutations" begin
        @var x[1:2] a b c
        p₁ = x - [2; 0]
        p₂ = x - [-2; 0]
        c₁ = p₁[1]^2 + p₁[2]^2 - 1
        c₂ = p₂[1]^2 + p₂[2]^2 - 1

        F = System([c₁ * c₂; a * x[1] + b * x[2] - c]; parameters = [a, b, c])
        S = monodromy_solve(F, [[1.0, 0.0]], [1, 1, 1], permutations = true)
        A = permutations(S)
        B = permutations(S, reduced = false)

        @test size(A) == (2, 2)
        @test A == [1 2; 2 1] || A == [2 1; 1 2]
        @test size(B, 1) == 2
        @test size(B, 2) > 2

        start_sols = solve(F; target_parameters = [1, 1, 1], start_system = :total_degree)
        S = monodromy_solve(F, solutions(start_sols), [1, 1, 1]; permutations = true)
        C = permutations(S)
        @test size(C, 1) == 4
    end
end
