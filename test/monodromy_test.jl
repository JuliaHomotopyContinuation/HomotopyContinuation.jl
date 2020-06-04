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
        result = monodromy_solve(F, [rand(6) for _ = 1:10], p₀)
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

        R = monodromy_solve(
            f;
            group_action = relabeling,
            show_progress = false,
            max_loops_no_progress = 20,
            target_solutions_count = 225,
        )
        @test nsolutions(R) == 225
    end
end
