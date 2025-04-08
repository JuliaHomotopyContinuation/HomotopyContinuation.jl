@testset "Monodromy" begin
    Random.seed!(0x8b868323)

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
            max_loops_no_progress = 50,
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
            max_loops_no_progress = 50,
            threading = false,
            seed = result.seed,
        )
        @test result2.statistics.tracked_loops[] == result.statistics.tracked_loops[]

        result = monodromy_solve(
            F,
            target_solutions_count = 21,
            max_loops_no_progress = 50,
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
            max_loops_no_progress = 50,
        )
        @test length(solutions(result)) == 21

        # # test with false input
        # @test_logs (:warn, "None of the provided solutions is a valid start solution.") (
        #     result = monodromy_solve(F, [rand(6) for _ = 1:10], p₀)
        # )
        # @test result.returncode == :invalid_startvalue

        result =
            monodromy_solve(F.expressions, [x₀, rand(6)], p₀, parameters = F.parameters)
        @test length(solutions(result)) == 21

        # distance function that satisfies triangle inequality
        result = monodromy_solve(F, x₀, p₀, distance = (x, y) -> 0.0)
        @test length(solutions(result)) == 1

        # distance function that does not satisfy triangle inequality
        result = monodromy_solve(F, x₀, p₀, distance = (x, y) -> norm(x - y, 2)^2)
        @test length(solutions(result)) == 21

        # don't use triangle inequality
        result = monodromy_solve(F, x₀, p₀, triangle_inequality = false)
        @test length(solutions(result)) == 21

        # use triangle inequality
        result = monodromy_solve(F, x₀, p₀, triangle_inequality = true)
        @test length(solutions(result)) == 21

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
            max_loops_no_progress = 50,
            group_actions = roots_of_unity,
        )
        @test nresults(result) == 7
        # Test that equivalence classes are on by default if we supply a group action
        result = monodromy_solve(
            F,
            x₀,
            p₀,
            group_action = roots_of_unity,
            max_loops_no_progress = 50,
        )
        @test nsolutions(result) == 7

        # AbstractSystem as input
        result = monodromy_solve(
            InterpretedSystem(F),
            group_action = roots_of_unity,
            target_solutions_count = 7,
            max_loops_no_progress = 50,
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
        S = monodromy_solve(
            F,
            [[1.0, 0.0]],
            [1, 1, 1],
            permutations = true,
            max_loops_no_progress = 20,
        )
        A = HC.permutations(S)
        B = HC.permutations(S, reduced = false)

        @test size(A) == (2, 2)
        @test A == [1 2; 2 1] || A == [2 1; 1 2]
        @test size(B, 1) == 2
        @test size(B, 2) > 2

        start_sols = solve(F; target_parameters = [1, 1, 1], start_system = :total_degree)
        S = monodromy_solve(F, solutions(start_sols), [1, 1, 1]; permutations = true)
        C = permutations(S)
        @test size(C, 1) == 4
    end

    @testset "Linear subspaces" begin
        @var x[1:4]
        f1 = rand_poly(x, 6)
        F = System([f1], x)
        res = monodromy_solve(F; dim = 3, threading = false)
        @test nsolutions(res) == 6
        @test trace(res) < 1e-6
        @test is_success(res)

        res = monodromy_solve(F; dim = 3, trace_test = false, threading = false)
        @test nsolutions(res) == 6
        @test is_heuristic_stop(res)

        res = monodromy_solve(
            F;
            dim = 3,
            trace_test = true,
            trace_test_tol = 1e-50,
            threading = false,
        )
        @test nsolutions(res) == 6
        @test is_heuristic_stop(res)

        @var x[1:4]
        f1 = rand_poly(x, 6; homogeneous = true)
        F = System([f1], x)
        res = monodromy_solve(F; dim = 2, compile = false)
        @test nsolutions(res) == 6
        @test trace(res) < 1e-6
        @test is_success(res)

        F = System([f1, rand_poly(x, 3), rand_poly(x, 4)], x)
        res = monodromy_solve(F; codim = 3, compile = false)
        @test nsolutions(res) == 72
        @test trace(res) < 1e-6
        @test is_success(res)
    end

    @testset "Monodromy rational functions" begin
        @var x y
        @var u[1:4]
        f1 = u[1] / x^2 + u[2]
        f2 = u[3] / y^2 + u[4]
        F = System([f1, f2], variables = [x, y], parameters = u[1:4])
        r1 = monodromy_solve(
            F,
            compile = true,
            target_solutions_count = 4,
            max_loops_no_progress = 100,
        )
        r2 = monodromy_solve(
            F,
            compile = false,
            target_solutions_count = 4,
            max_loops_no_progress = 100,
        )
        @test nsolutions(r1) == nsolutions(r2) == 4

        @var A1[1:3, 1:4] A2[1:3, 1:4]
        @var x[1:3] u1[1:2] u2[1:2]
        y1 = A1 * [x; 1]
        y2 = A2 * [x; 1]
        f = sum((u1 - y1[1:2] ./ y1[3]) .^ 2) + sum((u2 - y2[1:2] ./ y2[3]) .^ 2)
        F = System(
            differentiate(f, x),
            variables = x,
            parameters = [u1; u2; vec(A1); vec(A2)],
        )
        r = monodromy_solve(
            F,
            compile = false,
            max_loops_no_progress = 100,
            target_solutions_count = 6,
        )
        @test nsolutions(r) == 6
    end

    @testset "parameter homotopy with monodromy result" begin
        F = toric_ed([3 2 1 0; 0 1 2 3])
        mres = monodromy_solve(
            F,
            target_solutions_count = 21,
            max_loops_no_progress = 20,
            threading = false,
        )
        @test nsolutions(solve(F, mres, target_parameters = randn(ComplexF64, 4))) == 21
    end

    @testset "Verify solution_completeness" begin
        @var x y a b c
        f = x^2 + y^2 - 1
        l = a * x + b * y + c
        sys = System([f, l]; parameters = [a, b, c])
        res = solve(sys, [-0.6 - 0.8im, -1.2 + 0.4im]; target_parameters = [1, 2, 3])
        @test verify_solution_completeness(sys, solutions(res), [1, 2, 3])
        @test verify_solution_completeness(
            sys,
            solutions(res),
            [1, 2, 3],
            show_progress = false,
        )
        @test verify_solution_completeness(
            sys,
            solutions(res),
            [1, 2, 3],
            monodromy_options = (compile = false,),
            parameter_homotopy_options = (compile = false,),
            show_progress = false,
        )
        @test false == verify_solution_completeness(
            sys,
            solutions(res)[1:1],
            [1, 2, 3],
            monodromy_options = (compile = false,),
            parameter_homotopy_options = (compile = false,),
            show_progress = true,
        )
        @test false == verify_solution_completeness(
            sys,
            solutions(res),
            [1, 2, 3],
            monodromy_options = (compile = false,),
            parameter_homotopy_options = (compile = false,),
            show_progress = false,
            trace_tol = 1e-60,
        )
    end

    @testset "monodromy with predefined tolerance for unique_points" begin
        d = 3
        n = 4
        M = binomial(d - 1 + 2, 2)
        D = (n - 1) * M + 3
        N = binomial(n - 1 + d, d)
        # Example from
        # https://www.juliahomotopycontinuation.org/examples/symmetroids/

        @var x[0:n-1] a[1:D]
        A = map(0:n-2) do ℓ
            Aᵢ = []
            for i = 1:d
                k = ℓ * M + sum(d - j for j = 0:i-1)
                push!(Aᵢ, [zeros(i - 1); a[k-(d-i):k]])
            end
            Aᵢ = hcat(Aᵢ...)
            (Aᵢ + transpose(Aᵢ)) ./ 2
        end

        A₀ = [a[D-2] 0 0; 0 a[D-1] 0; 0 0 a[D]]
        μ = x[1] .* A₀ + sum(x[i+1] .* A[i] for i = 1:n-1)
        f = System(coefficients(det(μ), x))

        J₀ = jacobian(InterpretedSystem(f), randn(ComplexF64, D))
        dimQ = rank(J₀)
        dim_fibers = D - dimQ
        S1 = qr(randn(ComplexF64, D, D)).Q
        S = S1[:, 1:dimQ]
        s = S1[:, dimQ+1]

        @var b[1:dimQ]
        L₁ = System(S * b + s)
        @var k[1:N+1] f₀[1:length(f)]
        b₁ = randn(ComplexF64, dimQ)
        f₁ = f(L₁(b₁))

        R1 = qr(randn(ComplexF64, N, N)).Q
        R = [transpose(k[1:N]); R1[1:(dimQ-1), :]]
        r = [k[N+1]; R[2:end, :] * f₁]

        L₂ = System(R * f₀ - r, variables = f₀; parameters = k)


        p₁ = randn(ComplexF64, N)
        params = [p₁; transpose(p₁) * f₁]

        dist(x, y) = norm(f(L₁(x)) - f(L₁(y)), Inf)

        points = monodromy_solve(
            L₂ ∘ f ∘ L₁,
            [b₁],
            params,
            distance = dist,
            compile = false,
            unique_points_rtol = 1e-8,
            unique_points_atol = 1e-14,
            target_solutions_count = 305,
        )

        UP = unique_points(solutions(points), metric = dist, rtol = 1e-8, atol = 1e-14)

        @test length(solutions(points)) == length(UP)
    end
end
