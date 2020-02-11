@testset "Monodromy" begin
    function toric_ed(A)
        d, n = size(A)
        @polyvar t[1:d] y[1:n] u[1:n]

        φ = map(j -> prod(i -> t[i]^A[i, j], 1:d), 1:n)
        Dφ = [differentiate(φ[i], t[j]) for i = 1:n, j = 1:d]

        F = [φ + y - u; Dφ' * y]

        # We sample a random t, find a normal vector on this point and then assemble
        # an u where we know the solution to.
        t_rand = randn(Float64, size(A, 1))
        null = nullspace(map(fi -> fi(t => t_rand), Dφ'))
        y_rand = null * randn(Float64, size(null, 2))
        u₀ = map(fi -> fi(t => t_rand), φ) + y_rand

        x₀ = complex.([t_rand; y_rand])

        F, u, u₀, x₀
    end

    @testset "monodromy_solve" begin
        F, p, p₀, x₀ = toric_ed([3 2 1 0; 0 1 2 3])

        result = monodromy_solve(
            F,
            x₀,
            p₀,
            parameters = p,
            target_solutions_count = 21,
            max_loops_no_progress = 20,
            threading = false,
        )
        @test is_success(result)
        @test length(solutions(result)) == 21
        @test length(solutions(result, only_real = true)) >= 1
        @test result.statistics.ntrackedpaths ≥ 21
        @test result.statistics.nparametergenerations ≥ 1
        @test length(HC.UniquePoints(result.solutions).points) == 21
        @test isempty(sprint(show, result)) == false

        # test seed
        result2 = monodromy_solve(
            F,
            x₀,
            p₀,
            parameters = p,
            target_solutions_count = 21,
            max_loops_no_progress = 20,
            threading = false,
            seed = result.seed,
        )
        @test result2.statistics.ntrackedpaths == result.statistics.ntrackedpaths

        result = monodromy_solve(
            F,
            x₀,
            p₀,
            parameters = p,
            target_solutions_count = 21,
            max_loops_no_progress = 20,
            threading = true,
        )
        @test is_success(result)
        @test length(solutions(result)) == 21

        # test that timeout works
        result_timeout = monodromy_solve(
            F,
            x₀,
            p₀,
            parameters = p,
            target_solutions_count = 21,
            timeout = 1e-12,
        )
        @test length(result_timeout.solutions) < 21

        # test input of length > 1
        result = monodromy_solve(F, [x₀ for _ = 1:30], p₀, parameters = p)
        @test length(solutions(result)) == 21

        # test with false input
        result = monodromy_solve(F, [rand(6) for _ = 1:10], p₀, parameters = p)
        @test result.returncode == :invalid_startvalue
        result = monodromy_solve(F, [x₀, rand(6)], p₀, parameters = p)
        @test length(solutions(result)) == 21

        # different distance function
        result = monodromy_solve(F, x₀, p₀, parameters = p, distance = (x, y) -> 0.0)
        @test length(solutions(result)) == 1

        # Test stop heuristic using too high target_solutions_count
        result = monodromy_solve(F, x₀, p₀, parameters = p, target_solutions_count = 25)
        @test is_heuristic_stop(result)
        # Test stop heuristic with no target solutions count
        result = monodromy_solve(F, x₀, p₀, parameters = p)
        @test is_heuristic_stop(result)
        # Test stop heuristic with no target solutions count
        result = monodromy_solve(
            F,
            x₀,
            p₀,
            parameters = p,
            strategy = Triangle(useweights = true),
        )
        @test is_heuristic_stop(result)


        # By group_actions=nothing we force that complex conjugation is not used.
        result2 = monodromy_solve(
            F,
            x₀,
            p₀,
            parameters = p,
            target_solutions_count = 21,
            complex_conjugation = false,
            max_loops_no_progress = 100,
        )
        @test is_success(result2)

        result = monodromy_solve(
            F,
            x₀,
            p₀,
            parameters = p,
            target_solutions_count = 21,
            done_callback = ((_, _) -> true),
        )
        @test length(result.solutions) == 2

        roots_of_unity(s) = begin
            t = cis(π * 2 / 3)
            t² = t * t
            (vcat(t * s[1], t * s[2], s[3:end]), vcat(t² * s[1], t² * s[2], s[3:end]))
        end

        result = monodromy_solve(
            F,
            x₀,
            p₀,
            parameters = p,
            target_solutions_count = 21,
            max_loops_no_progress = 100,
            equivalence_classes = false,
            group_action = roots_of_unity,
        )
        @test length(result.solutions) == 21

        result = monodromy_solve(
            F,
            x₀,
            p₀,
            parameters = p,
            target_solutions_count = 21,
            max_loops_no_progress = 100,
            complex_conjugation = false, # disable complex conjugation to test it as a group action.
            equivalence_classes = false,
            group_actions = (roots_of_unity, s -> (conj.(s),)),
        )
        @test length(result.solutions) == 21
        @test length(solutions(result)) == 21
        @test length(real_solutions(result)) < 21
        test_treeviews(result)

        # group_actions as a vector
        result = monodromy_solve(
            F,
            x₀,
            p₀,
            parameters = p,
            target_solutions_count = 21,
            max_loops_no_progress = 100,
            complex_conjugation = false,
            equivalence_classes = false,
            group_actions = [roots_of_unity, s -> (conj.(s),)],
        )
        @test length(result.solutions) == 21


        # equivalence classes
        result = monodromy_solve(
            F,
            x₀,
            p₀,
            parameters = p,
            equivalence_classes = true,
            target_solutions_count = 7,
            max_loops_no_progress = 100,
            group_actions = roots_of_unity,
        )
        @test length(result.solutions) == 7
        # Test that equivalence classes are on by default if we supply a group action
        result = monodromy_solve(
            F,
            x₀,
            p₀,
            parameters = p,
            group_action = roots_of_unity,
            max_loops_no_progress = 20,
        )
        @test length(result.solutions) == 7

        # Test affine tracking
        result = monodromy_solve(
            F,
            x₀,
            p₀,
            parameters = p,
            affine_tracking = true,
            group_action = roots_of_unity,
            target_solutions_count = 7,
            max_loops_no_progress = 200,
        )
        @test length(result.solutions) == 7

        # AbstractSystem as input
        F_p = SPSystem(F; parameters = p)
        result = monodromy_solve(
            F_p,
            x₀,
            p₀,
            affine_tracking = true,
            group_action = roots_of_unity,
            target_solutions_count = 7,
            max_loops_no_progress = 200,
        )
        @test length(result.solutions) == 7


        # test reuse strategies
        result = monodromy_solve(
            F_p,
            x₀,
            p₀,
            group_action = roots_of_unity,
            target_solutions_count = 7,
            reuse_loops = :all,
            max_loops_no_progress = 200,
        )
        @test length(result.solutions) == 7
        result = monodromy_solve(
            F_p,
            x₀,
            p₀,
            group_action = roots_of_unity,
            target_solutions_count = 7,
            reuse_loops = :random,
            max_loops_no_progress = 200,
        )
        @test length(result.solutions) == 7
        result = monodromy_solve(
            F_p,
            x₀,
            p₀,
            group_action = roots_of_unity,
            target_solutions_count = 7,
            reuse_loops = :none,
            max_loops_no_progress = 200,
        )
        @test length(result.solutions) == 7
    end

    @testset "Method of Moments" begin
        Random.seed!(130793)
        @polyvar a[1:3] x[1:3] s[1:3]

        f0 = a[1] + a[2] + a[3]
        f1 = a[1] * x[1] + a[2] * x[2] + a[3] * x[3]
        f2 = a[1] * (x[1]^2 + s[1]) + a[2] * (x[2]^2 + s[2]) + a[3] * (x[3]^2 + s[3])
        f3 = a[1] * (x[1]^3 + 3 * s[1] * x[1]) +
            a[2] * (x[2]^3 + 3 * s[2] * x[2]) +
            a[3] * (x[3]^3 + 3 * s[3] * x[3])
        f4 = a[1] * (x[1]^4 + 6 * s[1] * x[1]^2 + 3 * s[1]^2) +
            a[2] * (x[2]^4 + 6 * s[2] * x[2]^2 + 3 * s[2]^2) +
            a[3] * (x[3]^4 + 6 * s[3] * x[3]^2 + 3 * s[3]^2)
        f5 = a[1] * (x[1]^5 + 10 * s[1] * x[1]^3 + 15 * x[1] * s[1]^2) +
            a[2] * (x[2]^5 + 10 * s[2] * x[2]^3 + 15 * x[2] * s[2]^2) +
            a[3] * (x[3]^5 + 10 * s[3] * x[3]^3 + 15 * x[3] * s[3]^2)
        f6 = a[1] * (x[1]^6 + 15 * s[1] * x[1]^4 + 45 * x[1]^2 * s[1]^2 + 15 * s[1]^3) +
            a[2] * (x[2]^6 + 15 * s[2] * x[2]^4 + 45 * x[2]^2 * s[2]^2 + 15 * s[2]^3) +
            a[3] * (x[3]^6 + 15 * s[3] * x[3]^4 + 45 * x[3]^2 * s[3]^2 + 15 * s[3]^3)
        f7 = a[1] * (
            x[1]^7 + 21 * s[1] * x[1]^5 + 105 * x[1]^3 * s[1]^2 + 105 * x[1] * s[1]^3
        ) +
            a[2] * (
            x[2]^7 + 21 * s[2] * x[2]^5 + 105 * x[2]^3 * s[2]^2 + 105 * x[2] * s[2]^3
        ) +
            a[3] * (
            x[3]^7 + 21 * s[3] * x[3]^5 + 105 * x[3]^3 * s[3]^2 + 105 * x[3] * s[3]^3
        )
        f8 = a[1] * (
            x[1]^8 +
                28 * s[1] * x[1]^6 +
                210 * x[1]^4 * s[1]^2 +
                420 * x[1]^2 * s[1]^3 +
                105 * s[1]^4
        ) +
            a[2] * (
            x[2]^8 +
                28 * s[2] * x[2]^6 +
                210 * x[2]^4 * s[2]^2 +
                420 * x[2]^2 * s[2]^3 +
                105 * s[2]^4
        ) +
            a[3] * (
            x[3]^8 +
                28 * s[3] * x[3]^6 +
                210 * x[3]^4 * s[3]^2 +
                420 * x[3]^2 * s[3]^3 +
                105 * s[3]^4
        )
        f = [f0, f1, f2, f3, f4, f5, f6, f7, f8]

        S₃ = SymmetricGroup(3)
        relabeling(v) = map(p -> (v[1:3][p]..., v[4:6][p]..., v[7:9][p]...), S₃)

        @polyvar p[1:9]

        function sample_moments(_ = nothing)
            a₀ = rand(3)
            a₀ ./= sum(a₀)
            x₀ = randn(3)
            s₀ = rand(3) .* 2

            y₀ = [a₀; x₀; s₀]
            p₀ = [fᵢ(a => a₀, x => x₀, s => s₀) for fᵢ in f]
            y₀, p₀
        end

        y₀, p₀ = sample_moments()

        R = monodromy_solve(
            f - p,
            y₀,
            p₀;
            parameters = p,
            group_action = relabeling,
            parameter_sampler = last ∘ sample_moments,
            show_progress = false,
            max_loops_no_progress = 5,
        )
        @test length(solutions(R)) ≤ 225

        R = monodromy_solve(
            f - p,
            y₀,
            p₀;
            parameters = p,
            group_action = relabeling,
            target_solutions_count = 225,
            affine_tracking = true,
            parameter_sampler = last ∘ sample_moments,
            show_progress = false,
            max_loops_no_progress = 20,
        )
        @test length(solutions(R)) == 225
    end

    @testset "Trace Test" begin

        @testset "Simple" begin
            @polyvar x y a b c
            f = x^2 + y^2 - 1
            l = a * x + b * y + c
            res = monodromy_solve(
                [f, l],
                [-0.6 - 0.8im, -1.2 + 0.4im],
                [1, 2, 3];
                parameters = [a, b, c],
            )
            @test verify_solution_completeness([f, l], res; parameters = [a, b, c])
        end

        @testset "ED Example" begin
            Random.seed!(1222)
            # define f
            @polyvar x y
            f = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
            # define new variables u₁, u₂ and λ₁
            @polyvar u[1:2] λ[1:1]
            # define the jacobian of F
            J = differentiate([f], [x, y])
            # J' defines the transpose of J
            C = [[x, y] - u - J' * λ; f]

            # create a random linear space
            l = sum(randn(ComplexF64, 3) .* [x, y, 1])
            # compute intersection points of the V(f) and V(l) and take the first solution
            s₀ = first(solutions(solve([f, l])))
            # The normal space is now generated by
            J_s₀ = [p([x, y] => s₀) for p in J]
            # create normal vector
            λ₀ = cis(2π * rand()) / norm(J_s₀)
            # x₀ - u₀ = J_s₀' * λ₀ => u = x -  J_s₀' * λ₀
            u₀ = s₀ - transpose(J_s₀) * λ₀ |> vec

            # This is not as realiable as I would like it to be, so just make a couple runs...
            is_complete = false
            for i = 1:5
                res = monodromy_solve(C, [s₀; λ₀], u₀; parameters = u)
                if verify_solution_completeness(C, res; parameters = u)
                    is_complete = true
                    break
                end
            end
            @test is_complete
        end
    end


    @testset "monodromy - automatic start pair" begin
        @polyvar x[1:4]
        # declare the variables for the 3 points
        @polyvar p[1:3, 1:3]
        p₁, p₂, p₃ = p[:, 1], p[:, 2], p[:, 3]
        quartic_monomials = monomials(sum(x)^4)
        @polyvar c[1:35]
        f = sum(c .* quartic_monomials)
        ∇f = differentiate(f, x)
        @polyvar u v
        F = [
            subs(f, x => [p₁; 1])
            subs(f, x => [p₂; 1])
            subs(f, x => [p₃; 1])
            (subs(∇f, x => [p₁; 1]) - u * subs(∇f, x => [p₂; 1]))
            (subs(∇f, x => [p₁; 1]) - v * subs(∇f, x => [p₃; 1]))
        ]

        x₀, c₀ = HC.find_start_pair(F, c)
        @test norm([f([p₁; p₂; p₃; u; v] => x₀, c => c₀) for f in F]) < 1e-8

        F1 = [subs(f, c[1] => 0.232) for f in F]
        x₀, c₀ = HC.find_start_pair(F1, c[2:end])
        @test norm([f([p₁; p₂; p₃; u; v] => x₀, c[2:end] => c₀) for f in F1]) < 1e-8

        # only use a small target_solutions_count, this is just to test that we actually
        # compute something
        res = monodromy_solve(F; parameters = c, target_solutions_count = 50)
        @test is_success(res)

        # Example where we cannot compute a start pair naively.
        @polyvar x y
        f = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        # define new variables u₁, u₂ and λ₁
        @polyvar u[1:2] λ[1:1]
        # define the jacobian of F
        J = differentiate([f], [x, y])
        # J' defines the transpose of J
        C = [[x, y] - u - J' * λ; f]

        @test_throws ArgumentError monodromy_solve(C; parameters = u)
        @test_throws ArgumentError HC.find_start_pair(C, u)
    end

    @testset "projective" begin
        @polyvar x y z a b c
        f = x^2 + y^2 - z^2
        l = a * x + b * y + c * z

        v = PVector([-0.6 - 0.8im, -1.2 + 0.4im, 1])
        res = monodromy_solve(
            [f, l],
            v,
            [1, 2, 3];
            parameters = [a, b, c], max_loops_no_progress = 20,
        )
        @test is_heuristic_stop(res)
        @test nsolutions(res) == 2
        @test isempty(real_solutions(res))

        @polyvar x y u v a b
        f = [x * y - a * u * v, x^2 - b * u^2]
        res = monodromy_solve(
            f,
            PVector([√2, 1], [√2, 1]),
            [2, 2];
            variable_groups = [[x, u], [y, v]],
            parameters = [a, b],
            max_loops_no_progress = 20,
        )
        @test is_heuristic_stop(res)
        @test nsolutions(res) == 2
        @test nreal(res) == 2
        @test sum(fubini_study.(real_solutions(res), solutions(res))) ≈ 0.0 atol = 1e-6
    end

    @testset "monodromy group" begin
        @polyvar x[1:2]  a b c
        p₁ = x - [2;0]
        p₂ = x - [-2;0]
        c₁ = p₁[1]^2 + p₁[2]^2 - 1
        c₂ = p₂[1]^2 + p₂[2]^2 - 1

        F = [c₁ * c₂; a * x[1] + b * x[2] - c]
        S = monodromy_solve(F, [[1.0, 0.0]], [1, 1, 1], parameters = [a;b;c])
        A = permutations(S)
        B = permutations(S, reduced = false)

        @test size(A) == (2,2)
        @test A == [1 2; 2 1] || A == [2 1; 1 2]
        @test size(B, 1) == 2
        @test size(B, 2) > 2

        F₀ = [c₁ * c₂; x[1] + x[2] - 1]
        start_sols=solutions(solve(F₀))
        S = monodromy_solve(F, start_sols, [1, 1, 1], parameters = [a;b;c])
        C = permutations(S)

        @test size(C, 1) == 4
    end
end
