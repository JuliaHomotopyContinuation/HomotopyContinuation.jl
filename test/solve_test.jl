@testset "solve" begin

    @testset "total degree (simple)" begin
        @var x y
        affine_sqr = System([
            2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3,
            2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5,
        ])
        @test count(is_success, track.(total_degree(affine_sqr; compile = false)...)) == 2

        @var x y z
        proj_square = System([
            2.3 * x^2 + 1.2 * y^2 + 3x * z - 2y * z + 3 * z^2,
            2.3 * x^2 + 1.2 * y^2 + 5x * z + 2y * z - 5 * z^2,
        ])
        @test count(is_success, track.(total_degree(proj_square; compile = false)...)) == 4

        @var x y
        affine_ov = System([
            (x^2 + y^2 + x * y - 3) * (x + 3),
            (x^2 + y^2 + x * y - 3) * (y - x + 2),
            2x + 5y - 3,
        ])
        @test count(is_success, track.(total_degree(affine_ov; compile = false)...)) == 2

        @var x y
        affine_ov_reordering = System([
            (x^2 + y^2 + x * y - 3) * (x + 3),
            2x + 5y - 3,
            (x^2 + y^2 + x * y - 3) * (y^2 - x + 2),
        ])
        tracker, starts = total_degree(affine_ov_reordering; compile = false)
        @test length(starts) == 4 * 3
        @test count(is_success, track.(tracker, starts)) == 2

        @var x y z
        proj_ov = System([
            (x^2 + y^2 + x * y - 3 * z^2) * (x + 3z),
            (x^2 + y^2 + x * y - 3 * z^2) * (y - x + 2z),
            2x + 5y - 3z,
        ])
        @test count(is_success, track.(total_degree(proj_ov; compile = false)...)) == 2

        @var x y
        proj_ov_reordering = System([
            (x^2 + y^2 + x * y - 3 * z^2) * (x + 3z),
            2x + 5y - 3z,
            (x^2 + y^2 + x * y - 3 * z^2) * (y^2 - x * z + 2 * z^2),
        ])
        tracker, starts = total_degree(proj_ov_reordering; compile = false)
        @test length(starts) == 4 * 3
        @test count(is_success, track.(tracker, starts)) == 2

        @var x y
        affine_underdetermined = System([2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3])
        @test_throws HC.FiniteException total_degree(affine_underdetermined)

        @var x y z
        proj_underdetermined = System([2.3 * x^2 + 1.2 * y^2 + 3x * z])
        @test_throws HC.FiniteException total_degree(proj_underdetermined)
    end

    @testset "total degree (variable groups)" begin
        @var x y v w
        affine_sqr = System([x * y - 2, x^2 - 4], variable_groups = [[x], [y]])
        tracker, starts = total_degree(affine_sqr; compile = false)
        @test length(collect(starts)) == 2
        @test count(is_success, track.(tracker, starts)) == 2

        @var x y v w
        proj_sqr =
            System([x * y - 2v * w, x^2 - 4 * v^2], variable_groups = [[x, v], [y, w]])
        tracker, starts = total_degree(proj_sqr)
        @test length(collect(starts)) == 2
        @test count(is_success, track.(tracker, starts)) == 2

        @var x y v w
        affine_ov = System(
            [(x^2 - 4) * (x * y - 2), x * y - 2, x^2 - 4],
            variable_groups = [[x], [y]],
        )
        tracker, starts = total_degree(affine_ov; compile = false)
        @test count(is_success, track.(tracker, starts)) == 2
        @var x y v w
        proj_ov = System(
            [(x^2 - 4 * v^2) * (x * y - v * w), x * y - v * w, x^2 - v^2],
            variable_groups = [[x, v], [y, w]],
        )
        tracker, starts = total_degree(proj_ov; compile = false)
        @test count(is_success, track.(tracker, starts)) == 2
    end

    @testset "polyhedral" begin
        @var x y
        affine_sqr = System([
            2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3,
            2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5,
        ])
        @test count(is_success, track.(polyhedral(affine_sqr; compile = false)...)) == 2

        @var x y z
        proj_square = System([
            2.3 * x^2 + 1.2 * y^2 + 3x * z - 2y * z + 3 * z^2,
            2.3 * x^2 + 1.2 * y^2 + 5x * z + 2y * z - 5 * z^2,
        ])
        @test count(is_success, track.(polyhedral(proj_square; compile = false)...)) == 4

        @var x y
        affine_ov = System([
            (x^2 + y^2 + x * y - 3) * (x + 3),
            (x^2 + y^2 + x * y - 3) * (y - x + 2),
            2x + 5y - 3,
        ])
        @test count(is_success, track.(polyhedral(affine_ov; compile = false)...)) == 2

        @var x y z
        proj_ov = System([
            (x^2 + y^2 + x * y - 3 * z^2) * (x + 3z),
            (x^2 + y^2 + x * y - 3 * z^2) * (y - x + 2z),
            2x + 5y - 3z,
        ])
        @test count(is_success, track.(polyhedral(proj_ov; compile = false)...)) == 2

        @var x y
        affine_underdetermined = System([2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3])
        @test_throws HC.FiniteException polyhedral(affine_underdetermined)

        @var x y z
        proj_underdetermined = System([2.3 * x^2 + 1.2 * y^2 + 3x * z])
        @test_throws HC.FiniteException polyhedral(proj_underdetermined)
    end

    @testset "overdetermined" begin
        @testset "3 by 5 minors" begin
            res = solve(
                minors();
                start_system = :total_degree,
                compile = false,
                show_progress = false,
            )
            @test count(is_success, res) == 80
            @test count(is_excess_solution, res) == 136
        end
    end

    @testset "Result" begin
        d = 2
        @var x y a[1:6]
        F = System(
            [
                (a[1] * x^d + a[2] * y) * (a[3] * x + a[4] * y) + 1,
                (a[1] * x^d + a[2] * y) * (a[5] * x + a[6] * y) + 1,
            ];
            parameters = a,
        )
        res = solve(F; target_parameters = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32])

        @test startswith(sprint(show, res), "Result with 3 solutions")
        @test seed(res) isa UInt32
        test_treeviews(res)

        seeded_res = solve(
            F;
            target_parameters = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32],
            seed = seed(res),
        )
        @test seed(seeded_res) == seed(res)
        test_treeviews(res)

        @test length(path_results(res)) == ntracked(res) == 7
        @test length(results(res)) == nresults(res) == 3
        @test length(solutions(res)) == 3
        @test length(findall(is_success, res)) == 3
        @test real_solutions(res) isa Vector{Vector{Float64}}
        @test length(real_solutions(res)) == nreal(res) == 1
        @test length(nonsingular(res)) == nnonsingular(res) == 3
        @test isempty(singular(res))
        @test nsingular(res) == 0
        @test length(at_infinity(res)) == nat_infinity(res) == 4
        @test isempty(failed(res))
        @test nfailed(res) == 0
        @test nexcess_solutions(res) == 0
        @test !isempty(sprint(show, statistics(res)))

        # singular result
        @var x y
        g = System([(29 / 16) * x^3 - 2 * x * y, x^2 - y])
        res = solve(g; start_system = :total_degree)
        @test startswith(sprint(show, res), "Result with 1 solution")
        @test seed(res) isa UInt32
        test_treeviews(res)
        @test !isempty(sprint(show, statistics(res)))
        @test !isempty(sprint(show, res))
    end

    @testset "composition" begin
        @var a b c x y z u v
        e = System([u + 1, v - 2])
        f = System([a * b - 2, a * c - 1])
        g = System([x + y, y + 3, x + 2])

        res = solve(e ∘ f ∘ g; start_system = :total_degree, compile = false)
        @test nsolutions(res) == 2

        res = solve(e ∘ f ∘ g; start_system = :polyhedral)
        @test nsolutions(res) == 2
    end

    @testset "paths to track" begin
        @var x y
        f = System([2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y])

        @test paths_to_track(f; start_system = :total_degree) == 16
        @test paths_to_track(f; start_system = :total_degree) == 16
        @test paths_to_track(f; start_system = :polyhedral) == 8
        @test paths_to_track(f; start_system = :polyhedral, only_non_zero = true) == 3
        @test paths_to_track(f) == 8
        @test_deprecated bezout_number(f) == 16
        @test mixed_volume(f) == 3

        @var x y a
        g = System([2y + a * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y], parameters = [a])
        @test paths_to_track(g; start_system = :total_degree) == 16
        @test paths_to_track(g; start_system = :polyhedral) == 8
    end

    @testset "solve (parameter homotopy)" begin
        # affine
        @var x a y b
        F = System([x^2 - a, x * y - a + b], [x, y], [a, b])
        s = [1, 1]
        res = solve(F, [s]; start_parameters = [1, 0], target_parameters = [2, 4])
        @test nsolutions(res) == 1
        res = solve(
            InterpretedSystem(F),
            [s];
            start_parameters = [1, 0],
            target_parameters = [2, 4],
            threading = false,
        )
        @test nsolutions(res) == 1

        @var x a y b
        F = System([x^2 - a], [x, y], [a, b])
        s = [1, 1]
        @test_throws FiniteException(1) solve(
            F,
            [s];
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )

        # proj
        @var x a y b z
        F_proj = System([x^2 - a * z^2, x * y + (b - a) * z^2], [x, y, z], [a, b])
        s = [1, 1, 1]
        res = solve(F_proj, [s]; start_parameters = [1, 0], target_parameters = [2, 4])
        @test nsolutions(res) == 1
        res = solve(InterpretedSystem(F_proj), [s]; p₁ = [1, 0], p₀ = [2, 4])
        @test nsolutions(res) == 1

        F_proj_err = System([x * y + (b - a) * z^2], [x, y, z], [a, b])
        @test_throws FiniteException solve(F_proj_err, [s]; p₁ = [1, 0], p₀ = [2, 4])

        # multi-proj
        @var x y v w a b
        F_multi_proj = System(
            [x * y - a * v * w, x^2 - b * v^2],
            parameters = [a, b],
            variable_groups = [[x, v], [y, w]],
        )
        S = [
            [
                -1.808683149843597 + 0.2761582523875564im,
                -0.9043415749217985 + 0.1380791261937782im,
                -0.0422893850686111 - 0.7152002569359284im,
                -0.0422893850686111 - 0.7152002569359283im,
            ],
            [
                -0.36370464807054353 + 0.6777414371333245im,
                0.18185232403527177 - 0.33887071856666223im,
                -0.3348980281838583 - 0.7759382220656511im,
                0.3348980281838583 + 0.7759382220656511im,
            ],
        ]
        res = solve(F_multi_proj, S; start_parameters = [2, 4], target_parameters = [3, 5])
        @test nsolutions(res) == 2
        res = solve(InterpretedSystem(F_multi_proj), S; p₁ = [2, 4], p₀ = [3, 5])
        @test nsolutions(res) == 2

        F_multi_proj_err = System(
            [x * y - a * v * w],
            parameters = [a, b],
            variable_groups = [[x, v], [y, w]],
        )
        @test_throws FiniteException(1) solve(F_multi_proj_err, S; p₁ = [2, 4], p₀ = [3, 5])
    end

    @testset "solve (Homotopy)" begin
        @var x a y b
        F = System([x^2 - a, x * y - a + b], [x, y], [a, b])
        s = [1, 1]
        H = ParameterHomotopy(F, [1, 0], [2, 4])
        res = solve(H, [s])
        @test nsolutions(res) == 1
    end

    @testset "solve (start target)" begin
        @var x a y b
        f = System([x^2 - a, x * y - a + b], parameters = [a, b])
        s = [1, 1]
        res = solve(f, f, [s]; start_parameters = [1, 0], target_parameters = [2, 4])
        @test nsolutions(res) == 1

        G = FixedParameterSystem(f, [1, 0])
        F = FixedParameterSystem(f, [2, 4])
        res = solve(G, F, [s];)
        @test nsolutions(res) == 1
    end

    @testset "solve (Vector{Expression})" begin
        @var x a y b
        F = [x^2 - a, x * y - a + b]
        s = [1, 1]
        res = solve(
            F,
            [s];
            parameters = [a, b],
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        @test nsolutions(res) == 1
    end

    @testset "solve (DynamicPolynomials)" begin
        @polyvar x y
        # define the polynomials
        f₁ = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        f₂ = x^2 + 2x * y^2 - 2 * y^2 - 1 / 2
        result = solve([f₁, f₂])
        @test nsolutions(result) == 18

        @polyvar x a y b
        F = [x^2 - a, x * y - a + b]
        s = [1, 1]
        res = solve(
            F,
            [s];
            variables = [x, y],
            parameters = [a, b],
            start_parameters = [1, 0],
            target_parameters = [2, 4],
            compile = false,
        )
        @test nsolutions(res) == 1
        res2 = solve(
            F,
            [s];
            variable_ordering = [y, x],
            parameters = [a, b],
            start_parameters = [1, 0],
            target_parameters = [2, 4],
            compile = false,
        )
        s = solutions(res)[1]
        s2 = solutions(res2)[1]
        @test s ≈ [s2[2], s2[1]]
    end

    @testset "change parameters" begin
        @var x a y b
        F = System([x^2 - a, x * y - a + b]; parameters = [a, b])
        s = [1.0, 1.0 + 0im]
        S = solver(F, generic_parameters = [2.2, 3.2])
        start_parameters!(S, [1, 0])
        target_parameters!(S, [2, 4])
        @test is_success(track(S, s))
    end

    @testset "solve (threading)" begin
        res = solve(cyclic(7), threading = true, show_progress = false)
        @test nsolutions(res) == 924
    end

    @testset "stop early callback" begin
        @var x
        first_result = nothing
        results = solve(
            [(x - 3) * (x + 6) * (x + 2)];
            stop_early_cb = r -> begin
                first_result = r
                true
            end,
            threading = false,
            show_progress = false,
            start_system = :total_degree,
        )
        @test length(results) == 1
        @test first(results) === first_result

        nresults = 0
        result = solve(
            [(x - 3) * (x + 6) * (x + 2)],
            stop_early_cb = r -> (nresults += 1) == 2,
            start_system = :total_degree,
            show_progress = false,
            threading = false,
        )
        @test length(result) == 2

        # threading
        @var x y z
        first_result = nothing
        # this has 5^3 = 125 solutions, so we should definitely stop early if we
        # have less than 64 threads
        results = solve(
            [
                (x - 3) * (x + 6) * (x + 2) * (x - 2) * (x + 2.5),
                (y + 2) * (y - 2) * (y + 3) * (y + 5) * (y - 1),
                (z + 2) * (z - 2) * (z + 3) * (z + 5) * (z - 2.1),
            ];
            stop_early_cb = r -> begin
                first_result = r
                true
            end,
            threading = true,
            show_progress = false,
            start_system = :total_degree,
        )
        @test length(results) < 125
    end

    @testset "Many parameters solver" begin
        # Setup
        @var x y
        f = x^2 + y^2 - 1

        @var a b c
        l = a * x + b * y + c
        F = [f, l]

        # Compute start solutions S₀ for given start parameters p₀
        p₀ = randn(ComplexF64, 3)
        S₀ = solutions(solve(subs(F, [a, b, c] => p₀)))
        # The parameters we are intersted in
        params = [rand(3) for i = 1:100]

        result1 = solve(
            F,
            S₀,
            ;
            start_parameters = p₀,
            target_parameters = params,
            parameters = [a, b, c],
            threading = true,
        )
        @test typeof(result1) == Vector{Tuple{Result,Vector{Float64}}}
        result1 = solve(
            F,
            S₀,
            ;
            start_parameters = p₀,
            target_parameters = params,
            parameters = [a, b, c],
            show_progress = false,
            threading = false,
        )
        @test typeof(result1) == Vector{Tuple{Result,Vector{Float64}}}

        # Only keep real solutions
        result2 = solve(
            F,
            S₀,
            ;
            start_parameters = p₀,
            target_parameters = params,
            parameters = [a, b, c],
            transform_result = (r, p) -> real_solutions(r),
            threading = true,
        )
        @test typeof(result2) == Vector{Vector{Vector{Float64}}}
        @test !isempty(result2)

        # Now instead of an Array{Array{Array{Float64,1},1},1} we want to have an
        # Array{Array{Float64,1},1}
        result3 = solve(
            F,
            S₀,
            ;
            start_parameters = p₀,
            target_parameters = params,
            parameters = [a, b, c],
            transform_result = (r, p) -> real_solutions(r),
            flatten = true,
            threading = false,
        )
        @test typeof(result3) == Vector{Vector{Float64}}
        @test !isempty(result3)

        # The passed `params` do not directly need to be the target parameters.
        # Instead they can be some more concrete informations (e.g. an index)
        # and we can them by using the `transform_parameters` method
        result4 = solve(
            F,
            S₀,
            ;
            start_parameters = p₀,
            target_parameters = 1:100,
            parameters = [a, b, c],
            transform_result = (r, p) -> (real_solutions(r), p),
            transform_parameters = _ -> rand(3),
        )
        @test typeof(result4) == Vector{Tuple{Vector{Vector{Float64}},Int64}}
    end
end
