@testset "solve" begin

    @testset "total degree (simple)" begin
        @var x y
        affine_square = System([
            2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3,
            2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5,
        ])
        @test count(is_success, track.(total_degree(affine_square)...)) == 2

        @var x y z
        proj_square = System([
            2.3 * x^2 + 1.2 * y^2 + 3x * z - 2y * z + 3 * z^2,
            2.3 * x^2 + 1.2 * y^2 + 5x * z + 2y * z - 5 * z^2,
        ])
        @test count(is_success, track.(total_degree(proj_square)...)) == 4

        @var x y
        affine_overdetermined = System([
            (x^2 + y^2 + x * y - 3) * (x + 3),
            (x^2 + y^2 + x * y - 3) * (y - x + 2),
            2x + 5y - 3,
        ])
        @test count(is_success, track.(total_degree(affine_overdetermined)...)) == 2

        @var x y
        affine_overdetermined_reordering = System([
            (x^2 + y^2 + x * y - 3) * (x + 3),
            2x + 5y - 3,
            (x^2 + y^2 + x * y - 3) * (y^2 - x + 2),
        ])
        tracker, starts = total_degree(affine_overdetermined_reordering)
        @test length(starts) == 4 * 3
        @test count(is_success, track.(tracker, starts)) == 2

        @var x y z
        proj_overdetermined = System([
            (x^2 + y^2 + x * y - 3 * z^2) * (x + 3z),
            (x^2 + y^2 + x * y - 3 * z^2) * (y - x + 2z),
            2x + 5y - 3z,
        ])
        @test count(is_success, track.(total_degree(proj_overdetermined)...)) == 2

        @var x y
        proj_overdetermined_reordering = System([
            (x^2 + y^2 + x * y - 3 * z^2) * (x + 3z),
            2x + 5y - 3z,
            (x^2 + y^2 + x * y - 3 * z^2) * (y^2 - x * z + 2 * z^2),
        ])
        tracker, starts = total_degree(proj_overdetermined_reordering)
        @test length(starts) == 4 * 3
        @test count(is_success, track.(tracker, starts)) == 2

        @var x y
        affine_underdetermined = System([2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3])
        @test_throws HC2.FiniteException total_degree(affine_underdetermined)

        @var x y z
        proj_underdetermined = System([2.3 * x^2 + 1.2 * y^2 + 3x * z])
        @test_throws HC2.FiniteException total_degree(proj_underdetermined)
    end

    @testset "total degree (variable groups)" begin
        @var x y v w
        affine_sqr = System([x * y - 2, x^2 - 4], variable_groups = [[x], [y]])
        tracker, starts = total_degree(affine_sqr)
        @test length(collect(starts)) == 2
        @test count(is_success, track.(tracker, starts)) == 2

        @var x y v w
        proj_sqr =
            System([x * y - 2v * w, x^2 - 4 * v^2], variable_groups = [[x, v], [y, w]])
        tracker, starts = total_degree(proj_sqr)
        @test length(collect(starts)) == 2
        @test count(is_success, track.(tracker, starts)) == 2

        @var x y v w
        affine_overdetermined = System(
            [(x^2 - 4) * (x * y - 2), x * y - 2, x^2 - 4],
            variable_groups = [[x], [y]],
        )
        tracker, starts = total_degree(affine_overdetermined)
        @test count(is_success, track.(tracker, starts)) == 2
        @var x y v w
        proj_overdetermined = System(
            [(x^2 - 4 * v^2) * (x * y - v * w), x * y - v * w, x^2 - v^2],
            variable_groups = [[x, v], [y, w]],
        )
        tracker, starts = total_degree(proj_overdetermined)
        @test count(is_success, track.(tracker, starts)) == 2
    end

    @testset "polyhedral" begin
        @var x y
        affine_square = System([
            2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3,
            2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5,
        ])
        @test count(is_success, track.(polyhedral(affine_square)...)) == 2

        @var x y z
        proj_square = System([
            2.3 * x^2 + 1.2 * y^2 + 3x * z - 2y * z + 3 * z^2,
            2.3 * x^2 + 1.2 * y^2 + 5x * z + 2y * z - 5 * z^2,
        ])
        @test count(is_success, track.(polyhedral(proj_square)...)) == 4

        @var x y
        affine_overdetermined = System([
            (x^2 + y^2 + x * y - 3) * (x + 3),
            (x^2 + y^2 + x * y - 3) * (y - x + 2),
            2x + 5y - 3,
        ])
        @test count(is_success, track.(polyhedral(affine_overdetermined)...)) == 2

        @var x y z
        proj_overdetermined = System([
            (x^2 + y^2 + x * y - 3 * z^2) * (x + 3z),
            (x^2 + y^2 + x * y - 3 * z^2) * (y - x + 2z),
            2x + 5y - 3z,
        ])
        @test count(is_success, track.(polyhedral(proj_overdetermined)...)) == 2

        @var x y
        affine_underdetermined = System([2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3])
        @test_throws HC2.FiniteException polyhedral(affine_underdetermined)

        @var x y z
        proj_underdetermined = System([2.3 * x^2 + 1.2 * y^2 + 3x * z])
        @test_throws HC2.FiniteException polyhedral(proj_underdetermined)
    end

    @testset "overdetermined" begin
        @testset "3 by 5 minors" begin
            res = track.(total_degree(minors())...)
            @test count(is_success, res) == 80
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

        seeded_res = solve(
            F;
            target_parameters = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32],
            seed = seed(res),
        )
        @test seed(seeded_res) == seed(res)

        @test length(path_results(res)) == ntracked(res) == 7
        @test length(results(res)) == nresults(res) == 3
        @test length(solutions(res)) == 3
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
    end

    @testset "composition" begin
        @var a b c x y z u v
        e = System([u + 1, v - 2])
        f = System([a * b - 2, a * c - 1])
        g = System([x + y, y + 3, x + 2])

        res = solve(e ∘ f ∘ g; start_system = :total_degree)
        @test nsolutions(res) == 2

        res = solve(e ∘ f ∘ g; start_system = :polyhedral)
        @test nsolutions(res) == 2
    end

    @testset "paths to track" begin
        @var x y
        f = System([2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y])

        @test paths_to_track(f; start_system = :total_degree) == 16
        @test paths_to_track(f; start_system = :polyhedral) == 8
        @test paths_to_track(f; start_system = :polyhedral, only_non_zero = true) == 3
        @test paths_to_track(f) == 8
        @test_deprecated bezout_number(f) == 16
        @test mixed_volume(f) == 3
    end

    @testset "solve (parameter homotopy)" begin
        # affine
        @var x a y b
        F = System([x^2 - a, x * y - a + b], [x, y], [a, b])
        s = [1, 1]
        res = solve(F, [s]; start_parameters = [1, 0], target_parameters = [2, 4])
        @test nsolutions(res) == 1
        res = solve(
            ModelKitSystem(F),
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
        res = solve(ModelKitSystem(F_proj), [s]; p₁ = [1, 0], p₀ = [2, 4])
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
        res = solve(ModelKitSystem(F_multi_proj), S; p₁ = [2, 4], p₀ = [3, 5])
        @test nsolutions(res) == 2

        F_multi_proj_err = System(
            [x * y - a * v * w],
            parameters = [a, b],
            variable_groups = [[x, v], [y, w]],
        )
        @test_throws FiniteException(1) solve(F_multi_proj_err, S; p₁ = [2, 4], p₀ = [3, 5])
    end

    @testset "Solve (Homotopy)" begin
        @var x a y b
        F = System([x^2 - a, x * y - a + b], [x, y], [a, b])
        s = [1, 1]
        H = ParameterHomotopy(F, [1, 0], [2, 4])
        res = solve(H, [s])
        @test nsolutions(res) == 1
    end

    @testset "Solve (threading)" begin
        res = solve(cyclic(7), threading = true, show_progress = false)
        @test nsolutions(res) == 924
    end
end
