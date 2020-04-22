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
            [(x^2 - 4v^2) * (x * y - v * w), x * y - v*w, x^2 - v^2],
            variable_groups = [[x,v], [y,w]],
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

    # @testset "Automatic start systems (solve)" begin
    #     @var x y
    #     affine_square = System([
    #         2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3,
    #         2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5,
    #     ])
    #
    #     solve(affine_square; start_system = :polyhedral)
    #     solve(affine_square; start_system = :total_degree)
    #
    #     @var x y z
    #     proj_square = System([
    #         2.3 * x^2 + 1.2 * y^2 + 3x * z - 2y * z + 3 * z^2,
    #         2.3 * x^2 + 1.2 * y^2 + 5x * z + 2y * z - 5 * z^2,
    #     ])
    #     solve(proj_square, start_system = :total_degree)
    #     solve(proj_square, start_system = :polyhedral)
    #
    #     #TODO: Multi-proj square
    #
    #     @var x y
    #     affine_overdetermined = System([
    #         (x^2 + y^2 + x * y - 3) * (x + 3),
    #         (x^2 + y^2 + x * y - 3) * (y - x + 2),
    #         2x + 5y - 3,
    #     ])
    #
    #     solve(affine_overdetermined, start_system = :total_degree)
    #     solve(affine_overdetermined, start_system = :polyhedral)
    #
    #     # TODO: Test reordering of degrees
    #     @var x y z
    #     proj_overdetermined = System([
    #         (x^2 + y^2 + x * y - 3 * z^2) * (x + 3z),
    #         (x^2 + y^2 + x * y - 3 * z^2) * (y - x + 2z),
    #         2x + 5y - 3z,
    #     ])
    #
    #     solve(proj_overdetermined, start_system = :total_degree)
    #     solve(proj_overdetermined, start_system = :polyhedral)
    #
    #     @var x y
    #     affine_underdetermined = System([2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3])
    #     @test_throws HC2.FiniteException solve(affine_underdetermined)
    # end
end
