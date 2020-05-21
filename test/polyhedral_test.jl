@testset "Polyhedral" begin
    @testset "affine + torus solutions" begin
        @var x y
        f = System([2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y])

        tracker, starts = polyhedral(f; only_torus = false)
        @test length(starts) == 8
        @test count(is_success, track.(tracker, starts)) == 6

        tracker, starts = polyhedral(f; only_torus = true)
        @test length(starts) == 3
        @test count(is_success, track.(tracker, starts)) == 3

        tracker, starts = polyhedral(
            f;
            only_torus = true,
            tracker_options = (automatic_differentiation = 3,),
        )
        @test length(starts) == 3
        @test count(is_success, track.(tracker, starts)) == 3
    end

    @testset "cyclic" begin
        tracker, starts = polyhedral(cyclic(5))
        res = track.(tracker, starts)
        @test count(is_success, res) == 70

        tracker, starts = polyhedral(cyclic(7))
        res = track.(tracker, starts)
        @test count(is_success, res) == 924
    end
end
