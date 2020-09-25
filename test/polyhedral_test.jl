@testset "Polyhedral" begin
    @testset "affine + torus solutions" begin
        @var x y
        f = System([2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y])

        tracker, starts = polyhedral(f; only_torus = false)
        @test length(collect(starts)) == 8
        @test count(is_success, track.(tracker, starts)) == 6

        tracker, starts = polyhedral(f; only_torus = true)
        @test length(collect(starts)) == 3
        @test count(is_success, track.(tracker, starts)) == 3

        tracker, starts = polyhedral(
            f;
            only_torus = true,
            tracker_options = (automatic_differentiation = 3,),
        )
        @test length(collect(starts)) == 3
        @test count(is_success, track.(tracker, starts)) == 3
    end

    @testset "only torus" begin
        @var x₁ x₂ s

        HC_I = [
            x₁^3 * s^15 + x₁ * x₂ * s + x₂^3 + s^12,
            x₁^2 * s^9 + x₁ * x₂^2 + x₂ * s^3,
            x₁^2 * x₂ * s^5 + x₁ * s^8 + x₂^2,
        ]
        F = System(HC_I)
        _, starts = solver_startsolutions(F; start_system = :polyhedral, only_torus = false)
        @test length(collect(starts)) == paths_to_track(F; only_torus = false) == 92
        _, starts = solver_startsolutions(F; start_system = :polyhedral, only_torus = true)
        @test length(collect(starts)) == paths_to_track(F; only_torus = true) == 54
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
