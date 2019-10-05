@testset "Polyhedral" begin
    @testset "Start solutions iterator" begin
        f = equations(cyclic(5))
        iter = HC.PolyhedralStartSolutionsIterator(f)
        @test isnothing(iter.mixed_cells)
        @test isnothing(iter.lifting)
        @test length(iter) == 70
        @test !isnothing(iter.mixed_cells)
        @test !isnothing(iter.lifting)
        starts = collect(iter)
        @test length(starts) == 70
        @test all(starts) do (cell, x)
            all(1:length(cell.indices)) do i
                a, b = cell.indices[i]
                v_i = iter.start_coefficients[i][a] * prod(x .^ iter.support[i][:, a]) +
                      iter.start_coefficients[i][b] * prod(x .^ iter.support[i][:, b])
                abs(v_i) < 1e-12
            end
        end
        @test length(collect(HC.PolyhedralStartSolutionsIterator(f))) == 70

        Random.seed!(130793) # this seed should result in the use of 128 bit hnf
        f = equations(PolynomialTestSystems.tritangents())
        iter = HC.PolyhedralStartSolutionsIterator(f)
        @test length(iter) == 12636
        @test length(collect(iter)) == 12636

        # Catch overflow in mixed_volume
        f = equations(PolynomialTestSystems.cyclooctane())
        F = [f; randn(2, 18) * [HC.variables(f); 1]]
        iter = HC.PolyhedralStartSolutionsIterator(F)
        @test_throws OverflowError length(iter)
    end

    @testset "Tracking" begin
        tracker, starts = pathtracker_startsolutions(
            equations(cyclic(5));
            start_system = :polyhedral,
            seed = 123122,
        )
        S = collect(starts)
        HC.prepare!(tracker, starts)
        @test count(s -> is_success(track!(tracker, s)), S) == 70

        @polyvar x y
        tracker, starts = pathtracker_startsolutions(
            [(x - 3), (y - 2)];
            start_system = :polyhedral,
            system_scaling = nothing,
            seed = 23121,
        )
        S = collect(starts)
        HC.prepare!(tracker, starts)
        @test solution(track(tracker, S[1])) ≈ [3, 2] atol = 1e-8

        tracker, starts = pathtracker_startsolutions(
            [(x - 3), (y - 2)];
            start_system = :polyhedral,
            system_scaling = :equations_and_variables,
            seed = 23121,
        )
        S = collect(starts)
        HC.prepare!(tracker, starts)
        @test solution(track(tracker, S[1])) ≈ [3, 2] atol = 1e-8
    end

    @testset "All affine solutions" begin
        @polyvar x y
        f = [2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y]
        tracker, starts = pathtracker_startsolutions(
            f,
            start_system = :polyhedral,
            only_torus = true,
        )
        S = collect(starts)
        HC.prepare!(tracker, starts)
        results = map(s -> track(tracker, s), S)
        @test count(is_success, results) == 3
        @test length(results) == 3

        tracker, starts = pathtracker_startsolutions(
            f,
            start_system = :polyhedral,
            only_torus = false,
        )
        S = collect(starts)
        HC.prepare!(tracker, starts)
        results = map(s -> track(tracker, s), S)
        @test count(is_success, results) == 6
        @test length(results) == 8
    end

end
