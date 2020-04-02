@testset "PathTracker" begin
    @testset "Tracking and PathResult" begin
        @var x y
        f = [2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3, 2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5]
        H, starts = total_degree_homotopy(f, [x, y]; gamma = 1.3im + 0.4)
        S = collect(starts)
        tracker = PathTracker(Tracker(H))
        @test !isempty(sprint(show, tracker.options))

        res = track.(tracker, S)

        @test is_success(res[1])
        @test all(isfinite, solution(res[1]))
        @test is_finite(res[1])
        @test winding_number(res[1]) == nothing
        @test !is_failed(res[1])
        @test accuracy(res[1]) < 1e-12
        @test steps(res[1]) < 10
        @test accepted_steps(res[1]) < 10
        @test rejected_steps(res[1]) == 0
        @test is_success(res[2])
        @test is_at_infinity(res[3])
        @test is_at_infinity(res[4])
        @test !isempty(sprint(show, res[1]))
        @test valuation(res[3]) ≈ [-1, -1] rtol = 1e-3

        # singular stuff
        @var x
        f = [(x - 10)^2]
        H, starts = total_degree_homotopy(f, [x])
        tracker = PathTracker(Tracker(H))
        res = track.(tracker, starts)

        @test winding_number(res[1]) == 2
        @test last_path_point(res[1]) isa Tuple{Vector{ComplexF64},Float64}
        @test 0 < last(last_path_point(res[1])) < 0.1
    end

    @testset "Overdetermined tracking" begin
        @var x y z

        p₁ = (x^2 + y^2 + z^2 - 1) * (x - 0.5)
        p₂ = (x^2 + y^2 + z^2 - 1) * (y - 0.5)
        p₃ = (z - x^2 - 2) * (x^2 + y^2 + z^2 - 1) * (z - 0.5)
        F = System([p₁, p₂, p₃])

        L₁ = rand_affine_subspace(3; codim = 2)
        L₂ = rand_affine_subspace(3; codim = 2)

        H, S =
            total_degree_homotopy(System([x^2 + y^2 + z^2 - 1]) ∩ extrinsic(L₁)([x, y, z]))
        res = track.(Tracker(H), S)
        @test all(is_success, res)

        H = AffineSubspaceHomotopy(F, L₁, L₂)
        res2 = track.(PathTracker(H), solution.(res))
        @test all(is_success, res2)
    end
end
