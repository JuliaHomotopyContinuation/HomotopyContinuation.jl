@testset "solve" begin
    @polyvar x
    @test nfinite(solve([x - 1])) == 1
    F = equations(katsura(5))
    @test nfinite(solve(F, threading=false)) == 32
    @test nfinite(solve(F, system=Systems.SPSystem, threading=false)) == 32
    @test nfinite(solve(F, system=Systems.FPSystem, threading=false)) == 32

    @test nfinite(solve(F, homotopy=Homotopies.StraightLineHomotopy)) == 32
    result = solve(F, predictor=Predictors.Euler(), homotopy=Homotopies.StraightLineHomotopy)
    @test nresults(result) == 32
    @test nfinite(solve(F, tol=1e-5)) == 32

    result = solve(F)
    @test nfinite(result) == 32
    @test string(result) isa String

    @test nfinite(solve(F, patch=AffinePatches.RandomPatch())) == 32
    @test nfinite(solve(F, patch=AffinePatches.EmbeddingPatch())) ≤ 32
    @test nfinite(solve(F, patch=AffinePatches.OrthogonalPatch())) == 32
end

@testset "solve - no endgame" begin
    F = equations(katsura(5))
    # no endgame
    @test nfinite(solve(F, endgame_start=0.0)) == 32

    @test nfinite(solve(F, endgame_start=0.0, threading=false)) == 32
end

@testset "solve - kwargs" begin
    @test_throws ErrorException solve(equations(cyclic(5)), def=0.4, abc=23)
end

@testset "Singular solutions" begin
    @polyvar x y z
    z = 1
    F = [x^2 + 2*y^2 + 2*im*y*z, (18 + 3*im)*x*y + 7*im*y^2 - (3 - 18*im)*x*z - 14*y*z - 7*im*z^2]
    result = solve(F)
    @test nsingular(result) == 3
    @test all(r -> r.windingnumber == 3, singular(result))
end

@testset "Path Crossing" begin
    # This tests that we indeed detect path crossings
    srand(2337)
    F = equations(cyclic(6))
    P = Problems.ProjectiveStartTargetProblem(Problems.TotalDegreeProblem(F))
    start_sols = Utilities.totaldegree_solutions(F) |> collect
    x₁ = Problems.embed(P, start_sols[1])
    H = Homotopies.PatchedHomotopy(P.homotopy, AffinePatches.OrthogonalPatch(), x₁)
    tracker = PathTracking.PathTracker(H, x₁, 1.0, 0.1, tol=1e-3)
    tracked_paths = map(start_sols) do x
        PathTracking.track(tracker, Problems.embed(P, x), 1.0, 0.1)
    end

    crossed_path_indices = Solving.check_crossed_paths(tracked_paths, 1e-2)
    @test length(crossed_path_indices) > 0

    tracker = PathTracking.PathTracker(H, x₁, 1.0, 0.1, tol=1e-8)
    tracked_paths = map(start_sols) do x
        PathTracking.track(tracker, Problems.embed(P, x), 1.0, 0.1)
    end
    crossed_path_indices = Solving.check_crossed_paths(tracked_paths, 1e-7)
    @test isempty(crossed_path_indices)

    # Test that we resolve path crossings
    F = equations(cyclic(6))
    # this will have three crossed paths
    @test nfinite(solve(F, tol=1e-3, seed=2337, threading=false)) ≤ 156
    @test nfinite(solve(F, tol=1e-3, seed=2337)) ≤ 156
end
