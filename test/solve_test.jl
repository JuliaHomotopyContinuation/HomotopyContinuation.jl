@testset "solve - invalid input" begin
    @polyvar x y z
    @test_throws AssertionError solve([x-2y+2, 0])
    @test_throws AssertionError solve([x-2, y-2], [x-2, y-2,y+2], [[2, -3]])

    # non homogenous overdetermiend
    @test_throws AssertionError solve([x-2z, y^2+3z, z^3+x^3], homvar=z)
    @test_throws AssertionError solve([x-2z, y^2+3z, z^3+x, z+x])
    # homogenous overdetermiend
    @test_throws AssertionError solve([x-2z, y^2+3z, z^3+x, z+x])
    @test_throws AssertionError solve(Systems.FPSystem([x-2z, y^2+3z^2, z^3+x^3, z+x]))

    @test_throws AssertionError solve([x-2z, y^2+3z^2, z^3+x^3])
end

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

    @polyvar w
    F = equations(cyclic(5))
    result = solve(Utilities.homogenize(F, w), homvar=w)
    @test result isa AffineResult

    @polyvar x y z

    G = [x-2, y+3]
    F = [x+2, y-3]
    @test nfinite(solve(G, F, [[2, -3]])) == 1
    @test nfinite(solve(G, F, [[2+0.0im, -3.0+0im]])) == 1

    F = Systems.FPSystem(Utilities.homogenize(equations(cyclic(5))))
    result = solve(F, homvar=6)
    @test nfinite(result) == 70
    @test natinfinity(result) == 50


    F = equations(katsura(5))
    prob, startsolutions = Problems.problem_startsolutions(Input.TotalDegree(F))

    result = solve(prob.homotopy, map(startsolutions) do s
        ProjectiveVectors.raw(Problems.embed(prob, s))
    end)
    @test result isa ProjectiveResult
    @test nnonsingular(result) == 32

    result = solve(prob.homotopy, map(startsolutions) do s
        ProjectiveVectors.raw(Problems.embed(prob, s))
    end, homvar=prob.homogenization.homvaridx)
    @test result isa AffineResult
    @test nnonsingular(result) == 32
    @test nfinite(result) == 32


    # test numerical homogenous check fails
    @polyvar x y z
    G = Systems.FPSystem([x-2z, y^2+3z])
    @test_throws ErrorException solve(G, homvar=3)

    # invalid homvar
    @test_throws AssertionError solve([x-2z, y^2+3z], homvar=z)
    # cyclic_hom = Utilities.homogenize(equations(cyclic(5)))
    # F = Systems.FPSystem(cyclic_hom)
    # G = Systems.FPSystem(Problems.totaldegree(cyclic_hom, Utilities.allvariables(cyclic_hom), Utilities.allvariables(cyclic_hom)[end]))
    # result = solve(G, F)
end


@testset "solve - random seed" begin
    R = solve(equations(katsura(5)), seed=1234)
    @test seed(R) == 1234

    @polyvar x y

    G = [x-2, y+3]
    F = [x+2, y-3]
    @test seed(solve(G, F, [[2, -3]], seed=222)) == 222
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
    P, start_sols = Problems.problem_startsolutions(Input.TotalDegree(F))
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
    @test nfinite(solve(F, tol=1e-3, seed=2337, threading=false)) < 158
    @test nfinite(solve(F, tol=1e-3, seed=2337)) < 158
end

@testset "Affine vs projective" begin
    @polyvar x y z
    f = [x-2y, y-2z]
    g = [x-2, y-2]

    F = solve(f)
    G = solve(g)

    @test length(F[1].solution) == 3
    @test length(G[1].solution) == 2
    @test F[1].solution_type == :projective
    @test G[1].solution_type == :affine
end

@testset "Parameter Homotopies" begin
    @polyvar x a y b
    F = [x^2-a, x*y-a+b]
    p = [a, b]
    S = solve(F, p, [1, 0], [2, 4], [[1.0, 1.0 + 0.0*im]])

    @test S[1].solution ≈ [complex(√2), -complex(√2)]
    @test nfinite(S) == 1

    @polyvar x a y b z
    F = [x^2-a*z^2, x*y-(a-b)*z^2]
    p = [a, b]
    S = solve(F, p, [1, 0], [2, 4], [[1.0, 1.0 + 0.0*im, 1.0]])
    @test S isa Solving.ProjectiveResult
    @test solution(S[1])[1:2] / solution(S[1])[3] ≈ [complex(√2), -complex(√2)]
    @test nnonsingular(S) == 1

    S2 = solve(F, p, [1, 0], [2, 4], [[1.0, 1.0 + 0.0*im, 1.0]], homvar=z)
    @test solution(S2[1]) ≈ [complex(√2), -complex(√2)]
    @test nfinite(S2) == 1
end

@testset "Overdetermined" begin
    @polyvar x y z w
    a = [0.713865+0.345317im, 0.705182+0.455495im, 0.9815+0.922608im, 0.337617+0.508932im]

    f = [x*z-y^2, y*w-z^2, x*w-y*z]
    L₁ = [1, -1, 1, -1] ⋅ [x, y, z, w]
    L₂ = rand(Complex128, 4) ⋅ [x, y, z, w]
    S = solve([f; L₁], [f; L₂], [[1, 1, 1, 1]])

    @test nnonsingular(S) == 1


    f = [x*z-y^2, y-z^2, x-y*z]
    L₁ = [1, -1, 1, -1] ⋅ [x, y, z, 1]
    L₂ = rand(Complex128, 4) ⋅ [x, y, z, 1]
    S = solve([f; L₁], [f; L₂], [[1, 1, 1]])

    @test nfinite(S) == 1
end
