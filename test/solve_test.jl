@testset "solve" begin
    @testset "Invalid input" begin
        @polyvar x y z
        @test_throws ErrorException solve([x-2y+2, 0])
        @test_throws ErrorException solve([x-2, y-2], [x-2, y-2,y+2], [[2, -3]])

        # non homogenous overdetermiend
        @test_throws ErrorException solve([x-2z, y^2+3z, z^3+x^3], homvar=z)
        @test_throws ErrorException solve([x-2z, y^2+3z, z^3+x, z+x])
        # homogenous overdetermiend
        @test_throws ErrorException solve([x-2z, y^2+3z^2, z^3+x^3, z+x])
        @test_throws ErrorException solve([x-2z, y^2+3z^2, z^3+x^3, z+x], homvar=z)
        @test_throws ErrorException solve(Systems.FPSystem([x-2z, y^2+3z^2, z^3+x^3, z+x]))
        @test_throws ErrorException solve(Systems.FPSystem([x-2z, y^2+3z^2, z^3+x^3, z+x]), homvar=4)

        @test_throws ErrorException solve([x-2z, y^2+3z^2, z^3+x^3])

        @test_throws ErrorException solve([x-2z, y^2+3z], homvar=z)

        # invalid kwargs
        @test_throws ErrorException solve(equations(cyclic(5)), def=0.4, abc=23)

        # test numerical homogenous check fails
        @polyvar x y z
        G = Systems.FPSystem([x-2z, y^2+3z])
        @test_throws ErrorException solve(G, homvar=3)
    end

    @testset "solve" begin
        @polyvar x
        @test nfinite(solve([x - 1])) == 1
        F = equations(katsura(5))
        @test nfinite(solve(F, threading=false)) == 32
        @test nfinite(solve(F, system=Systems.SPSystem, threading=false)) == 32
        @test nfinite(solve(F, system=Systems.FPSystem, threading=false)) == 32

        # Simple step size
        F = equations(katsura(5))
        @test nfinite(solve(F, simple_step_size=true, threading=false)) == 32

        # scaling
        F = equations(katsura(5))
        @test nfinite(solve(F, scale_systems=false, threading=false)) == 32

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
        result = solve(Utilities.homogenize(F, w), threading=false, homvar=w)
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
            Problems.embed(prob, s).data
        end)
        @test result isa ProjectiveResult
        @test nnonsingular(result) == 32

        result = solve(prob.homotopy, map(startsolutions) do s
            Problems.embed(prob, s).data
        end, homvar=Problems.homvars(prob)[1])
        @test result isa AffineResult
        @test nnonsingular(result) == 32
        @test nfinite(result) == 32

        # Composition
        @polyvar a b c x y z u v
        e = [u + 1, v - 2]
        f = [a * b - 2, a*c- 1]
        g = [x+y, y + 3, x + 2]
        res = solve(e ∘ f ∘ g)
        @test nnonsingular(res) == 2
        @test nnonsingular(solve(e ∘ f ∘ g, scale_systems=false)) == 2

        res = solve(e ∘ f ∘ g, system=SPSystem)
        @test nnonsingular(res) == 2
    end


    @testset "Random seed" begin
        R = solve(equations(katsura(5)), seed=1234)
        @test seed(R) == 1234

        @polyvar x y

        G = [x-2, y+3]
        F = [x+2, y-3]
        @test seed(solve(G, F, [[2, -3]], seed=222)) == 222
    end

    @testset "No Endgame" begin
        F = equations(katsura(5))
        # no endgame
        @test nfinite(solve(F, endgame_start=0.0)) == 32

        @test nfinite(solve(F, endgame_start=0.0, threading=false)) == 32
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
        F = equations(cyclic(6))
        tracker, start_sols = pathtracker_startsolutions(F, tol=1e-3, corrector_maxiters=5, seed=123512)
        tracked_paths = map(start_sols) do x
            PathTracking.track(tracker, x, 1.0, 0.1)
        end

        crossed_path_indices = Solving.check_crossed_paths(tracked_paths, 1e-2)
        @test length(crossed_path_indices) > 0

        tracker, start_sols = pathtracker_startsolutions(F, seed=123512)
        tracked_paths = map(start_sols) do x
            PathTracking.track(tracker, x, 1.0, 0.1)
        end
        crossed_path_indices = Solving.check_crossed_paths(tracked_paths, 1e-5)
        @test isempty(crossed_path_indices)
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
        S = solve(F, [[1.0, 1.0 + 0.0*im]], parameters=p, p₁=[1, 0], p₀=[2, 4])
        # S = solve(F, p, [1, 0], [2, 4], [[1.0, 1.0 + 0.0*im]])

        @test S[1].solution ≈ [complex(√2), -complex(√2)]
        @test nfinite(S) == 1

        @polyvar x a y b z
        F = [x^2-a*z^2, x*y-(a-b)*z^2]
        S = solve(F, [[1.0, 1.0 + 0.0*im, 1.0]], parameters=[a, b], startparameters=[1, 0], targetparameters=[2, 4])
        @test S isa Solving.ProjectiveResult
        @test solution(S[1])[1:2] / solution(S[1])[3] ≈ [complex(√2), -complex(√2)]
        @test nnonsingular(S) == 1

        S2 = solve(F, [[1.0, 1.0 + 0.0*im, 1.0]], parameters=[a, b], p₁=[1, 0], p₀=[2, 4], homvar=z)
        @test solution(S2[1]) ≈ [complex(√2), -complex(√2)]
        @test nfinite(S2) == 1

        γ₁, γ₀ =randn(ComplexF64, 2)
        S2 = solve(F, [[1.0, 1.0 + 0.0*im, 1.0]], parameters=[a, b], p₁=[1, 0], p₀=[2, 4], γ₁=γ₁, γ₀=γ₀, homvar=z)
        @test solution(S2[1]) ≈ [complex(√2), -complex(√2)]
        @test nfinite(S2) == 1
    end

    @testset "ParameterHomotopy with Composition" begin
        @polyvar p q a b c x y z u v
        f = [a * b - 2, a*c- 1]
        g = [x + y, y + 3, x + 2]
        res = solve(f ∘ g, system=SPSystem)

        # parameters at the end
        f2 = [a * b - q, a*c- p]
        g = [x + y, y + 3, x + 2]
        r = solve(f2 ∘ g, solutions(res), parameters=[p, q], p₁=[1, 2], p₀=[2, 3])
        @test HomotopyContinuation.nnonsingular(r) == 2

        # parameters at the beginning
        f = [a * b - 2, a*c- 1]
        g2 = [x + y, y + u, x + v]
        r = solve(f ∘ g2, solutions(res), parameters=[u, v], p₁=[3, 2], p₀=[-2, 3])
        @test HomotopyContinuation.nnonsingular(r) == 2

        # parameter in the middle
        e = [u + 1, v - 2]
        res2 = solve(e ∘ f ∘ g, system=SPSystem)
        f2 = [a * b - q, a * c- p]
        r = solve(e ∘ f2 ∘ g, solutions(res2), parameters=[p, q], p₁=[1, 2], p₀=[2, 3])
        @test HomotopyContinuation.nnonsingular(r) == 2
    end

    @testset "Overdetermined" begin
        @polyvar x y z w
        a = [0.713865+0.345317im, 0.705182+0.455495im, 0.9815+0.922608im, 0.337617+0.508932im]

        f = [x*z-y^2, y*w-z^2, x*w-y*z]
        L₁ = [1, -1, 1, -1] ⋅ [x, y, z, w]
        L₂ = rand(ComplexF64, 4) ⋅ [x, y, z, w]
        S = solve([f; L₁], [f; L₂], [[1, 1, 1, 1]])
        @test nnonsingular(S) == 1

        f = [x*z-y^2, y-z^2, x-y*z]
        L₁ = [1, -1, 1, -1] ⋅ [x, y, z, 1]
        L₂ = rand(ComplexF64, 4) ⋅ [x, y, z, 1]
        S = solve([f; L₁], [f; L₂], [[1, 1, 1]])
        @test nfinite(S) == 1
    end
end
