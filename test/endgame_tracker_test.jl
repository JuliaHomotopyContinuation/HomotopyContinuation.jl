@testset "EndgameTracker" begin
    @testset "Tracking and PathResult" begin
        @var x y
        f = [2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3, 2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5]
        tracker, starts = total_degree(System(f); gamma = 1.3im + 0.4)
        S = collect(starts)
        @test !isempty(sprint(show, tracker.options))

        res = map(i -> track(tracker, S[i]; path_number = i), 1:4)

        @test is_success(res[1])
        @test path_number.(res) == 1:4
        @test start_solution.(res) == S
        @test all(isfinite, solution(res[1]))
        @test is_finite(res[1])
        @test isfinite(res[1])
        @test multiplicity(res[1]) == 1
        @test winding_number(res[1]) === nothing
        @test !is_failed(res[1])
        @test accuracy(res[1]) < 1e-12
        @test residual(res[1]) < 1e-12
        @test steps(res[1]) < 20
        @test accepted_steps(res[1]) < 20
        @test rejected_steps(res[1]) == 0
        @test is_success(res[2])
        @test is_nonsingular(res[2])
        @test cond(res[2]) < 1e3
        @test is_at_infinity(res[3])
        @test is_at_infinity(res[4])
        @test !isempty(sprint(show, res[1]))
        @test valuation(res[3]) ≈ [-1, -1] rtol = 1e-3
        # singular stuff
        @var x
        f = [(x - 10)^2]
        tracker, starts = total_degree(System(f))
        res = track.(tracker, starts)
        @test winding_number(res[1]) == 2
        @test last_path_point(res[1]) isa Tuple{Vector{ComplexF64},Float64}
        @test 0 < last(last_path_point(res[1])) < 0.1
        @test all(is_real, res)
        @test all(isreal, res)
        @test all(is_singular, res)
    end

    @testset "Overdetermined tracking" begin
        @var x y z

        p₁ = (x^2 + y^2 + z^2 - 1) * (x - 0.5)
        p₂ = (x^2 + y^2 + z^2 - 1) * (y - 0.5)
        p₃ = (z - x^2 - 2) * (x^2 + y^2 + z^2 - 1) * (z - 0.5)
        F = System([p₁, p₂, p₃])

        L₁ = rand_subspace(3; codim = 2)
        L₂ = rand_subspace(3; codim = 2)
        F_L₁ = System([x^2 + y^2 + z^2 - 1]) ∩ extrinsic(L₁)([x, y, z])
        res = track.(total_degree(F_L₁)...)
        @test all(is_success, res)

        H = IntrinsicSubspaceHomotopy(F, L₁, L₂)
        res2 = track.(EndgameTracker(H), solution.(res))
        @test all(is_success, res2)
    end
end
