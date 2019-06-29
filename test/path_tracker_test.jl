@testset "PathTracker" begin

    @testset "Regauging" begin
        A = rand(3, 3)
        b = rand(3)

        @polyvar x[1:3]
        F = A * x - b

        tracker, start_sols = pathtracker_startsolutions(F, max_steps=10)
        s = first(start_sols)
        result = track(tracker, s)
        @test result.return_code == :success
        @test result.solution ≈ A \ b atol=1e-6

        @polyvar x y
        result = solve([x^2+y^2-1, 3x+2y-1], threading=false)
        u, v = solutions(result)
        @test norm(u) ≈ 1.0 atol=1e-7
        @test norm(v) ≈ 1.0 atol=1e-7
    end

    @testset "Affine tracking" begin
        @polyvar x a y b
        F = [x^2-a, x*y-a+b]
        p = [a, b]

        # parameter homotopy
        tracker = pathtracker(F, parameters=p, p₁=[1, 0], p₀=[2, 4], projective_tracking=true)
        @test tracker.problem isa Problem{ProjectiveTracking}

        tracker = pathtracker(F, parameters=p, p₁=[1, 0], p₀=[2, 4], projective_tracking=false)
        @test tracker.problem isa Problem{AffineTracking}

        tracker = pathtracker(F, parameters=p, p₁=[1, 0], p₀=[2, 4], projective_tracking=true, affine_tracking=true)
        @test tracker.problem isa Problem{AffineTracking}

        tracker = pathtracker(F, parameters=p, p₁=[1, 0], p₀=[2, 4], projective_tracking=false, affine_tracking=false)
        @test tracker.problem isa Problem{ProjectiveTracking}

        tracker = pathtracker(F, parameters=p, p₁=[1, 0], p₀=[2, 4], affine_tracking=false)
        @test tracker.problem isa Problem{ProjectiveTracking}

        tracker = pathtracker(F, parameters=p, p₁=[1, 0], p₀=[2, 4], affine_tracking=true)
        @test tracker.problem isa Problem{AffineTracking}

        @test affine_tracking(tracker.core_tracker) == true
        @test HC.type_of_x(tracker) == Vector{ComplexF64}
        res = track(tracker, [1, 1])
        @test res.return_code == :success
        @test isa(solution(res), Vector{ComplexF64})
        @test length(solution(res)) == 2
        @test !is_projective(res)
        @test is_affine(res)

        x = @SVector [1.0, 1.0 + 0.0*im]
        tracker, starts = pathtracker_startsolutions(F, x;
                                parameters=p, p₁=[1, 0], p₀=[2, 4], affine_tracking=true)
        @test tracker.problem isa Problem{AffineTracking}
        @test length(starts) == 1
        res = track(tracker, starts[1])
        @test isa(solution(res), Vector{ComplexF64})
        @test length(solution(res)) == 2

        # start target homotopy
        F₁ = subs.(F, Ref([a, b] => [1, 0]))
        F₀ = subs.(F, Ref([a, b] => [2, 4]))
        tracker = pathtracker(F₁, F₀, affine_tracking=true)
        @test tracker.problem isa Problem{AffineTracking}
        res = track(tracker, [1, 1])
        @test res.return_code == :success
        @test isa(solution(res), Vector{ComplexF64})
        @test length(solution(res)) == 2
        @test !is_projective(res)
        @test is_affine(res)

        # total degree
        tracker, starts = pathtracker_startsolutions(equations(katsura(5)), affine_tracking=true)
        res = track(tracker, first(starts))
        @test res.return_code == :success
        @test isa(solution(res), Vector{ComplexF64})
        @test length(solution(res)) == 6
        @test !is_projective(res)
        @test is_affine(res)
    end

    @testset "Details Level" begin
        @polyvar x a y b
        F = [x^2-a, x*y-a+b]
        p = [a, b]

        tracker = pathtracker(F, parameters=p, p₁=[1, 0], p₀=[2, 4])

        res = track(tracker, [1, 1]; details=:minimal)
        test_show_juno(res)
        @test res.residual === nothing

        res = track(tracker, [1, 1]; details=:default)
        @test res.residual !== nothing

        res = track(tracker, [1, 1]; details=:extensive)
        @test res.valuation !== nothing

        @test is_nonsingular(res)
    end

    @testset "Parameter homotopy with AbstractSystem" begin
        # affine
        @polyvar x a y b
        F = SPSystem([x^2-a, x*y-a+b]; parameters=[a, b])

        tracker, starts = pathtracker_startsolutions(F, [1.0, 1.0 + 0.0*im],
                        p₁=[1, 0], p₀=[2, 4])
        res = track(tracker, starts[1])
        @test res.return_code == :success

        # affine without explicit flag (numerical degree check fails)
        tracker, starts = pathtracker_startsolutions(F, [1.0, 1.0 + 0.0*im], p₁=[1, 0], p₀=[2, 4])
        res = track(tracker, starts[1])
        @test res.return_code == :success

        # projective
        @polyvar x y z a b
        F = SPSystem([x^2-a*z^2, x*y-a*z^2+b*z^2]; parameters=[a, b])

        tracker, starts = pathtracker_startsolutions(F, [1.0, 1.0 + 0.0*im],
                    p₁=[1, 0], p₀=[2, 4])
        res = track(tracker, starts[1])
        @test res.return_code == :success
        @test isa(res.solution, ProjectiveVectors.PVector)

        # projective with homvar
        tracker, starts = pathtracker_startsolutions(F, [1.0, 1.0 + 0.0*im],
                    p₁=[1, 0], p₀=[2, 4], homvar=z)
        @test isa(tracker.core_tracker.affine_patch, EmbeddingPatch)
        res = track(tracker, starts[1])
        @test res.return_code == :success
        @test isa(res.solution, Vector)
    end

    @testset "Parameter homotopy, change parameters" begin
        @polyvar x y z p[1:3]

        F = [
            x + 3 + 2y + 2y^2 - p[1],
            (x - 2 + 5y)*z + 4 - p[2] * z,
            (x + 2 + 4y)*z + 5 - p[3] * z
        ]
        # Generate generic parameters by sampling complex numbers from the normal distribution
        p₀ = randn(ComplexF64, 3)
        # Substitute p₀ for p
        F_p₀ = subs(F, p => p₀)
        # Compute all solutions for F_p₀
        result_p₀ = solve(F_p₀; affine_tracking=true)

        # Let's store the generic solutions
        S_p₀ = solutions(nonsingular(result_p₀))
        # Construct the PathTracker
        tracker = pathtracker(F; parameters=p, generic_parameters=p₀)

        q = randn(3)
        result = track(tracker, first(S_p₀); target_parameters=q)
        @test is_success(result)

        result = track(tracker, first(S_p₀); start_parameters=p₀, target_parameters=q)
        @test is_success(result)

        # Construct the PathTracker
        tracker = pathtracker(F; parameters=p, generic_parameters=p₀)
        set_parameters!(tracker; target_parameters=q)
        result2 = track(tracker, first(S_p₀))
        @test result.solution ≈ result2.solution atol=1e-7
    end

    @testset "Real Jacobian" begin
        @polyvar v w
        @test nfinite(solve([v - 2w , -v + 2w, v + w - 3])) == 1
        @test nfinite(solve([v - 2w , -v + 2w, v + w - 3], affine_tracking=false)) == 1
    end

    @testset "accuracy assigned" begin
        @polyvar x y
        f₁ = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        f₂ = x^2+2x*y^2 - 2y^2 - 1/2
        result = solve([f₁, f₂])
        acc = results(result; only_real=true)[1].accuracy
        @test acc !== nothing
        @test !isnan(acc)
    end
end
