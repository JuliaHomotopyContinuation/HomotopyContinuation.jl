@testset "PathTracker" begin

    @testset "EarlyAtInfinity" begin
        f = equations(griewank_osborne())

        tracker, startsolutions = pathtracker_startsolutions(f, seed=130793)
        S = collect(startsolutions)
        # This path has the special case that we track towards a non-singular solution
        # at infinity and the valuation doesn't stabilize fast enough
        result = track(tracker, S[3])
        @test result.return_code == :at_infinity
    end


    @testset "Affine tracking" begin
        @polyvar x a y b
        F = [x^2-a, x*y-a+b]
        p = [a, b]

        # parameter homotopy
        tracker = pathtracker(F, parameters=p, p₁=[1, 0], p₀=[2, 4], affine_tracking=true)
        @test tracker.problem isa Problem{AffineTracking}
        @test affine_tracking(tracker.core_tracker) == true
        @test HC.type_of_x(tracker) == Vector{ComplexF64}
        res = track(tracker, [1, 1])
        @test res.return_code == :success
        @test isa(solution(res), Vector{ComplexF64})
        @test length(solution(res)) == 2
        @test !isprojective(res)
        @test isaffine(res)

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
        @test !isprojective(res)
        @test isaffine(res)

        # total degree
        tracker, starts = pathtracker_startsolutions(equations(katsura(5)), affine_tracking=true)
        res = track(tracker, first(starts))
        @test res.return_code == :success
        @test isa(solution(res), Vector{ComplexF64})
        @test length(solution(res)) == 6
        @test !isprojective(res)
        @test isaffine(res)
    end

    @testset "Details Level" begin
        @polyvar x a y b
        F = [x^2-a, x*y-a+b]
        p = [a, b]

        tracker = pathtracker(F, parameters=p, p₁=[1, 0], p₀=[2, 4], affine_tracking=true)

        res = track(tracker, [1, 1]; details=:minimal)
        test_show_juno(res)
        @test res.condition_jacobian === nothing

        res = track(tracker, [1, 1]; details=:default)
        @test res.condition_jacobian !== nothing

        res = track(tracker, [1, 1]; details=:extensive)
        @test res.valuation !== nothing

        @test isnonsingular(res)
    end

    @testset "Parameter homotopy with AbstractSystem" begin
        # affine
        @polyvar x a y b
        F = SPSystem([x^2-a, x*y-a+b]; parameters=[a, b])

        tracker, starts = pathtracker_startsolutions(F, [1.0, 1.0 + 0.0*im],
                        p₁=[1, 0], p₀=[2, 4], affine_tracking=true)
        res = track(tracker, starts[1])
        @test res.return_code == :success

        # affine without explicit flag (numerical degree check fails)
        tracker, starts = pathtracker_startsolutions(F, [1.0, 1.0 + 0.0*im], p₁=[1, 0], p₀=[2, 4])
        res = track(tracker, starts[1])
        @test res.return_code == :success

        # projective
        @polyvar x y z a b
        F = SPSystem([x^2-a*z^2, x*y-a*z^2+b*z^2]; parameters=[a, b])

        tracker, starts = pathtracker_startsolutions(F, [1.0, 1.0 + 0.0*im], p₁=[1, 0], p₀=[2, 4], affine_tracking=false)
        res = track(tracker, starts[1])
        @test res.return_code == :success
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
        result_p₀ = solve(F_p₀)

        # Let's store the generic solutions
        S_p₀ = solutions(result_p₀)
        # Construct the PathTracker
        tracker = pathtracker(F; parameters=p, generic_parameters=p₀)

        q = randn(3)
        result = track(tracker, first(S_p₀); target_parameters=q)
        @test issuccess(result)

        result = track(tracker, first(S_p₀); start_parameters=p₀, target_parameters=q)
        @test issuccess(result)

        # Construct the PathTracker
        tracker = pathtracker(F; parameters=p, generic_parameters=p₀)
        set_parameters!(tracker; target_parameters=q)
        result2 = track(tracker, first(S_p₀))
        @test result.solution ≈ result2.solution atol=1e-7
    end
end
