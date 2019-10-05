@testset "CoreTracker" begin
    @testset "General" begin
        F = equations(katsura(5))
        # test construction
        t1, start_sols = coretracker_startsolutions(F)

        test_show_juno(t1)
        @test_nowarn show(devnull, t1)
        @test !isempty(string(t1))
        test_show_juno(options(t1))
        @test_nowarn show(devnull, options(t1))
        test_show_juno(t1.state)
        @test_nowarn show(devnull, t1.state)
        @test !isempty(string(t1.state))

        @test_nowarn current_Δt(t1)

        caccuracy = accuracy(t1)
        @test_nowarn set_accuracy!(t1, 5e-5)
        @test accuracy(t1) == 5e-5
        set_accuracy!(t1, caccuracy)

        cmaxiter = max_corrector_iters(t1)
        @test_nowarn set_max_corrector_iters!(t1, 11)
        @test max_corrector_iters(t1) == 11
        set_max_corrector_iters!(t1, cmaxiter)

        @test t1 isa CoreTracker
        @test isa(current_x(t1), Vector)
        @test length(current_x(t1)) == 6

        t3, start_sols = coretracker_startsolutions(F, predictor = Euler())
        @test t3.predictor isa HC.EulerCache

        @test_deprecated setup!(t1, first(start_sols), 1.0, 0.4)
        init!(t1, first(start_sols), 1.0, 0.4)
        @test is_tracking(status(t1))
        @test current_t(t1) == 1.0

        setup!(t1, first(start_sols), 0.5, 0.4)
        @test is_invalid_startvalue(status(t1))
        @test is_terminated(status(t1))
        @test current_t(t1) == 0.5
        @test real(current_Δt(t1)) < 0
        @test cond(t1) ≥ 1

        R = track(t1, first(start_sols), 1.0, 0.0)
        @test R isa CoreTrackerResult
        @test is_success(R.returncode)
        @test R.accuracy < 1e-7
        @test_nowarn show(devnull, R)
        @test !isempty(string(R))
        test_show_juno(R)

        out = copy(first(start_sols))
        retcode = track!(out, t1, first(start_sols), 1.0, 0.0)
        @test is_success(retcode)
        @test out == solution(R)


        @test norm(t1) isa AbstractNorm
    end

    @testset "log_homotopy" begin
        F = equations(katsura(5))
        log_tracker, start_sols = coretracker_startsolutions(F; log_homotopy = true)
        @test log_tracker.homotopy.homotopy isa HC.LogHomotopy
        @test log_tracker.options.logarithmic_time_scale
        @test is_success(track!(log_tracker, first(start_sols), 0.0, 32.0))
    end

    @testset "deprecated" begin
        F = equations(katsura(5))
        # test construction
        t1 = coretracker(F)

        @test_warn(_ -> true, refinement_accuracy(t1))
        @test_warn(_ -> true, set_refinement_accuracy!(t1, 1e-4))
        @test_warn(_ -> true, max_refinement_iters(t1))
        @test_warn(_ -> true, set_max_refinement_iters!(t1, 11))
    end

    @testset "Affine tracking" begin
        @polyvar x a y b
        F = [x^2 - a, x * y - a + b]
        p = [a, b]

        tracker, starts = coretracker_startsolutions(
            F,
            [1.0, 1.0 + 0.0 * im],
            parameters = p,
            p₁ = [1, 0],
            p₀ = [2, 4],
            affine_tracking = true,
        )
        @test affine_tracking(tracker) == true
        res = track(tracker, starts[1], 1.0, 0.0)
        @test is_success(res)
        @test isa(res.x, Vector{ComplexF64})
        @test length(res.x) == 2

        x = @SVector [1.0, 1.0 + 0.0 * im]
        tracker, starts = coretracker_startsolutions(
            F,
            x,
            parameters = p,
            p₁ = [1, 0],
            p₀ = [2, 4],
            affine_tracking = true,
        )
        @test length(starts) == 1
        res = track(tracker, starts[1], 1.0, 0.0)
        @test isa(res.x, Vector{ComplexF64})
        @test length(res.x) == 2
    end

    @testset "Allocations" begin
        F = equations(cyclic(5))
        tracker, start_sols = coretracker_startsolutions(F; seed = 12356)
        s = first(start_sols)
        @allocated track!(tracker, s, 1.0, 0.01)
        # The path tracker should not allocate after the first tracked path
        @test (@allocated track!(tracker, s, 1.0, 0.01)) == 0

        # log homotopy as well
        log_tracker, log_start_sols = coretracker_startsolutions(
            F;
            seed = 12356,
            log_homotopy = true,
        )
        log_s = first(log_start_sols)
        @allocated track!(log_tracker, log_s, 0.0, -log(0.01))
        # The path tracker should not allocate after the first tracked path
        @test (@allocated track!(log_tracker, log_s, 0.0, -log(0.01))) == 0
    end

    @testset "LinearSystem" begin
        A = rand(3, 3)
        b = rand(3)

        @polyvar x[1:3]
        F = A * x - b

        tracker, start_sols = coretracker_startsolutions(F; max_steps = 10)
        s = first(start_sols)
        result = track(tracker, s, 1.0, 0.0)
        @test is_success(result)
        @test result.t == 0.0
        x = result.x
        @test norm(x - A \ b) < 1e-6

        x_inter = copy(s)
        retcode = track!(x_inter, tracker, s, 1.0, 0.1)
        @test is_success(retcode)
        x_final = zero(x_inter)
        retcode = track!(x_final, tracker, x_inter, 0.1, 0.0)
        @test is_success(retcode)
        @test steps(tracker) < 4
        x = current_x(tracker)
        @test norm(x - A \ b) < 1e-6
    end

    @testset "FixedPatch continuation" begin
        F = homogenize(equations(katsura(5)))
        # test construction
        tracker, start_sols = coretracker_startsolutions(F, patch = OrthogonalPatch())
        patch = tracker.state.patch
        fixedtracker = coretracker(F, patch = FixedPatch())
        fixedpatch = fixedtracker.state.patch
        x1 = first(start_sols)

        r1 = track(fixedtracker, x1, 1.0, 0.1)
        v1 = copy(fixedpatch.v̄)
        track!(fixedtracker, r1.x, 0.1, 0.05, setup_patch = false)
        @test v1 == fixedpatch.v̄
        track!(fixedtracker, current_x(tracker), 0.05, 0.01, setup_patch = false)
        @test v1 == fixedpatch.v̄
        track!(fixedtracker, r1.x, 0.1, 0.05, setup_patch = true)
        @test v1 != fixedpatch.v̄

        r1 = track(tracker, x1, 1.0, 0.1)
        v1 = copy(patch.v̄)
        track!(tracker, r1.x, 0.1, 0.0)
        @test v1 != patch.v̄
    end

    @testset "coretracker_startsolutions" begin
        @polyvar x y
        tracker, S = coretracker_startsolutions([x^2 + y^2 - 4, x + y - 3])
        result = track(tracker, first(S), 1.0, 0.0)
        @test is_success(result)

        # test type promotion of startsolutions
        tracker, S = coretracker_startsolutions([x^2 + y^2 - 4, x + y - 3])
        result = track(tracker, [1, 1], 1.0, 0.0)
        @test is_success(result)

        tracker, S = coretracker_startsolutions(
            [x^2 + y^2 - 4, x + y - 3];
            simple_step_size_alg = true,
        )
        result = track(tracker, first(S), 1.0, 0.0)
        @test is_success(result)
    end

    @testset "iterator" begin
        @polyvar x a y b
        F = [x^2 - a, x * y - a + b]
        p = [a, b]
        s = [1.0, 1.0 + 0.0 * im]
        tracker = coretracker(F, parameters = p, p₁ = [1, 0], p₀ = [2, 4])
        init!(tracker, s, 1.0, 0.0)
        for tr in tracker
            # nothing
        end
        @test is_success(status(tracker))

        # path iterator
        typeof(first(iterator(tracker, s, 1.0, 0.0))) == Tuple{Vector{ComplexF64},Float64}

        @test max_step_size(tracker) == Inf
        set_max_step_size!(tracker, 0.01)
        @test max_step_size(tracker) == 0.01

        @test length(collect(iterator(tracker, s, 1.0, 0.0))) == 101

        @polyvar x a
        f = [x - a]
        ct = coretracker(
            f,
            [1.0],
            parameters = [a],
            p₁ = [1.0],
            p₀ = [2.0],
            max_step_size = 0.1,
        )
        Xs = Vector{ComplexF64}[]
        for (x, t) in iterator(ct, [1.0], 1.0, 0.0)
            push!(Xs, x)
        end
        @test round.(Int, real.(first.(Xs)) .* 10) == collect(10:20)
    end

    @testset "Limit Accuracy / Adaptive Precision" begin
        # We use Wilkinson polynomials of degree 12 and 15
        # 12 is fine for 64 bit precision, but in 15 is the evaluation error too large
        @polyvar x
        n = 12
        g = [x^n - 1]
        f = [prod(x - i for i = 1:n)]
        S = [[cis(i * 2π / n)] for i = 0:(n-1)]

        tracker = coretracker(
            g,
            f,
            S;
            log_homotopy = true,
            min_step_size = eps()^2,
            accuracy = 1e-7,
            precision = PRECISION_FIXED_64,
        )
        results = map(s -> track(tracker, s, 0.0, 70), S)
        @test all(is_success, results)

        n = 15
        g = [x^n - 1]
        f = [prod(x - i for i = 1:n)]
        S = [[cis(i * 2π / n)] for i = 0:(n-1)]

        tracker = coretracker(
            g,
            f,
            S;
            log_homotopy = true,
            min_step_size = eps()^2,
            accuracy = 1e-7,
            precision = :double,
        )
        results = map(s -> track(tracker, s, 0.0, 70), S)
        @test all(
            r -> is_success(r) || r.returncode == HC.CoreTrackerStatus.terminated_accuracy_limit,
            results,
        )


        tracker = coretracker(
            g,
            f,
            S;
            log_homotopy = true,
            min_step_size = eps()^2,
            accuracy = 1e-7,
            precision = :double_double,
        )
        results = map(s -> track(tracker, s, 0.0, 70), S)
        @test all(is_success, results)

        tracker = coretracker(
            g,
            f,
            S;
            log_homotopy = true,
            min_step_size = eps()^2,
            accuracy = 1e-7,
            precision = :adaptive,
        )
        results = map(s -> track(tracker, s, 0.0, 70), S)
        @test all(is_success, results)


        @test_throws ArgumentError coretracker(
            g,
            f,
            S;
            log_homotopy = true,
            precision = :NOT_A_PRECISION,
        )
    end

    @testset "Change parameters" begin
        @polyvar x a y b
        F = SPSystem([x^2 - a, x * y - a + b]; parameters = [a, b])

        tracker, starts = coretracker_startsolutions(
            F,
            [1.0, 1.0 + 0.0 * im],
            generic_parameters = [2.2, 3.2],
        )
        start_parameters!(tracker, [1, 0])
        target_parameters!(tracker, [2, 4])
        @test is_success(track(tracker, starts[1], 1.0, 0.0))
    end

    @testset "(Multi-) projective" begin
        @polyvar x y u v
        f = [x * y - 6, x^2 - 5]
        tracker, starts = coretracker_startsolutions(
            f,
            variable_groups = [[x], [y]],
            seed = 123456,
            predictor = Pade21(),
        )
        current_x(tracker) isa ProjectiveVectors.PVector{ComplexF64,2}
        S = collect(starts)
        R1 = track(tracker, S[1], 1.0, 0.0)
        @test is_success(R1)
        @test solution(R1) isa ProjectiveVectors.PVector{ComplexF64,2}

        R2 = track(tracker, S[2], 1.0, 0.0)
        @test is_success(R2)
    end
end
