@testset "Valuation" begin
    @testset "Correctness (zero valuation)" begin
        @polyvar x y
        f = [x^2 - 2, x + y - 1]

        tracker, starts = coretracker_startsolutions(
            f;
            homotopy = (g, f) -> StraightLineHomotopy(g, f; gamma = 1.0 + im),
            seed = 12345,
            min_step_size = eps()^2,
            log_homotopy = true,
        )

        S = collect(starts)
        state = tracker.state
        val = HC.Valuation(S[1])
        init!(tracker, S[1], 0.0, 28.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s)
            end
        end
        @test val.ν ≈ [0.0, 0.0] atol = 1e-8
        @test norm(val.ν̇) < 1e-2
        @test norm(val.ν̈) < 1e-2

        init!(val)
        init!(tracker, S[1], 0.0, 28.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s, tracker.predictor)
            end
        end
        @test val.ν ≈ [0.0, 0.0] atol = 1e-8
        @test norm(val.ν̇) < 1e-7
        @test norm(val.ν̈) < 1e-6

        init!(val)
        init!(tracker, S[2], 0.0, 28.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s, tracker.predictor)
            end
        end
        @test val.ν ≈ [0.0, 0.0] atol = 1e-8
        @test norm(val.ν̇) < 1e-7
        @test norm(val.ν̈) < 1e-6

        ## Multi-projective tracking

        # affine valuation
        @polyvar x y v w

        tracker, starts = coretracker_startsolutions(
            [x * y - 6 * v * w, x^2 - 5 * v^2],
            variable_groups = [(x, v), (y, w)],
            log_homotopy = true,
        )
        S = collect(starts)
        state = tracker.state
        val = HC.Valuation(S[1]; affine = true)
        init!(val)
        init!(tracker, S[1], 0.0, 30.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s, tracker.predictor)
            end
        end
        @test is_success(status(tracker))
        @test val.ν ≈ [0.0, 0.0] atol = 1e-6

        val = HC.Valuation(S[1]; affine = false)
        init!(val)
        init!(tracker, S[1], 0.0, 30.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s, tracker.predictor)
            end
        end
        @test is_success(status(tracker))
        @test val.ν ≈ [0.0, 0.0, 0.0, 0.0] atol = 1e-6
    end

    @testset "Correctness (non-zero valuation)" begin
        f = equations(griewank_osborne())
        tracker, starts = coretracker_startsolutions(
            f;
            homotopy = (g, f) -> StraightLineHomotopy(g, f; gamma = 1.0 + im),
            seed = 12345,
            min_step_size = eps()^2,
            log_homotopy = true,
        )

        S = collect(starts)
        state = tracker.state
        val = HC.Valuation(S[1])
        init!(tracker, S[1], 0.0, 20.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s)
            end
        end
        @test val.ν ≈ [-0.5, -1] atol = 1e-4
        @test norm(val.ν̇) < 1e-4
        @test norm(val.ν̈) < 1e-4

        # Use analytic estimates for ν̇ and ν̈
        init!(val)
        @test all(isnan, val.ν)
        init!(tracker, S[1], 0.0, 20.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s, tracker.predictor)
            end
        end
        @test val.ν ≈ [-0.5, -1] atol = 1e-4
        @test norm(val.ν̇) < 1e-4
        @test norm(val.ν̈) < 1e-4
        @test HC.judge(val; tol = 1e-4, tol_at_infinity = 1e-4) == HC.VAL_AT_INFINITY

        init!(val)
        init!(tracker, S[2], 0.0, 25.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s, tracker.predictor)
            end
        end
        @test val.ν ≈ [1 // 3, 2 // 3] atol = 1e-3
        @test norm(val.ν̇) < 1e-3
        @test norm(val.ν̈) < 1e-3
        @test HC.judge(val; tol = 1e-3) == HC.VAL_FINITE

        init!(val)
        init!(tracker, S[2], 0.0, 25.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s)
            end
        end
        @test val.ν ≈ [1 // 3, 2 // 3] atol = 1e-3
        @test norm(val.ν̇) < 1e-3
        @test norm(val.ν̈) < 1e-3
        @test HC.judge(val; tol = 1e-3) == HC.VAL_FINITE
        @test HC.judge(val; tol = 1e-10) == HC.VAL_INDECISIVE

        init!(val)
        init!(tracker, S[3], 0.0, 25.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s)
            end
        end
        @test val.ν ≈ [2, -1] atol = 1e-5
        @test norm(val.ν̇) < 1e-5
        @test norm(val.ν̈) < 1e-5
        @test HC.judge(val; tol_at_infinity = 1e-4) == HC.VAL_AT_INFINITY
    end

    @testset "Correctness projective (non-zero valuation)" begin
        f = homogenize(equations(griewank_osborne()))
        tracker, starts = coretracker_startsolutions(
            f;
            homotopy = (g, f) -> StraightLineHomotopy(g, f; gamma = 1.0 + im),
            min_step_size = eps()^2,
            log_homotopy = true,
        )

        S = collect(starts)
        state = tracker.state
        val = HC.Valuation(state.x; affine = true)
        init!(tracker, S[1], 0.0, 25.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s)
            end
        end
        @test val.ν ≈ [-0.5, -1] atol = 1e-4
        @test norm(val.ν̇) < 1e-4
        @test norm(val.ν̈) < 1e-4
        @test HC.judge(val; tol_at_infinity = 1e-4) == HC.VAL_AT_INFINITY

        # Use analytic estimates for ν̇ and ν̈
        init!(val)
        @test all(isnan, val.ν)
        init!(tracker, S[1], 0.0, 20.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s, tracker.predictor)
            end
        end
        @test val.ν ≈ [-0.5, -1] atol = 1e-4
        @test norm(val.ν̇) < 1e-4
        @test norm(val.ν̈) < 1e-4
        @test HC.judge(val; tol_at_infinity = 1e-4) == HC.VAL_AT_INFINITY
        @test HC.judge(val; tol_at_infinity = 1e-10) == HC.VAL_INDECISIVE

        # PROJECTIVE VALUATIONS

        val = HC.Valuation(state.x; affine = false)
        init!(tracker, S[1], 0.0, 25.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s)
            end
        end
        @test val.ν ≈ [0.5, 0, 1] atol = 1e-3
        @test norm(val.ν̇) < 1e-3
        @test norm(val.ν̈) < 1e-3
        @test HC.judge(val; tol = 1e-2) == HC.VAL_FINITE

        # Use analytic estimates for ν̇ and ν̈
        init!(val)
        @test all(isnan, val.ν)
        init!(tracker, S[1], 0.0, 25.0)
        for _ in tracker
            if !state.last_step_failed
                HC.update!(val, state.x, state.ẋ, state.s, tracker.predictor)
            end
        end
        @test val.ν ≈ [0.5, 0, 1] atol = 1e-3
        @test norm(val.ν̇) < 1e-3
        @test norm(val.ν̈) < 1e-3
        @test HC.judge(val; tol = 1e-2) == HC.VAL_FINITE
    end
end
