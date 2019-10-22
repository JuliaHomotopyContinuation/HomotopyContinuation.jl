@testset "Cauchy Endgame" begin
    # Case: 1 regular solution
    @testset "Regular solution" begin
        @polyvar x y
        f = [x^2 - 4, x + 2y - 5]
        tracker, S = coretracker_startsolutions(f; seed = 1234, log_homotopy = true)
        track!(tracker, first(S), 0.0, 20.0)
        eg = HC.CauchyEndgame(tracker.state.x)
        p = copy(tracker.state.x)
        code, m = HC.predict!(p, tracker, eg)
        @test code == HC.CAUCHY_SUCCESS
        @test m == 1
        @test p ≈ [2, 1.5] atol = 1e-8
    end

    # Case 2: Univariate Singular
    @testset "Univariate singular" begin
        @polyvar x
        f = [(x - 3)^5]
        tracker, S = coretracker_startsolutions(f; seed = 1234, log_homotopy = true)
        track!(tracker, first(S), 0.0, 20.0)
        s = copy(tracker.state.x)
        eg = HC.CauchyEndgame(s)
        p = copy(tracker.state.x)
        code, m = HC.predict!(p, tracker, eg)
        @test code == HC.CAUCHY_SUCCESS
        @test m == 5
        @test p ≈ [3] atol = 1e-8
        # test that we can continute tracking from there
        track!(tracker, s, 20.0, 23.0)
        p .= 0
        code, m = HC.predict!(p, tracker, eg)
        @test code == HC.CAUCHY_SUCCESS
        @test m == 5
        @test p ≈ [3] atol = 1e-8

        # use in iterator
        f = [(x - 3)^5]
        tracker, S = coretracker_startsolutions(f; seed = 1234, log_homotopy = true)
        init!(tracker, first(S), 0.0, 25.0)
        # use as an iterator to do cauchy endgame in between,
        # make sure that we can still continue the tracking
        for _ in tracker
            if tracker.state.s > 14
                code, m = HC.predict!(p, tracker, eg)
                @test code == HC.CAUCHY_SUCCESS
                @test m == 5
            end
        end
        @test is_success(tracker.state.status)
        @test current_t(tracker) == 25
    end

    # Case 3
    @testset "Griewank-Osborne" begin
        # affine
        f = equations(griewank_osborne())
        tracker, starts = coretracker_startsolutions(f; seed = 1234, log_homotopy = true)
        S = collect(starts)
        track!(tracker, S[3], 0.0, 20.0)
        s = copy(tracker.state.x)
        eg = HC.CauchyEndgame(s)
        p = copy(tracker.state.x)
        code, m = HC.predict!(p, tracker, eg)
        @test p ≈ [0, 0] atol = 1e-8
        @test m == 3

        # projective
        f = equations(griewank_osborne())
        # random patch
        tracker, starts = coretracker_startsolutions(
            f;
            seed = 1234,
            projective_tracking = true,
            log_homotopy = true,
            patch = RandomPatch(),
        )
        S = collect(starts)
        track!(tracker, S[3], 0.0, 20.0)
        s = copy(tracker.state.x)
        eg = HC.CauchyEndgame(s)
        p = copy(tracker.state.x)
        code, m = HC.predict!(p, tracker, eg)
        @test ProjectiveVectors.affine_chart(p) ≈ [0, 0] atol = 1e-8
        @test m == 3

        # OrthogonalPatch
        tracker, starts = coretracker_startsolutions(
            f;
            seed = 1234,
            projective_tracking = true,
            log_homotopy = true,
            patch = OrthogonalPatch(),
        )
        S = collect(starts)
        track!(tracker, S[3], 0.0, 20.0)
        s = copy(tracker.state.x)
        eg = HC.CauchyEndgame(s)
        p = copy(tracker.state.x)
        code, m = HC.predict!(p, tracker, eg)
        @test ProjectiveVectors.affine_chart(p) ≈ [0, 0] atol = 1e-8
        @test m == 3
    end

    @testset "Limit accuracy error" begin
        # affine
        @polyvar x
        f = [(x - 3)^5]
        tracker, starts = coretracker_startsolutions(
            f;
            seed = 1234,
            min_step_size = eps()^2,
            log_homotopy = true,
        )
        S = collect(starts)
        tracker.options.precision = PRECISION_ADAPTIVE
        track(tracker, S[3], 0.0, 30.0)
        tracker.options.precision = PRECISION_FIXED_64
        s = copy(tracker.state.x)
        eg = HC.CauchyEndgame(s)
        p = copy(tracker.state.x)
        code, m = HC.predict!(p, tracker, eg)
        @test code == HC.CAUCHY_TERMINATED_ACCURACY_LIMIT
    end
end
