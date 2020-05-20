@testset "Valuation" begin
    @testset "Example 1" begin
        # setup
        @var x
        f = [(x - 10)^5]
        endgame_tracker, starts =
            total_degree(System(f); tracker_options = (automatic_differentiation = 3,))
        S = collect(starts)
        tracker = endgame_tracker.tracker
        val = HC2.Valuation(1)
        HC2.init!(val)
        track!(tracker, S[1], 1, 1e-13)

        t = tracker.state.t
        HC2.update!(val, tracker.predictor, real(t))

        @test val.val_x[1] ≈ 0 atol = 10 * (1e-13)^(1 / 5)
        @test val.val_tẋ[1] ≈ 1 / 5 atol = 10 * (1e-13)^(1 / 5)
    end

    @testset "Example 2" begin
        @var x y
        f = [2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3, 2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5]
        endgame_tracker, starts = total_degree(
            System(f);
            gamma = 1.3im + 0.4,
            tracker_options = TrackerOptions(
                parameters = :conservative,
                automatic_differentiation = 3,
            ),
        )
        S = collect(starts)
        tracker = endgame_tracker.tracker
        val = HC2.Valuation(2)
        tf = 1e-10
        track!(tracker, S[3], 1, tf)

        t = tracker.state.t
        HC2.update!(val, tracker.predictor, real(t))

        @test val.val_x[1] ≈ -1 atol = 10 * tf^(1 / 2)
        @test val.val_x[2] ≈ -1 atol = 10 * tf^(1 / 2)
        @test val.val_tẋ[1] ≈ -1 atol = sqrt(tf)
        @test val.val_tẋ[2] ≈ -1 atol = sqrt(tf)
    end

    @testset "Example 3" begin
        a = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32]
        @var x y
        f1 = (a[1] * x^2 + a[2] * y) * (a[3] * x + a[4] * y) + 1
        f2 = (a[1] * x^2 + a[2] * y) * (a[5] * x + a[6] * y) + 1
        endgame_tracker, starts = total_degree(
            System([f1, f2]);
            gamma = 1.3im + 0.4,
            tracker_options = TrackerOptions(
                parameters = :conservative,
                automatic_differentiation = 3,
            ),
        )
        S = collect(starts)
        tracker = endgame_tracker.tracker
        val = HC2.Valuation(2)
        tf = 1e-10
        HC2.track!(tracker, S[3], 1, tf)

        t = tracker.state.t
        HC2.update!(val, tracker.predictor, real(t))

        @test val.val_x[1] ≈ -1 / 6 atol = tf^(1 / 6)
        @test val.val_x[2] ≈ -2 / 6 atol = tf^(1 / 6)
        @test val.val_tẋ[1] ≈ -1 / 6 atol = tf^(1 / 6)
        @test val.val_tẋ[2] ≈ -2 / 6 atol = tf^(1 / 6)
        @test !isempty(sprint(show, val))
    end
end
