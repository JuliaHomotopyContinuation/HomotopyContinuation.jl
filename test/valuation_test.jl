@testset "Valuation" begin
    @testset "Example 1" begin
        # setup
        @var x
        f = [(x - 10)^5]
        H, starts = total_degree_homotopy(f, [x])
        S = collect(starts)
        tracker = Tracker(H)
        val = HC2.Valuation(1)
        HC2.init!(val)
        track!(tracker, S[1], 1, 1e-13)

        @unpack x, t, x¹, x², x³ = tracker.state
        HC2.update!(val, x, x¹, x², x³, real(t))

        @test val.val_x[1] ≈ 0 atol = 10(1e-13)^(1/5)
        @test val.val_tx¹[1] ≈ 1/5 atol = 10(1e-13)^(1/5)
    end

    @testset "Example 2" begin
        @var x y
        f = [
            2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3,
            2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5,
        ]
        H, starts = total_degree_homotopy(f, [x, y]; gamma = 1.3im + 0.4)
        S = collect(starts)
        tracker = Tracker(H, β_τ = 0.5)
        val = HC2.Valuation(2)
        tf = 1e-10
        track!(tracker, S[3], 1, tf)

        @unpack x, t, x¹, x², x³ = tracker.state
        HC2.update!(val, x, x¹, x², x³, real(t))

        @test val.val_x[1] ≈ -1 atol = tf^(1/2)
        @test val.val_x[2] ≈ -1 atol = tf^(1/2)
        @test val.val_tx¹[1] ≈ -1 atol = tf
        @test val.val_tx¹[2] ≈ -1 atol = tf
    end

    @testset "Example 3" begin
        a = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32]
        @var x y
        f1 = (a[1] * x^2 + a[2] * y) * (a[3] * x + a[4] * y) + 1
        f2 = (a[1] * x^2 + a[2] * y) * (a[5] * x + a[6] * y) + 1
        H, starts = total_degree_homotopy([f1, f2], [x, y]; gamma = 1.3im + 0.4)
        S = collect(starts)
        tracker = Tracker(H, β_τ = 0.5)
        val = HC2.Valuation(2)
        tf = 1e-10
        HC2.track!(tracker, S[3], 1, tf)

        @unpack x, t, x¹, x², x³ = tracker.state
        HC2.update!(val, x, x¹, x², x³, real(t))

        @test val.val_x[1] ≈ -1/6 atol = tf^(1/6)
        @test val.val_x[2] ≈ -2/6 atol = tf^(1/6)
        @test val.val_tx¹[1] ≈ -1/6 atol = tf^(1/6)
        @test val.val_tx¹[2] ≈ -2/6 atol = tf^(1/6)
        @test !isempty(sprint(show, val))
    end
end
