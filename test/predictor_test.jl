@testset "Extended predictor solves" begin
    @var x y p q
    F = System(
        [x + y - 1, x + y - 1 + p * ((1 + q) * y - 1)];
        variables = [x, y],
        parameters = [p, q],
    )
    for compile in (false, true, :mixed)
        H = ParameterHomotopy(F, [2.0^-60, 1], [2.0^-60, 0]; compile)
        P = HC.Predictor(H)
        J = HC.Jacobian(zeros(ComplexF64, 2, 2))
        z = ComplexF64[1/3, 2/3]
        r = zeros(ComplexF64, 2)
        evaluate_and_jacobian!(r, HC.matrix(J), H, z, 0.5)
        HC.updated!(J)
        @test P.extended_workspace === nothing
        HC.update!(P, H, z, 0.5, J, InfNorm())
        @test P.extended_workspace !== nothing
        @test P.used_extended_precision
        @test P.tx³[2, 2] ≈ -4 / 9
        @test P.tx³[2, 3] ≈ 8 / 27
        @test P.tx³[2, 4] ≈ -16 / 81
        @test P.tx³[1, 2] ≈ 4 / 9
        @test all(isfinite, P.tx_norm)
        cache = P.extended_workspace
        @test HC.extended_predictor!(P, H, z, 0.5, 3)
        @test P.extended_workspace === cache

        disabled = HC.Predictor(H)
        HC.update!(disabled, H, z, 0.5, J, InfNorm(); extended_precision = false)
        @test disabled.extended_workspace === nothing
        @test !disabled.used_extended_precision
        tracker = Tracker(H)
        result = track(tracker, [0.5, 0.5])
        @test is_success(result)
        @test solution(result) ≈ [0, 1] atol = 1e-12
        @test tracker.state.used_extended_prec

        # Returning to a well-conditioned problem uses ordinary arithmetic;
        # retaining an allocated workspace must not count as using it.
        HC.parameters!(H, [1, 1], [1, 0])
        evaluate_and_jacobian!(r, HC.matrix(J), H, z, 0.5)
        HC.updated!(J)
        HC.update!(P, H, z, 0.5, J, InfNorm())
        @test !P.used_extended_precision
        @test P.extended_workspace === cache
        @test P.tx³[2, 3] ≈ 8 / 27
    end

    # A failed extended solve must leave the previous coefficients intact.
    H = InterpretedHomotopy(Homotopy([x + y - q, x + y - q], [x, y], q))
    P = HC.Predictor(H)
    fill!(P.tx³.data, 7)
    @test !HC.extended_predictor!(P, H, ComplexF64[0, 0], 0.5, 3)
    @test all(==(7), P.tx³.data)
end
