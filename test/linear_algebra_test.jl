@testset "LinearAlgebra" begin
    @testset "Constructor + basics" begin
        A = rand(ComplexF32, 12, 12)
        WS = HC2.MatrixWorkspace(A)
        @test WS isa HC2.MatrixWorkspace

        WS.factorized[] = true
        B = randn(ComplexF64, 12, 12)
        HC2.update!(WS, B)
        @test WS.A == B
        @test WS.A !== B
        @test WS.factorized[] == false

        # test indexing etc
        @test WS[2, 1] == B[2, 1]
        @test WS[2] == B[2]
        WS[2, 3] = 5.0
        @test WS[2, 3] == 5.0
        WS[5] = 2.32
        @test WS[5] == 2.32
        @test size(WS) == (12, 12)
        C = randn(ComplexF64, 12, 12)
        WS .= C
        @test WS.A == C

        #test struct array switch
        A = rand(ComplexF32, 30, 30)
        WS = HC2.MatrixWorkspace(A)
    end

    @testset "QR" begin
        A = randn(ComplexF64, 12, 7)
        MW = HC2.MatrixWorkspace(A)
        B = randn(ComplexF64, 12, 7)
        true_qr = LinearAlgebra.qrfactUnblocked!(copy(B))
        MW .= B

        HC2.updated!(MW)

        b = randn(ComplexF64, 12)
        b2 = copy(b)
        x = zeros(ComplexF64, 7)

        ldiv!(x, MW, b)

        @test norm(MW.qr.factors - true_qr.factors) /norm(MW.qr.factors) < 1e-14
        @test norm(b2 - b) == 0
        @test norm(x - (B \ b2)) / norm(x) < 1e-14
    end

    @testset "ldiv" begin
        for n in [3, 13, 31] # test that struct array also works
            A = randn(ComplexF64, n, n)
            b = randn(ComplexF64, n)
            x = zeros(ComplexF64, n)
            WS = HC2.MatrixWorkspace(A)
            HC2.update!(WS, A)
            ldiv!(x, WS, b)
            WS.factorized[] = false
            @test (@allocated ldiv!(x, WS, b)) == 0
            @test (lu(A) \ b) ≈ x rtol = cond(A) * 10
        end
    end

    @testset "Inf-norm estimator / cond" begin
        d_r = rand() * exp10.(range(-6, stop = 6, length = 6))
        D_R = diagm(d_r)
        d_l = rand() * exp10.(range(6, stop = -6, length = 6))
        D_L = diagm(d_l)

        A = randn(6, 6)

        B = A * inv(D_R)
        WB = HC2.MatrixWorkspace(B)
        @test 0.1 ≤ opnorm(inv(B), Inf) / HC2.inverse_inf_norm_est(WB) ≤ 10
        @test 0.1 ≤
              opnorm(inv(B * D_R), Inf) / HC2.inverse_inf_norm_est(WB, nothing, d_r) ≤
              10
        @test 0.1 ≤ opnorm(inv(D_R * B), Inf) / HC2.inverse_inf_norm_est(WB, d_r) ≤ 10

        @test 0.1 ≤ cond(B, Inf) / cond(WB) ≤ 10
        @test 0.1 ≤ cond(B * D_R, Inf) / cond(WB, nothing, d_r) ≤ 10
        @test 0.1 ≤ cond(D_R * B, Inf) / cond(WB, d_r, nothing) ≤ 10


        C = inv(D_L) * A
        WC = HC2.MatrixWorkspace(C)
        @test 0.1 ≤ opnorm(inv(C), Inf) / HC2.inverse_inf_norm_est(WC) ≤ 10
        @test 0.1 ≤ opnorm(inv(D_L * C), Inf) / HC2.inverse_inf_norm_est(WC, d_l) ≤ 10
        @test 0.1 ≤
              opnorm(inv(C * D_L), Inf) / HC2.inverse_inf_norm_est(WC, nothing, d_l) ≤
              10
        @test 0.1 ≤ cond(C, Inf) / cond(WC) ≤ 10
        @test 0.1 ≤ cond(C * D_L, Inf) / cond(WC, nothing, d_l) ≤ 10
        @test 0.1 ≤ cond(D_L * C, Inf) / cond(WC, d_l, nothing) ≤ 10


        D = inv(D_L) * A * inv(D_R)
        WD = HC2.MatrixWorkspace(D)
        @test 0.1 ≤ opnorm(inv(D), Inf) / HC2.inverse_inf_norm_est(WD) ≤ 10
        @test 0.1 ≤
              opnorm(inv(D_L * D * D_R), Inf) / HC2.inverse_inf_norm_est(WD, d_l, d_r) ≤
              10
        @test HC2.inverse_inf_norm_est(WD, d_r, d_l) >
              100 * HC2.inverse_inf_norm_est(WD, d_l, d_r)
        @test 0.1 ≤ cond(D, Inf) / cond(WD) ≤ 10
        @test 0.1 ≤ cond(D_L * D * D_R, Inf) / cond(WD, d_l, d_r) ≤ 10
    end

    @testset "Jacobian" begin
        A = randn(ComplexF64, 6, 6)
        b = randn(ComplexF64, 6)
        x̂ = similar(b)
        x = similar(b)
        JM = HC2.Jacobian(zeros(ComplexF64, 6, 6))
        @test JM isa HC2.Jacobian
        HC2.matrix(JM) .= A
        HC2.updated!(JM)
        ldiv!(x̂, JM, b)
        ldiv!(x, HC2.workspace(JM), b)
        @test x == x̂
    end
end
