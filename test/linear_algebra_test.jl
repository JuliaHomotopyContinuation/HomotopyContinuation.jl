@testset "LinearAlgebra" begin
    include("linear_algebra_triangular_test.jl")
    @testset "Constructor + basics" begin
        A = rand(ComplexF32, 12, 12)
        WS = HC.MatrixWorkspace(A)
        @test WS isa HC.MatrixWorkspace

        WS.factorized[] = true
        B = randn(ComplexF64, 12, 12)
        HC.update!(WS, B)
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
        WS = HC.MatrixWorkspace(A)
    end

    @testset "QR" begin
        A = randn(ComplexF64, 12, 7)
        MW = HC.MatrixWorkspace(A)
        B = randn(ComplexF64, 12, 7)
        true_qr = LinearAlgebra.qrfactUnblocked!(copy(B))
        MW .= B

        HC.updated!(MW)

        b = randn(ComplexF64, 12)
        b2 = copy(b)
        x = zeros(ComplexF64, 7)

        ldiv!(x, MW, b)

        # the QR factors are stored as separate real and imaginary parts and only
        # assembled on demand
        QR = HC.synced_qr(MW)
        @test norm(QR.factors - true_qr.factors) / norm(QR.factors) < 1e-14
        @test norm(QR.τ - true_qr.τ) / norm(QR.τ) < 1e-14
        @test norm(b2 - b) == 0
        @test norm(x - (B \ b2)) / norm(x) < 1e-14
    end

    @testset "split QR (overdetermined)" begin
        for (m, n) in [(2, 1), (3, 1), (3, 2), (5, 4), (6, 5), (12, 7), (21, 20), (40, 17), (41, 40)]
            A = randn(ComplexF64, m, n)
            b = randn(ComplexF64, m)
            x = zeros(ComplexF64, n)
            WS = HC.MatrixWorkspace(A)
            HC.update!(WS, A)
            b0 = copy(b)
            ldiv!(x, WS, b)
            @test b == b0
            @test x ≈ A \ b rtol = cond(A) * 1e-13
            # least squares solution: the residual is orthogonal to the range of A
            @test norm(A' * (A * x - b)) <= 1e-10 * norm(A) * norm(b)
            # no allocations
            WS.factorized[] = false
            ldiv!(x, WS, b)
            WS.factorized[] = false
            @test (@allocated ldiv!(x, WS, b)) == 0
            # the factors do not depend on the right-hand side
            x2 = zeros(ComplexF64, n)
            ldiv!(x2, WS, b)
            @test x2 == x
            # a second matrix in the same workspace
            A2 = randn(ComplexF64, m, n)
            HC.update!(WS, A2)
            ldiv!(x, WS, b)
            @test x ≈ A2 \ b rtol = cond(A2) * 1e-13
        end

        # a zero column (zero Householder reflector) and a rank-one perturbation
        A = randn(ComplexF64, 6, 4)
        A[:, 2] .= 0
        A[1, 1] = 0
        WS = HC.MatrixWorkspace(A)
        HC.update!(WS, A)
        HC.factorize!(WS)
        @test iszero(HC.synced_qr(WS).τ[2])
        @test cond(WS, nothing, ones(4)) == Inf

        # the condition numbers read the diagonal of R
        B = randn(ComplexF64, 9, 5)
        WS = HC.MatrixWorkspace(B)
        HC.update!(WS, B)
        d = abs.(diag(HC.synced_qr(WS).factors))
        @test cond(WS, nothing, ones(5)) ≈ maximum(d) / minimum(d)
    end

    @testset "ldiv" begin
        for n in [3, 13, 31] # test that struct array also works
            A = randn(ComplexF64, n, n)
            b = randn(ComplexF64, n)
            x = zeros(ComplexF64, n)
            WS = HC.MatrixWorkspace(A)
            HC.update!(WS, A)
            ldiv!(x, WS, b)
            WS.factorized[] = false
            @allocated ldiv!(x, WS, b)
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
        WB = HC.MatrixWorkspace(B)
        @test 0.1 ≤ opnorm(inv(B), Inf) / HC.inverse_inf_norm_est(WB) ≤ 10
        @test 0.1 ≤
              opnorm(inv(B * D_R), Inf) / HC.inverse_inf_norm_est(WB, nothing, d_r) ≤
              10
        @test 0.1 ≤ opnorm(inv(D_R * B), Inf) / HC.inverse_inf_norm_est(WB, d_r) ≤ 10

        @test 0.1 ≤ cond(B, Inf) / cond(WB) ≤ 10
        @test 0.1 ≤ cond(B * D_R, Inf) / cond(WB, nothing, d_r) ≤ 10
        @test 0.1 ≤ cond(D_R * B, Inf) / cond(WB, d_r, nothing) ≤ 10


        C = inv(D_L) * A
        WC = HC.MatrixWorkspace(C)
        @test 0.1 ≤ opnorm(inv(C), Inf) / HC.inverse_inf_norm_est(WC) ≤ 10
        @test 0.1 ≤ opnorm(inv(D_L * C), Inf) / HC.inverse_inf_norm_est(WC, d_l) ≤ 10
        @test 0.1 ≤
              opnorm(inv(C * D_L), Inf) / HC.inverse_inf_norm_est(WC, nothing, d_l) ≤
              10
        @test 0.1 ≤ cond(C, Inf) / cond(WC) ≤ 10
        @test 0.1 ≤ cond(C * D_L, Inf) / cond(WC, nothing, d_l) ≤ 10
        @test 0.1 ≤ cond(D_L * C, Inf) / cond(WC, d_l, nothing) ≤ 10


        D = inv(D_L) * A * inv(D_R)
        WD = HC.MatrixWorkspace(D)
        @test 0.1 ≤ opnorm(inv(D), Inf) / HC.inverse_inf_norm_est(WD) ≤ 10
        @test 0.1 ≤
              opnorm(inv(D_L * D * D_R), Inf) / HC.inverse_inf_norm_est(WD, d_l, d_r) ≤
              10
        @test HC.inverse_inf_norm_est(WD, d_r, d_l) >
              100 * HC.inverse_inf_norm_est(WD, d_l, d_r)
        @test 0.1 ≤ cond(D, Inf) / cond(WD) ≤ 10
        @test 0.1 ≤ cond(D_L * D * D_R, Inf) / cond(WD, d_l, d_r) ≤ 10
    end

    @testset "Jacobian" begin
        A = randn(ComplexF64, 6, 6)
        b = randn(ComplexF64, 6)
        x̂ = similar(b)
        x = similar(b)
        JM = HC.Jacobian(zeros(ComplexF64, 6, 6))
        @test JM isa HC.Jacobian
        HC.matrix(JM) .= A
        HC.updated!(JM)
        ldiv!(x̂, JM, b)
        ldiv!(x, HC.workspace(JM), b)
        @test x == x̂
    end
end
