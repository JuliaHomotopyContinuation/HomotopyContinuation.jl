import DoubleFloats: Double64, ComplexDF64

@testset "LinearAlgebra" begin
    @testset "Constructor + basics" begin
        A = rand(ComplexF32, 12, 12)
        WS = HC.MatrixWorkspace(A)
        @test WS isa HC.MatrixWorkspace{Float64}
        @test_throws ArgumentError HC.update!(WS, rand(12, 13))

        WS.factorized[] = true
        B = randn(ComplexF64, 12, 12)
        HC.update!(WS, B)
        @test WS.A == B
        @test WS.A !== B
        @test WS.factorized[] == false


        HC.factorization!(WS, HC.LU_FACT)
        @test WS.fact[] == HC.LU_FACT
        WS.factorized[] = true
        HC.factorization!(WS, HC.QR_FACT)
        @test WS.fact[] == HC.QR_FACT
        @test WS.factorized[] == false

        HC.factorization!(WS, HC.LU_FACT)
        @test (@allocated HC.factorize!(WS)) == 0
        HC.factorization!(WS, HC.QR_FACT)
        @test WS.factorized[] == false
        HC.factorize!(WS) # first factorize resizes the work buffers
        WS.factorized[] = false
        @test (@allocated HC.factorize!(WS)) == 0

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
    end

    @testset "ldiv" begin
        A = randn(ComplexF64, 13, 13)
        b = randn(ComplexF64, 13)
        x = zeros(ComplexF64, 13)
        WS = HC.MatrixWorkspace(A)
        HC.update!(WS, A)
        ldiv!(x, WS, b)
        WS.factorized[] = false
        @test (@allocated ldiv!(x, WS, b)) == 0
        @test (lu(A) \ b) ≈ x rtol = 1e-12

        HC.factorization!(WS, HC.QR_FACT)
        x .= 0
        ldiv!(x, WS, b)
        @test (qr(A, Val(true)) \ b) ≈ x rtol = 1e-12
        WS.factorized[] = false
        @test (@allocated ldiv!(x, WS, b)) == 0

		# overdetermined
        A = randn(ComplexF64, 13, 8)
        x = randn(ComplexF64, 8)
        b = A * x
        WS = HC.MatrixWorkspace(A)
        HC.update!(WS, A)
        @test_throws ArgumentError HC.factorization!(WS, HC.LU_FACT)
        ldiv!(x, WS, b)
        WS.factorized[] = false
        x .= 0
        @test (@allocated ldiv!(x, WS, b)) == 0
        @test (qr(A, Val(true)) \ b) ≈ x rtol = 1e-12
    end

    @testset "(Mixed precision) residual" begin
        A = randn(ComplexF64, 13, 13)
        x = randn(ComplexF64, 13)
        b = A * x
        û = convert.(ComplexF64, ComplexDF64.(A) * ComplexDF64.(x) - ComplexDF64.(b))

        u = zeros(ComplexF64, 13)
		# fixed precision residual should be of order unit roundoff
        HC.residual!(u, A, x, b)
        @test maximum(abs, ComplexDF64.(u) - ComplexDF64.(û)) < 1e-14
        @test maximum(abs, ComplexDF64.(u) - ComplexDF64.(û)) > 1e-18

        HC.residual!(u, A, x, b, Complex{Double64})
        @test maximum(abs, ComplexDF64.(u) - ComplexDF64.(û)) < 1e-26
    end


    @testset "(Mixed precision) iterative refinement" begin
        A = randn(ComplexF64, 4, 4)
        b = randn(ComplexF64, 4)
        x̂ = zeros(ComplexF64, 4)
        x = copy(x̂)
        WS = HC.MatrixWorkspace(A)
        U, S, VT = svd(A)
        S[end] *= 1e-4
        A2 = U * diagm(0 => S) * VT'
        HC.update!(WS, A2)
        ldiv!(x̂, WS, b)
        δx1 = HC.iterative_refinement_step!(x, WS, x̂, b, HC.InfNorm(), ComplexDF64)
        @test δx1 / norm(x, Inf) > 1e-14
        δx2 = HC.iterative_refinement_step!(WS, x, b, HC.InfNorm(), ComplexDF64)
        @test δx2 / norm(x, Inf) < eps()

        δx = HC.iterative_refinement_step!(x, WS, x̂, b)
        @test δx1 / norm(x, Inf) > 1e-14

		# overdetermined
        A = randn(ComplexF64, 4, 3)
        x̂ = zeros(ComplexF64, 3)
        x = copy(x̂)
        WS = HC.MatrixWorkspace(A)
        U, S, VT = svd(A)
        S[end] *= 1e-4
        A2 = U * diagm(0 => S) * VT'
        HC.update!(WS, A2)
        b = A2 * randn(ComplexF64, 3)
        ldiv!(x̂, WS, b)
        δx1 = HC.iterative_refinement_step!(x, WS, x̂, b, HC.InfNorm(), ComplexDF64)
        @test δx1 / norm(x, Inf) > 1e-14
        δx2 = HC.iterative_refinement_step!(WS, x, b, HC.InfNorm(), ComplexDF64)
        @test δx2 / norm(x, Inf) < eps()

        δx = HC.iterative_refinement_step!(x, WS, x̂, b)
        @test δx1 / norm(x, Inf) > 1e-14
    end

    @testset "rcond" begin
        A = randn(ComplexF64, 13, 13)
        WS = HC.MatrixWorkspace(A)
        U, S, VT = svd(A)
        S[end] *= 1e-5
        A2 = U * diagm(0 => S) * VT'
        HC.update!(WS, A2)
        HC.factorization!(WS, HC.LU_FACT)
        HC.factorize!(WS)
        @test cond(A2, Inf) ≈ inv(HC.rcond(WS))
        HC.factorization!(WS, HC.QR_FACT)
        @test 0.5 ≤ cond(A2, Inf) / inv(HC.rcond(WS)) ≤ 1.5
		# @test (@allocated rcond(WS)) == 0

        A = randn(ComplexF64, 13, 10)
        WS = HC.MatrixWorkspace(A)
        U, S, VT = svd(A)
        S[end] *= 1e-8
        A2 = U * diagm(0 => S) * VT'
        HC.update!(WS, A2)
        HC.factorize!(WS)
        @test 1 / HC.rcond(WS) ≈ cond(UpperTriangular(WS.qr.R), Inf) rtol = 1e-14
		# @test (@allocated rcond(WS)) == 0

		# meddle aroudn with the internals to check dead branch
        WS.fact[] = HC.LU_FACT
        @test isnan(HC.rcond(WS))
    end

    @testset "inf norm est" begin
        A = randn(ComplexF64, 5, 5)
        U, s, Vt = svd(A)
        s[end] *= 1e-6
        B = U * diagm(0 => s) * Vt'

        d = rand(5) * 1000
        D = diagm(0 => d)
        g = rand(5) .* 1000
        G = diagm(0 => g)

        @test 0.1 ≤ opnorm(inv(B), Inf) / HC.inf_norm_est(lu(B)) ≤ 10
        @test 0.1 ≤ opnorm(inv(D) * inv(B), Inf) / HC.inf_norm_est(lu(B), nothing, d) ≤ 10
        @test 0.1 ≤ opnorm(inv(B) * G, Inf) / HC.inf_norm_est(lu(B), g, nothing) ≤ 10
        @test 0.1 ≤ opnorm(inv(D) * inv(B) * G, Inf) / HC.inf_norm_est(lu(B), g, d) ≤ 10

        WS = HC.MatrixWorkspace(B)
        @test 0.1 ≤ opnorm(inv(B), Inf) / HC.inf_norm_est(WS) ≤ 10
        @test 0.1 ≤ opnorm(inv(D) * inv(B), Inf) / HC.inf_norm_est(WS, nothing, d) ≤ 10
        @test 0.1 ≤ opnorm(inv(B) * G, Inf) / HC.inf_norm_est(WS, g, nothing) ≤ 10
        @test 0.1 ≤ opnorm(inv(D) * inv(B) * G, Inf) / HC.inf_norm_est(WS, g, d) ≤ 10
    end

    @testset "JacobianMonitor" begin
        A = randn(ComplexF64, 6, 6)
        b = randn(ComplexF64, 6)
        x̂ = similar(b)
        x = similar(b)
        JM = HC.JacobianMonitor(zeros(ComplexF64, 6, 6))
        @test JM isa HC.JacobianMonitor
        HC.jacobian(JM) .= A
        HC.updated!(JM)
        ldiv!(x̂, JM, b)
        ldiv!(x, jacobian(JM), b)
        @test x == x̂
    end
    @testset "Forward error" begin
        A = [
            1.81451 + 0.46473im -0.473651 + 0.813117im 0.130822 - 0.254704im -0.203277 +
                                                                             0.397106im;
            -2.04468 - 1.59228im 0.441157 - 1.39994im -0.0984804 + 1.42619im -0.818292 -
                                                                             0.487898im;
            0.35501 + 0.0235574im 1.32618 + 0.577264im 1.12561 + 0.43383im 0.579533 +
                                                                           1.49484im;
            0.144381 + 0.104175im 0.665384 - 1.43389im 0.22097 + 0.258016im 0.374103 -
                                                                            0.635066im
        ]
        b = [
            -0.14581 - 0.868847im,
            0.0162454 + 0.91156im,
            -0.133769 + 0.510834im,
            -0.867022 + 0.752196im
        ]
        x̂ = similar(b)
        x = similar(b)
        JM = HC.JacobianMonitor(A)
        HC.updated!(JM)
        ldiv!(x̂, JM, b)

        ferr = HC.forward_err!(x, JM, x̂, b, HC.InfNorm())
        @test HC.forward_err(JM) == ferr
        @test ferr < eps() * cond(A)
        @test ferr > 1e-12
        ferrD64 = HC.forward_err!(x, JM, x̂, b, HC.InfNorm(), ComplexDF64)
        @test 0.5 ≤ ferrD64 / ferr ≤ 2

        ferr = HC.forward_err!(JM, x̂, b, HC.InfNorm(), ComplexF64)
        @test ferr < eps() * cond(A)
        @test ferr > 1e-12

        D = diagm(0 => [10.0^(2i) for i = 1:4])
        JM = HC.JacobianMonitor(D * A)
        HC.updated!(JM)
        ldiv!(x̂, JM, D * b)
        conderr = HC.cond(D * A) * eps()
        sferr = HC.strong_forward_err!(x, JM, x̂, D * b, HC.InfNorm())
        @test HC.forward_err(JM) == sferr
        ferr = HC.forward_err!(x, JM, x̂, D * b, HC.InfNorm())
        @test ferr < sferr < conderr * 1e-4 # check that error estimate is actually substantially less

        sferr2 = HC.strong_forward_err!(JM, x̂, D * b, HC.InfNorm())
        @test sferr == sferr2
        @test x̂ ≈ x

        DR = diagm(0 => [10.0^(2i - 4) for i = 1:4])
        JM = HC.JacobianMonitor(A * DR)
        ldiv!(x̂, JM, b)

        WN = HC.WeightedNorm(HC.InfNorm(), x̂)
        WN .= abs.(x̂)

        wferr = HC.forward_err!(x, JM, x̂, b, WN)
        ferr = HC.forward_err!(x, JM, x̂, b, HC.InfNorm())

        strong_wferr = HC.strong_forward_err!(x, JM, x̂, b, WN)
        strong_ferr = HC.strong_forward_err!(x, JM, x̂, b, HC.InfNorm())
        @test ferr ≤ strong_ferr ≤ eps() * inv(HC.rcond(jacobian(JM)))
        @test wferr ≤ strong_wferr ≤ eps() * inv(HC.rcond(jacobian(JM)))

		# overdetermined
        A = randn(ComplexF64, 3, 2)
        x = randn(ComplexF64, 2)
        b = A * x
        x̂ = similar(x)

        JM = HC.JacobianMonitor(A)
        ldiv!(x̂, JM, b)

        ferr = HC.forward_err!(x, JM, x̂, b, HC.InfNorm())
        strong_ferr = HC.strong_forward_err!(x, JM, x̂, b, HC.InfNorm())
        @test ferr < 1e-15
        @test ferr == strong_ferr # test fallback to ferr for qr
    end

    @testset "cond" begin
        A = randn(ComplexF64, 13, 13)
        D = diagm(0 => [2.0^(i - 5) for i = 1:13])
        JM = HC.JacobianMonitor(D * A)
        comp_cond = HC.cond!(JM)
        @test comp_cond === cond(JM)
        B = diagm(0 => inv.(JM.J.row_scaling)) * D * A
        @test abs(cond(JM) - cond(B, Inf)) / cond(JM) < 2

        A = randn(ComplexF64, 13, 12)
        JM = HC.JacobianMonitor(D * A)
        comp_cond = HC.cond!(JM)
        @test comp_cond === cond(JM)
    end

    @testset "JacobianMonitor ldiv updates" begin
        A = randn(ComplexF64, 6, 6)
        b = randn(ComplexF64, 6)
        x̂ = similar(b)
        x = similar(b)
        JM = HC.JacobianMonitor(zeros(ComplexF64, 6, 6))
        @test JM isa HC.JacobianMonitor
        ldiv!(x, JM, b, HC.InfNorm(), HC.JAC_MONITOR_UPDATE_FERR)
        @test HC.forward_err(JM) != 0.0
        ldiv!(x, JM, b, HC.InfNorm(), HC.JAC_MONITOR_UPDATE_COND)
        @test cond(JM) != 1.0
        JM = HC.JacobianMonitor(zeros(ComplexF64, 6, 6))
        ldiv!(x, JM, b, HC.InfNorm(), HC.JAC_MONITOR_UPDATE_ALL)
        @test cond(JM) != 1.0
        @test HC.forward_err(JM) != 0.0
    end

    @testset "Hermite Normal Form" begin
        A = [0 3 1 0; -2 2 -1 2; 1 -1 2 3; -3 3 3 2]
        H, U = HC.hnf(A)
        @test A * U == H

        D = zeros(Int, 1, 1)
        D[1, 1] = -2
        H, U = HC.hnf(D)
        @test H[1, 1] == 2
        @test U[1, 1] == -1

		# overflow A
        A = [
            2 4 4 2 -4 -3 -4 3 0 3;
            -4 -2 -5 -5 1 0 -3 0 -1 -2;
            -3 5 -3 1 4 3 -1 -4 -1 0;
            1 3 3 4 -2 -3 -2 -5 -4 -3;
            2 2 0 2 -4 4 3 -4 -2 -4;
            -3 -3 1 1 -4 4 0 4 4 -4;
            -3 -1 -1 -1 2 -4 -3 -4 4 -4;
            0 1 -2 4 5 4 3 1 -5 2;
            -5 -4 -5 3 1 5 0 -3 -3 -1;
            2 -5 -3 -1 -1 5 -4 5 -3 1
        ]
        @test_throws OverflowError HC.hnf(A)
        @test_throws OverflowError HC.hnf(A, Int128)
        H, U = HC.hnf(A, BigInt)
        @test A * U == H
    end
end
