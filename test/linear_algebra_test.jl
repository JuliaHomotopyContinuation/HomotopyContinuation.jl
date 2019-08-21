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


		HC.factorization!(WS, LU_FACT)
		@test WS.fact[] == LU_FACT
		WS.factorized[] = true
		HC.factorization!(WS, QR_FACT)
		@test WS.fact[] == QR_FACT
		@test WS.factorized[] == false

		HC.factorization!(WS, LU_FACT)
		@test (@allocated HC.factorize!(WS)) == 0
		HC.factorization!(WS, QR_FACT)
		@test WS.factorized[] == false
		HC.factorize!(WS) # first factorize resizes the work buffers
		WS.factorized[] = false
		@test (@allocated HC.factorize!(WS)) == 0
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
		@test (lu(A) \ b) ≈ x atol=1e-12

		HC.factorization!(WS, QR_FACT)
		x .= 0
		ldiv!(x, WS, b)
		@test (qr(A, Val(true)) \ b) ≈ x atol=1e-12
		WS.factorized[] = false
		@test (@allocated ldiv!(x, WS, b)) == 0

		# overdetermined
		A = randn(ComplexF64, 13, 8)
		x = randn(ComplexF64, 8)
		b = A*x
		WS = HC.MatrixWorkspace(A)
		HC.update!(WS, A)
		ldiv!(x, WS, b)
		WS.factorized[] = false
		x .= 0
		@test (@allocated ldiv!(x, WS, b)) == 0
		@test (qr(A, Val(true)) \ b) ≈ x atol=1e-12
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
		A2 = U * diagm(0=>S) * VT'
		update!(WS, A2)
		cond(A2)
		ldiv!(x̂, WS, b)
		δx1 = HC.iterative_refinement_step!(x, WS, x̂, b, ComplexDF64)
		@test δx1 / norm(x, Inf) > 1e-14
		δx2 = HC.iterative_refinement_step!(WS, x, b,  ComplexDF64)
		@test δx2/ norm(x, Inf) < 1e-16

		δx = HC.iterative_refinement_step!(x, WS, x̂, b)
		@test δx1 / norm(x, Inf) > 1e-14
	end

	@testset "rcond" begin
		A = randn(ComplexF64, 13, 13)
		WS = HC.MatrixWorkspace(A)
		U, S, VT = svd(A)
		S[end] *= 1e-5
		A2 = U * diagm(0=>S) * VT'
		HC.update!(WS, A2)
		HC.factorization!(WS, LU_FACT)
		HC.factorize!(WS)
		@test cond(A2, Inf) ≈ inv(HC.rcond(WS))
		# @test (@allocated rcond(WS)) == 0
		HC.factorization!(WS, QR_FACT)
		@test 0.5 ≤ cond(A2, Inf) / inv(HC.rcond(WS)) ≤ 1.5
		# @test (@allocated rcond(WS)) == 0

		A = randn(ComplexF64, 13, 10)
		WS = HC.MatrixWorkspace(A)
		U, S, VT = svd(A)
		S[end] *= 1e-8
		A2 = U * diagm(0=>S) * VT'
		HC.update!(WS, A2)
		HC.factorize!(WS)
		@test 1 / HC.rcond(WS) ≈ cond(UpperTriangular(WS.qr.R), Inf) rtol=1e-14
		# @test (@allocated rcond(WS)) == 0
	end

	@testset "inf norm est" begin
		A = randn(ComplexF64, 5, 5)
		U, s, Vt = svd(A)
		s[end] *= 1e-6
		B = U*diagm(0=>s)*Vt'

		d = rand(5)*0.001
		D = diagm(0=>d)
		g = rand(5) .* 1000
		G = diagm(0=>g)

		@test 0.1 ≤ opnorm(inv(B), Inf) / HC.inf_norm_est(lu(B)) ≤ 10
		@test 0.1 ≤ opnorm(D*inv(B), Inf) / HC.inf_norm_est(lu(B), nothing, d) ≤ 10
		@test 0.1 ≤ opnorm(inv(B)*G, Inf) / HC.inf_norm_est(lu(B), g, nothing) ≤ 10
		@test 0.1 ≤ opnorm(D*inv(B)*G, Inf) / HC.inf_norm_est(lu(B), g, d) ≤ 10

		WS = HC.MatrixWorkspace(B)
		@test 0.1 ≤ opnorm(inv(B), Inf) / HC.inf_norm_est(WS) ≤ 10
		@test 0.1 ≤ opnorm(D*inv(B), Inf) / HC.inf_norm_est(WS, nothing, d) ≤ 10
		@test 0.1 ≤ opnorm(inv(B)*G, Inf) / HC.inf_norm_est(WS, g, nothing) ≤ 10
		@test 0.1 ≤ opnorm(D*inv(B)*G, Inf) / HC.inf_norm_est(WS, g, d) ≤ 10
	end

	@testset "Hermite Normal Form" begin
        A = [0 3 1 0; -2 2 -1 2; 1 -1 2 3; -3 3 3 2]
        H, U = HC.hnf(A)
        @test A * U == H

        D = zeros(Int, 1,1)
        D[1,1] = -2
        H, U = HC.hnf(D)
        @test H[1,1] == 2
        @test U[1,1] == -1

		# overflow A
		A = [2 4 4 2 -4 -3 -4 3 0 3;
		     -4 -2 -5 -5 1 0 -3 0 -1 -2;
			 -3 5 -3 1 4 3 -1 -4 -1 0;
			 1 3 3 4 -2 -3 -2 -5 -4 -3;
			 2 2 0 2 -4 4 3 -4 -2 -4;
			 -3 -3 1 1 -4 4 0 4 4 -4;
			 -3 -1 -1 -1 2 -4 -3 -4 4 -4;
			 0 1 -2 4 5 4 3 1 -5 2;
			 -5 -4 -5 3 1 5 0 -3 -3 -1;
			 2 -5 -3 -1 -1 5 -4 5 -3 1]
		@test_throws OverflowError HC.hnf(A)
		@test_throws OverflowError HC.hnf(A, Int128)
		H, U = HC.hnf(A, BigInt)
		@test A * U == H
    end
end
