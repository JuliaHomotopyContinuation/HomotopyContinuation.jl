@testset "LinearAlgebra" begin
    @testset "Constructor + basics" begin
        A = rand(ComplexF32, 12, 12)
        WS = HC2.MatrixWorkspace(A)
        @test WS isa HC2.MatrixWorkspace{Float64}

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

    @testset "Jacobian" begin
        A = randn(ComplexF64, 6, 6)
        b = randn(ComplexF64, 6)
        x̂ = similar(b)
        x = similar(b)
        JM = HC2.Jacobian(zeros(ComplexF64, 6, 6))
        @test JM isa HC2.Jacobian
        HC2.jacobian(JM) .= A
        HC2.updated!(JM)
        ldiv!(x̂, JM, b)
        ldiv!(x, HC2.jacobian(JM), b)
        @test x == x̂
    end
end
