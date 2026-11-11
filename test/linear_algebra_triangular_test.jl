@testset "Triangular solve kernels" begin
    rng = MersenneTwister(0x6271)
    for T in (Float64, ComplexF64), n in (0, 1, 2, 3, 6, 7, 12, 13, 24, 25, 26, 31, 48)
        A = randn(rng, T, n, n) / max(n, 1)
        for i = 1:n
            A[i, i] += 2
        end
        rhs = randn(rng, T, n)
        for (solve!, triangular) in
            ((HC.ldiv_unit_lower!, UnitLowerTriangular), (HC.ldiv_upper!, UpperTriangular))
            expected = triangular(A) \ rhs
            for separate_output in (false, true), strided in (false, true)
                storage = zeros(T, 2n)
                b = strided ? view(storage, 1:2:2n) : zeros(T, n)
                copyto!(b, rhs)
                x = separate_output ? similar(b) : b
                solve!(A, b, x)
                @test x ≈ expected rtol = 5e-14
                @test triangular(A) * x ≈ rhs rtol = 5e-14
            end
        end
    end
end
