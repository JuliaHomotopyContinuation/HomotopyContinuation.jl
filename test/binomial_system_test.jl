import MixedSubdivisions
const MS = MixedSubdivisions

@testset "Binomial System Test" begin
    @testset "Hermite Normal Form" begin
        A = [0 3 1 0; -2 2 -1 2; 1 -1 2 3; -3 3 3 2]
        H, U = HC2.hnf(A)
        @test A * U == H

        D = zeros(Int, 1, 1)
        D[1, 1] = -2
        H, U = HC2.hnf(D)
        @test H[1, 1] == 2
        @test U[1, 1] == -1

        # overflow A
        A = [
            2 4 4 2 -4 -3 -4 3 0 3
            -4 -2 -5 -5 1 0 -3 0 -1 -2
            -3 5 -3 1 4 3 -1 -4 -1 0
            1 3 3 4 -2 -3 -2 -5 -4 -3
            2 2 0 2 -4 4 3 -4 -2 -4
            -3 -3 1 1 -4 4 0 4 4 -4
            -3 -1 -1 -1 2 -4 -3 -4 4 -4
            0 1 -2 4 5 4 3 1 -5 2
            -5 -4 -5 3 1 5 0 -3 -3 -1
            2 -5 -3 -1 -1 5 -4 5 -3 1
        ]
        @test_throws OverflowError HC2.hnf(A)
        H, U = HC2.hnf(A, BigInt)
        @test A * U == H
    end

    @testset "Solving (unit b)" begin
        A = [0 1 1 0 -1; 0 0 0 1 -1; -1 0 -1 0 -1; 0 -1 0 -1 -1; 1 0 0 0 -1]
        b = [
            0.9053223983046926 + 0.4247250347316951im,
            0.7000429487508004 - 0.7141007421255663im,
            0.018811552539476483 - 0.9998230470893611im,
            -0.9983533473373446 + 0.057363698105327716im,
            0.9999409355072137 - 0.010868555421874721im,
        ]
        X = HC2.solve!(HC2.BinomialSystemSolver(A, b))
        @test maximum(eachcol(X)) do x
            maximum(abs.([prod(x .^ a) for a in eachcol(A)] - b) ./ abs.(b))
        end < 1e-14

        A = [
            3 2 3 4 2 5
            6 -3 8 -3 8 7
            -2 -5 7 3 6 5
            1 2 3 4 5 6
            1 0 2 0 5 0
            1 2 0 -3 0 5
        ]
        b = randn(ComplexF64, 6)
        b ./= abs.(b)
        X = HC2.solve!(HC2.BinomialSystemSolver(A, b))
        @test maximum(eachcol(X)) do x
            maximum(abs.([prod(x .^ a) for a in eachcol(A)] - b) ./ abs.(b))
        end < 1e-12
    end

    @testset "General solving" begin
        A = [0 1 1 0 -1; 0 0 0 1 -1; -1 0 -1 0 -1; 0 -1 0 -1 -1; 1 0 0 0 -1]
        b = randn(ComplexF64, size(A, 1))
        X = HC2.solve!(HC2.BinomialSystemSolver(A, b))
        @test maximum(eachcol(X)) do x
            maximum(abs.([prod(x .^ a) for a in eachcol(A)] - b) ./ abs.(b))
        end < 1e-14

        A = [
            3 2 3 4 2 5
            6 -3 8 -3 8 7
            -2 -5 7 3 6 5
            1 2 3 4 5 6
            1 0 2 0 5 0
            1 2 0 -3 0 5
        ]
        b = randn(ComplexF64, 6)
        X = HC2.solve!(HC2.BinomialSystemSolver(A, b))
        @test maximum(eachcol(X)) do x
            maximum(abs.([prod(x .^ a) for a in eachcol(A)] - b) ./ abs.(b))
        end < 1e-12
    end

    @testset "MixedCell" begin
        f = cyclic(5)
        supp = first.(ModelKit.exponents_coefficients.(f, Ref(variables(f))))

        cells, lift = MS.fine_mixed_cells(supp)
        coeffs = map(l -> randn(ComplexF64, length(l)), lift)

        cell = first(cells)
        BSS = HC2.BinomialSystemSolver(5)
        X = HC2.solve(BSS, supp, coeffs, cell)
        @test maximum(eachcol(X)) do x
            maximum(abs.([prod(x .^ a) for a in eachcol(BSS.A)] - BSS.b) ./ abs.(BSS.b))
        end < 1e-12
    end
end
