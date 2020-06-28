@testset "Witness Sets" begin

    @testset "affine" begin
        @var x y

        F = System([x^2 + y^2 - 5], [x, y])

        W = witness_set(F)

        @test dim(W) == 1
        @test codim(W) == 1
        @test degree(W) == 2
        @test solutions(W) isa Vector{Vector{ComplexF64}}
        @test results(W) isa Vector{PathResult}

        L = LinearSubspace([1 1], [-1])

        W_L = witness_set(W, L)
        @test degree(W_L) == 2
        @test sort(real.(solutions(W_L))) ≈ [[-2, 1], [1, -2]]
        @test linear_subspace(W_L) == L

        @test trace_test(W) < 1e-8
    end

    @testset "projective" begin
        @var x y z

        F = System([x^2 + y^2 - 5z^2], [x, y,z])

        W = witness_set(F)

        @test dim(W) == 1
        @test codim(W) == 1
        @test degree(W) == 2
        @test solutions(W) isa Vector{Vector{ComplexF64}}
        @test results(W) isa Vector{PathResult}

        L = LinearSubspace([1 1 1])

        W_L = witness_set(W, L)
        @test degree(W_L) == 2

        @test trace_test(W) < 1e-8
    end
end
