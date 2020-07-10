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

        F = System([x^2 + y^2 - 5z^2], [x, y, z])

        W = witness_set(F; compile = false)

        @test dim(W) == 1
        @test codim(W) == 2
        @test degree(W) == 2
        @test solutions(W) isa Vector{Vector{ComplexF64}}
        @test results(W) isa Vector{PathResult}
        @test trace_test(W) < 1e-8

        L = LinearSubspace([1 1 1])
        W_L = witness_set(W, L; compile = false)
        @test degree(W_L) == 2

        L = rand_subspace(3; codim = 1, affine = false)
        W_L = witness_set(W, L; compile = false)
        @test degree(W_L) == 2

        L = rand_subspace(3; codim = 1, affine = true)
        @test_throws ErrorException witness_set(W, L; compile = false)
    end

    @testset "dim / codim" begin
        @var x[1:6]
        homogeneous = true

        f = System([
            rand_poly(x, 2; homogeneous = homogeneous),
            rand_poly(x, 2; homogeneous = homogeneous),
            rand_poly(x, 2; homogeneous = homogeneous),
            rand_poly(x, 2; homogeneous = homogeneous),
        ])

        # quadratic_surface = rand_poly([x,y,z,w], 2; homogeneous = true)
        # sextic_curve = System([cubic_surface, quadratic_surface])
        @test degree(witness_set(f; dim = 1, compile = false)) == 16
        @test degree(witness_set(f; codim = 4, compile = false)) == 16
        @test degree(witness_set(f; compile = false)) == 16

        homogeneous = false
        f = System([
            rand_poly(x, 2; homogeneous = homogeneous),
            rand_poly(x, 2; homogeneous = homogeneous),
            rand_poly(x, 2; homogeneous = homogeneous),
            rand_poly(x, 2; homogeneous = homogeneous),
        ])
        # quadratic_surface = rand_poly([x,y,z,w], 2; homogeneous = true)
        # sextic_curve = System([cubic_surface, quadratic_surface])
        @test degree(witness_set(f; dim = 2, compile = false)) == 16
        @test degree(witness_set(f; codim = 4, compile = false)) == 16
        @test degree(witness_set(f; compile = false)) == 16
    end
end
