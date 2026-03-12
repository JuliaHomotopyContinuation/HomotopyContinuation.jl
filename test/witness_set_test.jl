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

        L = rand_subspace([x, y, z]; codim = 1, affine = false)
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

    @var x, y, z
    p = (x * y - x^2) + 1 - z
    q = x^4 + x^2 - y - 1
    F = [
        p * q * (x - 3) * (x - 5)
        p * q * (y - 3) * (y - 5)
        p * (z - 3) * (z - 5)
    ]

    @testset "membership" begin
        W = witness_set(F; codim = 2)

        p = randn(3)
        q = solutions(W)[1]

        @test !membership(p, W)
        @test membership(q, W; show_progress = false)
        a = membership([p, q], W; show_progress = false)
        @test a == [false, true]
    end

    @testset "intersect" begin
        H = [witness_set(f) for f in F]
        B = intersect(H[1], H[2])
        C = vcat([intersect(Hi, H[3]; show_progress = false) for Hi in B]...)
        @test degree.(C) == [2, 8, 8]
    end

end
