@testset "AffinePatches" begin
    @testset "OrthogonalPatch" begin
        x = ProjectiveVectors.PVector(rand(Complex{Float64}, 3), 1)
        patch = AffinePatches.state(AffinePatches.OrthogonalPatch(), x)
        @test patch isa AffinePatches.OrthogonalPatchState
        @test AffinePatches.nequations(patch) == 1

        AffinePatches.setup!(patch, x)
        @test norm(patch.v̄) ≈ 1.0 atol=1e-15
        u = [0.0im]
        AffinePatches.evaluate!(u, patch, x)
        @test u[1] ≈ 0.0  atol=1e-15

        U = zeros(Complex{Float64}, 1, 3)
        AffinePatches.jacobian!(U, patch, x)
        @test U ≈ reshape(conj.(LinearAlgebra.normalize(x)), 1, 3)

        y = ProjectiveVectors.PVector(randn(ComplexF64, 3), 1)
        AffinePatches.onpatch!(y, patch)
        AffinePatches.evaluate!(u, patch, y)
        @test u[1] ≈ 0.0  atol=1e-15
    end

    @testset "EmbeddingPatch" begin
        x = ProjectiveVectors.PVector(rand(Complex{Float64}, 3), 1)
        patch = AffinePatches.state(AffinePatches.EmbeddingPatch(), x)
        @test patch isa AffinePatches.EmbeddingPatchState
        @test AffinePatches.nequations(patch) == 1

        AffinePatches.setup!(patch, x)
        @test x[1] ≈ 1.0 atol=1e-15
        u = [0.0im]
        AffinePatches.evaluate!(u, patch, x)
        @test u[1] ≈ 0.0 atol=1e-15

        U = zeros(Complex{Float64}, 1, 3)
        AffinePatches.jacobian!(U, patch, x)
        @test U == reshape([1.0, 0, 0], 1, 3)

        y = ProjectiveVectors.PVector(randn(ComplexF64, 3), 1)
        AffinePatches.onpatch!(y, patch)
        AffinePatches.evaluate!(u, patch, y)
        @test u[1] ≈ 0.0  atol=1e-15
    end


    @testset "RandomPatch" begin
        x = ProjectiveVectors.PVector(rand(Complex{Float64}, 3), 1)
        patch = AffinePatches.state(AffinePatches.RandomPatch(), x)
        @test patch isa AffinePatches.RandomPatchState
        @test AffinePatches.nequations(patch) == 1

        AffinePatches.setup!(patch, x)
        u = [0.0im]
        AffinePatches.evaluate!(u, patch, x)
        @test u[1] ≈ 0.0 atol=1e-15

        U = zeros(Complex{Float64}, 1, 3)
        AffinePatches.jacobian!(U, patch, x)
        @test U == reshape(conj.(patch.v), 1, 3)

        y = ProjectiveVectors.PVector(randn(ComplexF64, 3), 1)
        AffinePatches.onpatch!(y, patch)
        AffinePatches.evaluate!(u, patch, y)
        @test u[1] ≈ 0.0  atol=1e-15
    end

    @testset "FixedPatch" begin
        x = ProjectiveVectors.PVector(rand(Complex{Float64}, 3), 1)
        patch = AffinePatches.state(AffinePatches.FixedPatch(), x)
        @test patch isa AffinePatches.FixedPatchState
        @test AffinePatches.nequations(patch) == 1

        AffinePatches.setup!(patch, x)
        u = [0.0im]
        AffinePatches.evaluate!(u, patch, x)
        u
        @test u[1] ≈ 0.0 atol=1e-15

        U = zeros(Complex{Float64}, 1, 3)
        AffinePatches.jacobian!(U, patch, x)
        @test U ≈ reshape(patch.v̄, 1, 3)

        y = ProjectiveVectors.PVector(randn(ComplexF64, 3), 1)
        AffinePatches.onpatch!(y, patch)
        AffinePatches.evaluate!(u, patch, y)
        @test u[1] ≈ 0.0  atol=1e-15
    end
end
