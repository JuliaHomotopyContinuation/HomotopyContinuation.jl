@testset "AffinePatches.OrthogonalPatch" begin
    x = ProjectiveVectors.PVector(rand(Complex{Float64}, 3), 1)
    patch = AffinePatches.state(AffinePatches.OrthogonalPatch(), x)
    @test patch isa AffinePatches.OrthogonalPatchState
    @test AffinePatches.nequations(patch) == 1

    AffinePatches.precondition!(patch, x)
    @test norm(patch.v_conj) ≈ 1.0 atol=1e-15
    u = [0.0im]
    AffinePatches.evaluate!(u, patch, x)
    @test u[1] ≈ 0.0  atol=1e-15

    U = zeros(Complex{Float64}, 1, 3)
    AffinePatches.jacobian!(U, patch, x)
    @test U == reshape(conj.(normalize(x)), 1, 3)
end

@testset "AffinePatches.EmbeddingPatch" begin
    x = ProjectiveVectors.PVector(rand(Complex{Float64}, 3), 1)
    patch = AffinePatches.state(AffinePatches.EmbeddingPatch(), x)
    @test patch isa AffinePatches.EmbeddingPatchState
    @test AffinePatches.nequations(patch) == 1

    AffinePatches.precondition!(patch, x)
    @test x[1] ≈ 1.0 atol=1e-15
    u = [0.0im]
    AffinePatches.evaluate!(u, patch, x)
    @test u[1] ≈ 0.0 atol=1e-15

    U = zeros(Complex{Float64}, 1, 3)
    AffinePatches.jacobian!(U, patch, x)
    @test U == reshape([1.0, 0, 0], 1, 3)
end


@testset "AffinePatches.RandomPatch" begin
    x = ProjectiveVectors.PVector(rand(Complex{Float64}, 3), 1)
    patch = AffinePatches.state(AffinePatches.RandomPatch(), x)
    @test patch isa AffinePatches.RandomPatchState
    @test AffinePatches.nequations(patch) == 1

    AffinePatches.precondition!(patch, x)
    u = [0.0im]
    AffinePatches.evaluate!(u, patch, x)
    @test u[1] ≈ 0.0 atol=1e-15

    U = zeros(Complex{Float64}, 1, 3)
    AffinePatches.jacobian!(U, patch, x)
    @test U == reshape(conj.(patch.v), 1, 3)
end

@testset "AffinePatches.FixedPatch" begin
    x = ProjectiveVectors.PVector(rand(Complex{Float64}, 3), 1)
    patch = AffinePatches.state(AffinePatches.FixedPatch(), x)
    @test patch isa AffinePatches.FixedPatchState
    @test AffinePatches.nequations(patch) == 1

    AffinePatches.precondition!(patch, x)
    u = [0.0im]
    AffinePatches.evaluate!(u, patch, x)
    u
    @test u[1] ≈ 0.0 atol=1e-15

    U = zeros(Complex{Float64}, 1, 3)
    AffinePatches.jacobian!(U, patch, x)
    @test U == reshape(patch.v_conj, 1, 3)
end
