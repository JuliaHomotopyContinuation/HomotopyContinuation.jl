function Base.isapprox(xs::NTuple{N}, ys::NTuple{N}; kwargs...) where {N}
    for i = 1:N
        if !isapprox(xs[i], ys[i]; kwargs...)
            return false
        end
    end
    true
end

@testset "AffinePatches" begin
    @testset "OrthogonalPatch" begin
        x = ProjectiveVectors.embed(rand(ComplexF64, 2))
        patch = state(OrthogonalPatch(), x)
        @test patch isa HomotopyContinuation.OrthogonalPatchState
        @test nequations(patch) == 1

        init!(patch, x)
        @test norm(patch.v̄) ≈ (1.0,) atol = 1e-15
        u = [0.0im]
        evaluate!(u, patch, x)
        @test u[1] ≈ 0.0 atol = 1e-15

        U = zeros(Complex{Float64}, 1, 3)
        jacobian!(U, patch, x)
        @test U ≈ reshape(conj.(LinearAlgebra.normalize(x)), 1, 3)

        y = ProjectiveVectors.embed(rand(ComplexF64, 2))
        onpatch!(y, patch)
        evaluate!(u, patch, y)
        @test u[1] ≈ 0.0 atol = 1e-15
    end

    @testset "EmbeddingPatch" begin
        x = ProjectiveVectors.embed(rand(ComplexF64, 2))
        patch = state(EmbeddingPatch(), x)
        @test patch isa HomotopyContinuation.EmbeddingPatchState
        @test nequations(patch) == 1

        init!(patch, x)
        @test x[3] ≈ 1.0 atol = 1e-15
        u = [0.0im]
        evaluate!(u, patch, x)
        @test u[1] ≈ 0.0 atol = 1e-15

        U = zeros(Complex{Float64}, 1, 3)
        jacobian!(U, patch, x)
        @test U == reshape([0.0, 0, 1], 1, 3)

        y = ProjectiveVectors.embed(rand(ComplexF64, 2))
        onpatch!(y, patch)
        evaluate!(u, patch, y)
        @test u[1] ≈ 0.0 atol = 1e-15
    end


    @testset "RandomPatch" begin
        x = ProjectiveVectors.embed(rand(ComplexF64, 2))
        patch = state(RandomPatch(), x)
        @test patch isa HomotopyContinuation.RandomPatchState
        @test nequations(patch) == 1

        init!(patch, x)
        u = [0.0im]
        evaluate!(u, patch, x)
        @test u[1] ≈ 0.0 atol = 1e-14

        U = zeros(Complex{Float64}, 1, 3)
        jacobian!(U, patch, x)
        @test U == reshape(patch.v, 1, 3)

        y = ProjectiveVectors.embed(rand(ComplexF64, 2))
        onpatch!(y, patch)
        evaluate!(u, patch, y)
        @test u[1] ≈ 0.0 atol = 1e-14
    end

    @testset "FixedPatch" begin
        x = ProjectiveVectors.embed(rand(ComplexF64, 2))
        patch = state(FixedPatch(), x)
        @test patch isa HomotopyContinuation.FixedPatchState
        @test nequations(patch) == 1

        init!(patch, x)
        u = [0.0im]
        evaluate!(u, patch, x)
        u
        @test u[1] ≈ 0.0 atol = 1e-15

        U = zeros(Complex{Float64}, 1, 3)
        jacobian!(U, patch, x)
        @test U ≈ reshape(patch.v̄, 1, 3)

        y = ProjectiveVectors.embed(rand(ComplexF64, 2))
        onpatch!(y, patch)
        evaluate!(u, patch, y)
        @test u[1] ≈ 0.0 atol = 1e-15
    end
end
