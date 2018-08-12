using HomotopyContinuation.ProjectiveVectors

@testset "ProjectiveVectors" begin
    x = rand(Complex{Float64}, 4)

    z = PVector(x, 1)
    infinity_norm(z)
    @test z isa PVector{Complex{Float64}}
    @test z.homvar == 1
    @test PVector(x, 3).homvar == 3
    @test_throws AssertionError PVector(x, 5)

    @test embed([2, 2, 2], 1) == PVector([1, 2, 2, 2], 1)
    @test embed([2, 2, 2], 3) == PVector([2, 2, 1, 2], 3)

    @test length(z) == 4
    @test z[3] == x[3]
    @test z[2:4] == x[2:4]
    @test z[2:end] == x[2:end]

    @test all(affine(z) .≈ (z[2:end] ./ z[1]))
    LinearAlgebra.normalize!(z)

    inf_z = infinity_norm(z)
    affine!(z)
    @test unsafe_infinity_norm(z, z) ≈ 0.0
    @test infinity_norm(z, z) ≈ 0.0 atol=1e-15

    @test at_infinity(z, 1e5) == false
    @test at_infinity(PVector([2, 1e-7, 4], 2), 1e5) == true
    @test at_infinity(PVector([2, 1e-7, 4], 2), 5e7) == false
    @test at_infinity(PVector([2.3 + 0.2im, 1e-5, 4], 2), 1e5) == true
    @test at_infinity(PVector([2.3 + 0.2im, 1e-4, 4], 2), 1e5) == false

    z2 = PVector(x)
    @test z2 isa PVector{<:Complex, Nothing}
    @test_throws MethodError affine(z2)
    @test affine(z2, 1) == affine(z)

    # x1 = PVector(rand(Complex{Float64}, 3), 2)
    # x2 = PVector(rand(Complex{Float64}, 7), 3)
    # z2 = ProdPVector([x1, x2])
    # @test affine(z2) == [affine(x1), affine(x2)]
    #
    # @test length.(pvectors(z2)) == [3, 7]
    # @test homvar.(pvectors(z2)) == [2, 3]
    #
    # @test at_infinity(z2, 1e5) == false
    # @test at_infinity(z2, 1e-1) == true
    # normalize!(z2)
    # all(v -> norm(v) ≈ 1.0, pvectors(z2))
    #
    # infnorm = infinity_norm(z2)
    # affine!(z2)
    # @test sqrt(maximum(abs2, raw(z2))) ≈ infnorm
    # normalize!(z2)
    # @test infinity_norm(z2) ≈ infnorm
    #
    # @test z2 == copy(z2)
    #
    # @test abs(infinity_norm(z2, z2)) < 1e-14
    #
    #
    # converted = similar(ProdPVector([PVector(rand(Float64, 3), 2), PVector(rand(Float64, 7), 3)]), Complex{Float64})
    # @test converted isa ProdPVector{Complex{Float64}}
end
