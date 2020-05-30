@testset "Norms" begin
    x = [2.0im, 3.0 - 1im, 5.0 + 2.0im]
    y = [-2.0im, 3.0 - 1im, 5.0 + 2.0im]
    Norm = HC.InfNorm()
    @test Norm isa HC.AbstractNorm
    @test Norm(x) ≈ abs(5.0 + 2.0im)
    @test Norm(x, x) ≈ 0.0
    @test Norm(x, y) ≈ abs(4im)

    weighted_norm = HC.WeightedNorm(Norm, x)
    @test HC.weights(weighted_norm) == ones(3)
    weighted_norm .= 2.0
    @test HC.weights(weighted_norm) == [2.0, 2, 2]
    weighted_norm[1] = 4.0
    @test weighted_norm isa HC.AbstractNorm
    @test weighted_norm(x) ≈ 0.5 * abs(5 + 2.0im)
    @test weighted_norm(x, x) ≈ 0.0
    @test weighted_norm(x, y) ≈ 0.25 * abs(4im)

    Norm = HC.EuclideanNorm()
    @test Norm isa HC.AbstractNorm

    @test Norm(x) ≈ LinearAlgebra.norm(x, 2)
    @test Norm(x, x) ≈ 0.0
    @test Norm(x, y) ≈ abs(4im)

    weighted_norm = HC.WeightedNorm([4.0, 2.0, 2.0], Norm)
    @test weighted_norm isa HC.AbstractNorm
    @test weighted_norm(x) ≈ sqrt(abs2(x[1]) / 16 + abs2(x[2]) / 4 + abs2(x[3]) / 4)
    @test weighted_norm(x, x) ≈ 0.0
    @test weighted_norm(x, y) ≈ sqrt(abs2(4im) / 16)
end
