@testset "Norms" begin
    x = [2.0im, 3.0 - 1im, 5.0 + 2.0im]
    y = [-2.0im, 3.0 - 1im, 5.0 + 2.0im]
    Norm = HC.InfNorm()
    @test Norm isa HC.AbstractNorm
    @test Norm(x) ≈ abs(5.0 + 2.0im)
    @test Norm(x, x) ≈ 0.0
    @test Norm(x, y) ≈ abs(4im)

    WeightedNorm = HC.WeightedNorm(Norm, x)
    @test HC.weights(WeightedNorm) == ones(3)
    WeightedNorm .= 2.0
    @test HC.weights(WeightedNorm) == [2.0, 2, 2]
    WeightedNorm[1] = 4.0
    @test WeightedNorm isa HC.AbstractNorm
    @test WeightedNorm(x) ≈ 0.5 * abs(5 + 2.0im)
    @test WeightedNorm(x, x) ≈ 0.0
    @test WeightedNorm(x, y) ≈ 0.25 * abs(4im)

    Norm = HC.EuclideanNorm()
    @test Norm isa HC.AbstractNorm

    @test Norm(x) ≈ LinearAlgebra.norm(x, 2)
    @test Norm(x, x) ≈ 0.0
    @test Norm(x, y) ≈ abs(4im)

    WeightedNorm = HC.WeightedNorm([4.0, 2.0, 2.0], Norm)
    @test WeightedNorm isa HC.AbstractNorm
    @test WeightedNorm(x) ≈ sqrt(abs2(x[1]) / 16 + abs2(x[2]) / 4 + abs2(x[3]) / 4)
    @test WeightedNorm(x, x) ≈ 0.0
    @test WeightedNorm(x, y) ≈ sqrt(abs2(4im) / 16)
end
