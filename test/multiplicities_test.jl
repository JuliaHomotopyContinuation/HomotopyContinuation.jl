@testset "Multiplicities" begin
    V = [randn(3) + im.* randn(3) for i in 1:10]
    W = [map(v -> v + [1e-7 * v[1] ;1.0;1.0], V); V; map(v -> v + [1e-6 * v[1]; 0; 0], V)]

    M = Solving.multiplicities(V, 1e-5)
    @test length(M) == 10
    @test unique([length(m) for m in M]) == [1]
    @test M[1] == [1]

    N = Solving.multiplicities(W, 1e-5)
    @test length(N) == 20
    @test unique([length(m) for m in N]) == [1, 2]
    @test N[1] == [1]
    @test N[2] == [11, 21]
end
