@testset "Multiplicities" begin
    V = [randn(3) + im.* randn(3) for i in 1:10]
    W = [V; map(v -> v + 1e-6 .* randn(3), V)]

    M = Solving.multiplicities(V, 1e-5)
    @test length(M) == 10
    @test unique([length(m) for m in M]) == [1]
    @test M[1] == [1]

    N = Solving.multiplicities(W, 1e-5)
    @test length(N) == 10
    @test unique([length(m) for m in N]) == [2]
    @test N[1] == [1, 11]
end
