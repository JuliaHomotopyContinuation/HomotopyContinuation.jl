@testset "PathResult" begin
    R = solve(equations(katsura5()))
    @test R isa Solving.AffineResult
    @test natinfinity(R) == 0
    @test nfinite(R) == 32
    @test finite(R) isa Vector{<:Solving.PathResult}
    @test string(collect(R)[1]) isa String
end
