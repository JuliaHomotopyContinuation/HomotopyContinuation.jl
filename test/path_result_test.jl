@testset "PathResult" begin
    R = solve(equations(katsura5()))

    @test !isempty(string(R[1]))
end
