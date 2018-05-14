@testset "PathResult" begin
    R = solve(equations(katsura5())).PathResults

    @test !isempty(string(R[1]))
end
