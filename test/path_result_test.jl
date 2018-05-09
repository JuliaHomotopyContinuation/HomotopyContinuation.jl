@testset "PathResult" begin
    R = solve(equations(katsura5())).path_results

    @test !isempty(string(R[1]))
end
