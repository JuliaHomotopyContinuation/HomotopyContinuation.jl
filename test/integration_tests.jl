@testset "integration tests" begin
    @test all(r -> r.returncode == :success, solve(equations(katsura10())).path_results)
    @test all(r -> r.returncode == :success, solve(equations(katsura10()), report_progress=false).path_results)

    @test count(r -> r.returncode == :success, solve(equations(cyclic6())).path_results) == 156
    @test count(r -> r.returncode == :success, solve(equations(cyclic7())).path_results) == 924
end
