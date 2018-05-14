@testset "integration tests" begin
    @test all(r -> r.returncode == :success, solve(equations(katsura10())).PathResults)
    @test all(r -> r.returncode == :success, solve(equations(katsura10()), report_progress=false).PathResults)

    @test count(r -> r.returncode == :success, solve(equations(cyclic6())).PathResults) == 156
    @test count(r -> r.returncode == :success, solve(equations(cyclic7())).PathResults) == 924
end
