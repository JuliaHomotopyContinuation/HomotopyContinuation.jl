@testset "integration tests" begin
    @test all(r -> r.returncode == :success, solve(equations(katsura8())))

    @test count(r -> r.returncode == :success, solve(equations(cyclic6()))) == 156
    @test count(r -> r.returncode == :success, solve(equations(cyclic7()))) == 924
end
