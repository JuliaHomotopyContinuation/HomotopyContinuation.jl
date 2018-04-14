@testset "integration tests" begin
    @test all(r -> r.returncode == :success, solve(equations(katsura8())))
end
