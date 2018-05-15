@testset "integration tests" begin
    @test all(issuccess, solve(equations(katsura10())))
    @test all(issuccess, solve(equations(katsura10()), report_progress=false))

    @test nfinite(solve(equations(cyclic6()))) == 156
    @test nfinite(solve(equations(cyclic7()))) == 924
end
