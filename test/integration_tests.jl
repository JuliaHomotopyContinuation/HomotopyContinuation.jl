@testset "integration tests" begin
    @test all(issuccess, solve(equations(katsura(10))))
    @test all(issuccess, solve(equations(katsura(10)), report_progress=false))

    @test nfinite(solve(equations(cyclic(6)))) == 156
    R = solve(equations(cyclic(7)))
    @test nfinite(R) == 924
    @test nfailed(R) == 0
end
