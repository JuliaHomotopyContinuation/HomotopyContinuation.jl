@testset "Integration Tests" begin
    @test all(issuccess, solve(equations(katsura(10)), show_progress=false))

    @test nfinite(solve(equations(cyclic(6)), affine_tracking=true, show_progress=false)) == 156
    @test nfinite(solve(equations(cyclic(6)), affine_tracking=false, show_progress=false)) == 156
    R = solve(equations(cyclic(7)), seed=14113, show_progress=false)
    @test nfinite(R) == 924
    @test ntracked(R) == 5040
end
