@testset "solve" begin
    PolyImpl.@polyvar x
    res = solve([(x - 2.0) * (x - (2.5+ 4.0im))], endgame_start=0.0)

    @test issuccessfull(res[1])
    @test issuccessfull(res[2])

    res = solve((x - 2.0) * (x - (2.5+ 4.0im)), homotopytype=StraightLineHomotopy, endgame_start=0.0)
    @test issuccessfull(res[1])
    @test issuccessfull(res[2])
end
