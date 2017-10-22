@testset "solve" begin
    PolyImpl.@polyvar x
    res = solve([(x - 2.0) * (x - (2.5+ 4.0im))])

    @test res[1].retcode == :success
    @test res[2].retcode == :success

    res = solve((x - 2.0) * (x - (2.5+ 4.0im)), homotopytype=StraightLineHomotopy)
    @test res[1].retcode == :success
    @test res[2].retcode == :success
end
