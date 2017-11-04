@testset "solve" begin
    PolyImpl.@polyvar x
    # default
    res = solve([(x - 2.0) * (x - (2.5+ 4.0im))])

    @test res[1].retcode == :success
    @test res[2].retcode == :success

    # just with a single polynomial
    res = solve((x - 2.0) * (x - (2.5+ 4.0im)))
    @test res[1].retcode == :success
    @test res[2].retcode == :success

    # other homotopy type
    res = solve((x - 2.0) * (x - (2.5+ 4.0im)), homotopytype=GeodesicOnTheSphere)
    @test res[1].retcode == :success
    @test res[2].retcode == :success

    # non float polynomial, non complex
    res = solve((x - 2) * (x - 4), homotopytype=GeodesicOnTheSphere)
    @test res[1].retcode == :success
    @test res[2].retcode == :success

    # non float homotopy
    H = StraightLineHomotopy([(x - 4) * (x + (2 - 4im))],[(x - 2) * (x - (2 + 4im))])
    @test solve(H, [[4], [-2 + 4im]]) isa Result{Complex128}

    # non float homotopy, non complex
    H = StraightLineHomotopy([(x - 4) * (x + 4)], [(x - 2) * (x + 2)])
    @test solve(H, [[4], [-4]]) isa Result{Complex128}
end
