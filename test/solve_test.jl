@testset "solve" begin
    PolyImpl.@polyvar x
    # default
    res = solve([(x - 2.0) * (x - (2.5+ 4.0im))])
    @test res[1].returncode == :success && res[2].returncode == :success

    # just with a single polynomial
    res = solve((x - 2.0) * (x - (2.5+ 4.0im)))
    @test res[1].returncode == :success
    @test res[2].returncode == :success

    # other homotopy type
    res = solve((x - 2.0) * (x - (2.5+ 4.0im)), GeodesicOnTheSphere)
    @test res[1].returncode == :success
    @test res[2].returncode == :success

    res = solve((x - 2.0) * (x - (2.5+ 4.0im)), GeodesicOnTheSphere{Complex128}, SphericalPredictorCorrector())
    @test res[1].returncode == :success && res[2].returncode == :success

    res = solve((x - 2.0) * (x - (2.5+ 4.0im)), GeodesicOnTheSphere{Complex128}, SphericalPredictorCorrector(), CauchyEndgame())
    @test res[1].returncode == :success && res[2].returncode == :success

    res = solve((x - 2.0) * (x - (2.5+ 4.0im)), SphericalPredictorCorrector())
    @test res[1].returncode == :success && res[2].returncode == :success
    res = solve((x - 2.0) * (x - (2.5+ 4.0im)), SphericalPredictorCorrector(), CauchyEndgame())
    @test res[1].returncode == :success && res[2].returncode == :success

    # non float polynomial, non complex
    res = solve((x - 2) * (x - 4), GeodesicOnTheSphere)
    @test res[1].returncode == :success
    @test res[2].returncode == :success

    # non float homotopy
    H = StraightLineHomotopy([(x - 4) * (x + (2 - 4im))],[(x - 2) * (x - (2 + 4im))])
    @test solve(H, [[4], [-2 + 4im]]) isa Result{Complex128}

    # non float homotopy, non complex
    H = StraightLineHomotopy([(x - 4) * (x + 4)], [(x - 2) * (x + 2)])
    @test solve(H, [[4], [-4]]) isa Result{Complex128}

    H = StraightLineHomotopy([(x - 4) * (x + 4)], [(x - 2) * (x + 2)])
    @test solve(H, [[4], [-4]], SphericalPredictorCorrector()) isa Result{Complex128}
    @test solve(H, [[4], [-4]], SphericalPredictorCorrector(), CauchyEndgame()) isa Result{Complex128}
    @test solve(H, [[4], [-4]], SphericalPredictorCorrector(), CauchyEndgame(), BigFloat) isa Result{Complex128}



    H, s = randomhomotopy(StraightLineHomotopy, 3)

    @test first(solve(H, [rand(Complex128, 3)])).returncode == :invalid_startvalue
end
