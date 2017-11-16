@testset "Solver" begin
    # I ran into some problems with type instabilities due to Julia internal limits.
    # Hence we should check that everything is properly inferred.
    H, s = randomhomotopy(StraightLineHomotopy, 4)
    @inferred Solver(H)
    @inferred Solver(H, AffinePredictorCorrector())

    H, s = randomhomotopy(StraightLineHomotopy{Complex64}, 4)
    @inferred Solver(H, SphericalPredictorCorrector(), CauchyEndgame(), BigFloat)

    @test_throws ArgumentError Solver(H, wrong_kwarg=true)

    PolyImpl.@polyvar x
    # default
    @inferred solve([(x - 2.0) * (x - (2.5+ 4.0im))])

    # just with a single polynomial
    @inferred solve((x - 2.0) * (x - (2.5+ 4.0im)))


    # non float polynomial, non complex
    @inferred solve((x - 2) * (x - 4))

    # non float homotopy
    H = StraightLineHomotopy([(x - 4) * (x + (2 - 4im))],[(x - 2) * (x - (2 + 4im))])
    @inferred solve(H, [[4], [-2 + 4im]])

    # non float homotopy, non complex
    H = StraightLineHomotopy([(x - 4) * (x + 4)], [(x - 2) * (x + 2)])
    @inferred solve(H, [[4], [-4]])
end
