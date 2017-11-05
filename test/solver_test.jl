@testset "Solver" begin
    # I ran into some problems with type instabilities due to Julia internal limits.
    # Hence we should check that everything is properly inferred.
    H, s = randomhomotopy(StraightLineHomotopy, 4)
    @inferred Solver(H, s)
    @inferred Solver(H, s, AffinePredictorCorrector())

    H, s = randomhomotopy(StraightLineHomotopy{Complex64}, 4)
    @inferred Solver(H, s, SphericalPredictorCorrector(), CauchyEndgame(), BigFloat)

    @test_throws ArgumentError Solver(H, s, wrong_kwarg=true)

    PolyImpl.@polyvar x
    # default
    @inferred Solver([(x - 2.0) * (x - (2.5+ 4.0im))])

    # just with a single polynomial
    @inferred Solver((x - 2.0) * (x - (2.5+ 4.0im)))


    # non float polynomial, non complex
    @inferred Solver((x - 2) * (x - 4))

    # non float homotopy
    H = StraightLineHomotopy([(x - 4) * (x + (2 - 4im))],[(x - 2) * (x - (2 + 4im))])
    @inferred Solver(H, [[4], [-2 + 4im]])

    # non float homotopy, non complex
    H = StraightLineHomotopy([(x - 4) * (x + 4)], [(x - 2) * (x + 2)])
    @inferred Solver(H, [[4], [-4]])
end
