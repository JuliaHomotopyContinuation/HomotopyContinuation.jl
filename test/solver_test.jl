@testset "Solver" begin
    # I ran into some problems with type instabilities due to Julia internal limits.
    # Hence we should check that everything is properly inferred.
    H, s = randomhomotopy(StraightLineHomotopy, 4)
    @inferred Solver(H, s)

    H, s = randomhomotopy(StraightLineHomotopy{Complex64}, 4)
    @inferred Solver(H, s, BigFloat)


    @test_throws ArgumentError Solver(H, s, wrong_kwarg=true)
end
