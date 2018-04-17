@testset "solve" begin
    F = equations(katsura5())
    @test count(r -> r.returncode == :success, solve(F)) == 32

    @test count(r -> r.returncode == :success, solve(F, homotopy=Homotopies.StraightLineHomotopy)) == 32
    result = solve(F, predictor=Predictors.Euler(), homotopy=Homotopies.StraightLineHomotopy)
    @test count(r -> r.returncode == :success, result) == 32
    result = solve(F,  options=PathTracking.Options(tol=1e-5))
    @test count(r -> r.returncode == :success, result) == 32
end
