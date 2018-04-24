@testset "solve" begin
    F = equations(katsura5())
    @test count(r -> r.returncode == :success, solve(F)) == 32

    @test count(r -> r.returncode == :success, solve(F, homotopy=Homotopies.StraightLineHomotopy)) == 32
    result = solve(F, predictor=Predictors.Euler(), homotopy=Homotopies.StraightLineHomotopy)
    @test count(r -> r.returncode == :success, result) == 32
    result = solve(F,  tol=1e-5)
    @test count(r -> r.returncode == :success, result) == 32
end

@testset "Path Crossing" begin
    F = equations(katsura5())
    # this will have two crossed paths
    srand(120)
    @test count(r -> r.returncode == :success, solve(F, tol=1e-1)) == 32
end
