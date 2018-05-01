@testset "solve" begin
    F = equations(katsura5())
    @test count(r -> r.returncode == :success, solve(F, threading=false)) == 32

    @test count(r -> r.returncode == :success, solve(F, homotopy=Homotopies.StraightLineHomotopy)) == 32
    result = solve(F, predictor=Predictors.Euler(), homotopy=Homotopies.StraightLineHomotopy)
    @test count(r -> r.returncode == :success, result) == 32
    result = solve(F,  tol=1e-5)
    @test count(r -> r.returncode == :success, result) == 32

    result = solve(F)
    @test count(r -> r.returncode == :success, result) == 32

    @test string.(result) isa Vector{String}
end

@testset "solve - no endgame" begin
    F = equations(katsura5())
    # no endgame
    @test count(r -> r.returncode == :success, solve(F, endgame_start=0.0)) == 32

    @test count(r -> r.returncode == :success, solve(F, endgame_start=0.0, threading=false)) == 32
end

@testset "Path Crossing" begin
    F = equations(katsura5())
    # this will have two crossed paths
    srand(120)
    @test count(r -> r.returncode == :success, solve(F, tol=1e-1, threading=true)) == 32
    srand(120)
    @test count(r -> r.returncode == :success, solve(F, tol=1e-1, threading=false)) == 32
end
