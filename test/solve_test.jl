@testset "solve" begin
    @polyvar x
    @test count(r -> r.returncode == :success, solve([x - 1]).path_results) == 1
    F = equations(katsura5())
    @test count(r -> r.returncode == :success, solve(F, threading=false).path_results) == 32

    @test count(r -> r.returncode == :success, solve(F, homotopy=Homotopies.StraightLineHomotopy)).path_results == 32
    result = solve(F, predictor=Predictors.Euler(), homotopy=Homotopies.StraightLineHomotopy).path_results
    @test count(r -> r.returncode == :success, result) == 32
    result = solve(F,  tol=1e-5).path_results
    @test count(r -> r.returncode == :success, result) == 32

    result = solve(F).path_results
    @test count(r -> r.returncode == :success, result) == 32

    @test string.(result) isa Vector{String}
end

@testset "solve - no endgame" begin
    F = equations(katsura5())
    # no endgame
    @test count(r -> r.returncode == :success, solve(F, endgame_start=0.0).path_results) == 32

    @test count(r -> r.returncode == :success, solve(F, endgame_start=0.0, threading=false).path_results) == 32
end

@testset "Path Crossing" begin
    F = equations(katsura5())
    # this will have two crossed paths
    srand(120)
    @test count(r -> r.returncode == :success, solve(F, tol=1e-1, threading=true).path_results) == 32
    srand(120)
    @test count(r -> r.returncode == :success, solve(F, tol=1e-1, threading=false).path_results) == 32
end
