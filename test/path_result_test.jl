@testset "Result+PathResult" begin
    R = solve(equations(heart()), seed=9459)
    @test R isa Solving.AffineResult
    @test natinfinity(R) == 572
    @test nfinite(R) == 4
    @test finite(R) isa Vector{<:Solving.PathResult}

    @test length(finite(R, include_singular=false)) == 4
    @test isempty(failed(R))
    @test length(real(R, tol=1e-7)) == 2
    @test length(atinfinity(R)) == 572
    @test length(results(R, onlyreal=true, realtol=1e-8)) == 2
    @test length(results(R, includesingular=false, singulartol=1e9)) == 4
    @test length(results(R, includeatinfinity=true)) == 576

    @test_nowarn results(solution, R)
    @test_nowarn results(start_solution, R)
    @test_nowarn results(residual, R)
    @test_nowarn results(issingular, R)
    @test count(results(isreal, R)) == 2

    @test_nowarn string(R)
    @test_nowarn string(R[end])
end

@testset "PathResult Juno" begin
    R = solve(equations(katsura(5)))
    @test_nowarn Juno.render(Juno.Inline(), R)
    @test_nowarn Juno.render(Juno.Inline(), R[1])
end
