@testset "PathResult" begin
    R = solve(equations(katsura5()))
    @test R isa Solving.AffineResult
    @test natinfinity(R) == 0
    @test nfinite(R) == 32
    @test finite(R) isa Vector{<:Solving.PathResult}

    @test length(finite(R, include_singular=false)) == 32
    @test string(collect(R)[1]) isa String
    @test isempty(failed(R))
    @test length(real(R, tol=1e-7)) == 12
    @test length(atinfinity(R)) == 0
    @test length(results(R, onlyreal=true, realtol=1e-8)) == 12
    @test length(results(R, includesingular=false, singulartol=1e9)) == 32
    @test length(results(R, includeatinfinity=true)) == 32

    @test_nowarn Juno.render(Juno.Inline(), R)
    @test_nowarn Juno.render(Juno.Inline(), R[1])
end
