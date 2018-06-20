@testset "Result+PathResult" begin
    R = solve(equations(heart()), seed=9459)
    @test R isa Solving.AffineResult
    @test natinfinity(R) == 572
    @test nfinite(R) == 4
    @test length(collect(R)) == 576
    @test finite(R) isa Vector{<:Solving.PathResult}

    @test length(finite(R, onlynonsingular=false)) == 4
    @test length(finite(R, onlynonsingular=true)) == 4
    @test isempty(failed(R))
    @test length(real(R, tol=1e-7)) == 2
    @test length(atinfinity(R)) == 572
    @test length(results(R, onlyreal=true, realtol=1e-8)) == 2
    @test length(results(R, onlynonsingular=true, singulartol=1e9)) == 4
    @test length(results(R, onlyfinite=false)) == 576
    @test nnonsingular(R) == 4
    @test length(nonsingular(R)) == 4

    @test_nowarn results(solution, R)
    @test_nowarn results(start_solution, R)
    @test_nowarn results(residual, R)
    @test_nowarn results(issingular, R)
    @test count(results(isreal, R)) == 2

    @test length(solutions(R)) == 4
    @test solutions(R) isa Vector{Vector{Complex128}}
    @test solutions(R, Val{true}, realtol=1e-8) isa Vector{Vector{Float64}}
    @test length(solutions(R, Val{true}, realtol=1e-8)) == 2

    @test_nowarn string(R)
    @test_nowarn string(R[end])

    @test_nowarn Juno.render(Juno.Inline(), R)
    @test_nowarn Juno.render(Juno.Inline(), R[1])

    @polyvar x y z
    R = solve([(x-3)^3,(y-2)])
    @test R isa Solving.AffineResult
    @test_nowarn Juno.render(Juno.Inline(), R)
    @test_nowarn Juno.render(Juno.Inline(), R[1])
    @test nnonsingular(R) == 0
    @test nsingular(R) == 3

    @polyvar x y z
    R = solve([(x-3z),(y-2z)])
    @test R isa Solving.ProjectiveResult
    @test_nowarn string(R)
    @test_nowarn Juno.render(Juno.Inline(), R)
    @test_nowarn Juno.render(Juno.Inline(), R[1])
    @test nnonsingular(R) == 1

    @polyvar x y z
    R = solve([(x-3z)^3,(y-2z)])
    @test R isa Solving.ProjectiveResult
    @test_nowarn Juno.render(Juno.Inline(), R)
    @test_nowarn Juno.render(Juno.Inline(), R[1])
    @test nnonsingular(R) == 0
    @test nsingular(R) == 3
end
