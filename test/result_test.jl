@testset "Result" begin
    @testset "Result+PathResult" begin
        R = solve(equations(heart()), seed=506435)
        @test R isa AffineResult
        @test natinfinity(R) ≤ 572
        @test nfinite(R) == 4
        @test length(collect(R)) == 576
        @test finite(R) isa Vector{<:PathResult}

        @test length(finite(R, onlynonsingular=false)) == 4
        @test length(finite(R, onlynonsingular=true)) == 4
        @test length(finite(R, onlysingular=true)) == 0
        @test 572 - length(failed(R)) == natinfinity(R)
        @test length(real(R, tol=1e-6)) == 2
        @test nreal(R, tol=1e-6) == 2
        @test length(atinfinity(R)) ≤ 572
        @test length(results(R, onlyreal=true, realtol=1e-8)) == 2
        @test length(results(R, onlynonsingular=true, singulartol=1e9)) == 4
        @test length(finite(results(R, onlyreal=true))) == 2
        @test nresults(R, onlynonsingular=true, singulartol=1e9) == 4
        @test length(results(R, onlyfinite=false)) == 576
        @test nresults(R, onlyfinite=false) == 576
        @test nnonsingular(R) == 4
        @test length(nonsingular(R)) == 4
        @test length(singular(R)) == 0
        @test length(singular(real(R))) == 0

        @test_nowarn mapresults(solution, R)
        @test_nowarn mapresults(startsolution, R)
        @test_nowarn mapresults(residual, R)
        @test_nowarn mapresults(issingular, R)
        @test count(mapresults(isreal, R)) == 2
        # test fallback
        @test count(results(isreal, R)) == 2

        @test length(solutions(R)) == 4
        @test solutions(R) isa Vector{Vector{ComplexF64}}
        @test realsolutions(R, realtol=1e-8) isa Vector{Vector{Float64}}
        @test length(realsolutions(R, realtol=1e-8)) == 2

        @test_nowarn string(R)
        @test_nowarn string(R[end])
        @test_nowarn show(IOContext(devnull, :compact=>true), R)
        @test_nowarn show(IOContext(devnull, :compact=>true), R[end])

        test_treeviews(R)
        @test_nowarn TreeViews.treelabel(devnull, R[1], MIME"application/prs.juno.inline"())

        @polyvar x y z
        R = solve([(x-3)^3,(y-2)])
        @test R isa AffineResult
        test_treeviews(R)
        @test_nowarn TreeViews.treelabel(devnull, R[1], MIME"application/prs.juno.inline"())
        @test nnonsingular(R) == 0
        @test nsingular(R) == 3
        @test_nowarn sprint(show, R)

        @polyvar x y z
        R = solve([(x-3z),(y-2z)])
        @test R isa ProjectiveResult
        @test_nowarn string(R)
        test_treeviews(R)
        @test_nowarn TreeViews.treelabel(devnull, R[1], MIME"application/prs.juno.inline"())
        @test nnonsingular(R) == 1

        @polyvar x y z
        R = solve([(x-3z)^3,(y-2z)])
        @test R isa ProjectiveResult
        test_treeviews(R)
        @test_nowarn TreeViews.treelabel(devnull, R[1], MIME"application/prs.juno.inline"())
        @test nnonsingular(R) == 0
        @test nsingular(R) == 3

        @polyvar x y z
        V = x^2 + y^2 - z^2
        L = randn(1,3) * [x, y, z]
        R = solve([V; L])
        @test_nowarn sprint(show, R)
    end

    @testset "multiplicities" begin
        @polyvar x y z
        f = (x-1*y)*(x-1.01*y)
        g = x - y - z
        S = solve([f,g], threading=false)
        @test length(multiplicities(S, tol=1e-1)) == 1
        @test length(multiplicities(S, tol=1e-5)) == 0
        @test length(multiplicities(S.pathresults, tol=1e-1)) == 1
        @test length(multiplicities(S.pathresults, tol=1e-5)) == 0

        @polyvar x z
        f = (x-1)*(x-1.01)
        g = x - 1 - z
        S = solve([f,g], threading=false)
        @test length(multiplicities(S, tol=1e-1)) == 1
        @test length(multiplicities(S, tol=1e-5)) == 0
        @test length(multiplicities(S.pathresults, tol=1e-1)) == 1
        @test length(multiplicities(S.pathresults, tol=1e-5)) == 0
        @test norm(solution(S[1]) - solution(S[2])) < 1e-1
    end

    @testset "uniquesolutions" begin
        @polyvar x
        f = (x-3)^3*(x-2)
        R = solve([f], seed=171090)
        @test length(uniquesolutions(R)) == 2
        @test uniquesolutions(R) isa Vector{Vector{Complex{Float64}}}
        @test length(uniquesolutions(R, multiplicities=true)) == 2
        a, b = uniquesolutions(R, multiplicities=true)
        @test (a[2] == 3 && b[2] == 1) || (b[2] == 3 && a[2] == 1)

        @polyvar x y
        f = (x-3y)^3*(x-2y)
        R = solve([f], seed=171090)
        @test length(uniquesolutions(R)) == 2
        @test uniquesolutions(R) isa Vector{Vector{Complex{Float64}}}
        @test length(uniquesolutions(R, multiplicities=true)) == 2
        a, b = uniquesolutions(R, multiplicities=true)
        @test (a[2] == 3 && b[2] == 1) || (b[2] == 3 && a[2] == 1)
    end
end
