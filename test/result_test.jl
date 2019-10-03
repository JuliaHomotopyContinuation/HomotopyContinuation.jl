@testset "Result" begin
    @testset "Result+PathResult" begin
        R = solve(
            equations(heart()),
            save_all_paths = true,
            seed = 13186,
            show_progress = false,
        )
        @test R isa Result
        @test nat_infinity(R) ≤ 572
        @test nfinite(R) == 4
        @test length(collect(R)) == 576
        @test finite(R) isa Vector{<:PathResult}
        @test seed(R) == 13186

        @test length(finite(R, only_nonsingular = false)) == 4
        @test length(finite(R, only_nonsingular = true)) == 4
        @test length(finite(R, only_singular = true)) == 0
        @test 572 - length(failed(R)) == nat_infinity(R)
        @test length(real(R, tol = 1e-6)) == 2
        @test nreal(R, tol = 1e-6) == 2
        @test length(real_solutions(R)) == 2
        @test_deprecated realsolutions(R)
        @test length(at_infinity(R)) ≤ 572
        @test length(results(R, only_real = true, real_tol = 1e-8)) == 2
        @test length(results(R, only_nonsingular = true, singular_tol = 1e9)) == 4
        @test length(finite(results(R, only_real = true))) == 2
        @test nresults(R, only_nonsingular = true, singular_tol = 1e9) == 4
        @test length(results(R, onlyfinite = false)) == 576
        @test nresults(R, onlyfinite = false) == 576
        @test nnonsingular(R) == 4
        @test length(nonsingular(R)) == 4
        @test length(singular(R)) == 0
        @test length(singular(real(R))) == 0

        @test_nowarn mapresults(solution, R)
        @test_nowarn mapresults(start_solution, R)
        @test_nowarn mapresults(residual, R)
        @test_nowarn mapresults(is_singular, R)
        @test count(mapresults(is_real, R)) == 2
        # test fallback
        @test count(results(is_real, R)) == 2

        @test length(solutions(R)) == 4
        @test solutions(R) isa Vector{Vector{ComplexF64}}
        @test real_solutions(R, real_tol = 1e-8) isa Vector{Vector{Float64}}
        @test length(real_solutions(R, real_tol = 1e-8)) == 2

        @test_nowarn string(R)
        @test_nowarn string(R[end])
        @test_nowarn show(IOContext(devnull, :compact => true), R)
        @test_nowarn show(IOContext(devnull, :compact => true), R[end])

        test_treeviews(R)

        @polyvar x y z
        R = solve([(x - 3)^3, (y - 2)], affine_tracking = true)
        @test R isa Result
        test_treeviews(R)
        @test nnonsingular(R) == 0
        @test nsingular(R) == 1
        @test nsingular(R, counting_multiplicities = true) == 3
        @test_nowarn sprint(show, R)

        R = solve([(x - 3)^3, (y - 2)], affine_tracking = true, system_scaling = nothing)
        @test nnonsingular(R) == 0
        @test nsingular(R) == 1
        @test nsingular(R, counting_multiplicities = true) == 3
        @test_nowarn sprint(show, R)

        @polyvar x y z
        R = solve([(x - 3z), (y - 2z)])
        @test R isa Result{<:ProjectiveVectors.PVector}
        @test_nowarn string(R)
        test_treeviews(R)
        @test nnonsingular(R) == 1

        @polyvar x y z
        R = solve([(x - 3z)^3, (y - 2z)])
        @test R isa Result{<:ProjectiveVectors.PVector}
        test_treeviews(R)
        @test nnonsingular(R) == 0
        @test nsingular(R) == 1
        @test nsingular(R, counting_multiplicities = true) == 3

        @polyvar x y z
        V = x^2 + y^2 - z^2
        L = randn(1, 3) * [x, y, z]
        R = solve([V; L])
        @test_nowarn sprint(show, R)
    end

    @testset "singular result" begin
        @polyvar x
        f = (x - 3)^3 * (x - 2)
        R = solve([f])
        @test length(solutions(R)) == 2
        @test sort(results(multiplicity, R)) == [1, 3]

        @polyvar x y
        f = (x - 3y)^3 * (x - 2y)
        R = solve([f])
        @test nsolutions(R) == 2
        @test sort(results(multiplicity, R)) == [1, 3]
    end

    @testset "multiplicities" begin
        @polyvar x y z
        f = (x - 1 * y) * (x - 1.01 * y)
        g = x - y - z
        S = solve([f, g], threading = false)
        @test length(multiplicities(S, tol = 1e-1)) == 1
        @test length(multiplicities(S, tol = 1e-5)) == 0
        @test length(multiplicities(S.pathresults, tol = 1e-1)) == 1
        @test length(multiplicities(S.pathresults, tol = 1e-5)) == 0

        @polyvar x z
        f = (x - 1) * (x - 1.01)
        g = x - 1 - z
        S = solve([f, g], threading = false)
        @test length(multiplicities(S, tol = 1e-1)) == 1
        @test length(multiplicities(S, tol = 1e-5)) == 0
        @test length(multiplicities(S.pathresults, tol = 1e-1)) == 1
        @test length(multiplicities(S.pathresults, tol = 1e-5)) == 0
        @test norm(solution(S[1]) - solution(S[2])) < 1e-1
    end
end
