@testset "Result" begin
    @testset "Basic functionality of Result" begin
        d = 2
        @var x y a[1:6]
        F = System(
            [
                (a[1] * x^d + a[2] * y) * (a[3] * x + a[4] * y) + 1,
                (a[1] * x^d + a[2] * y) * (a[5] * x + a[6] * y) + 1,
            ];
            parameters = a,
        )
        res = solve(F; target_parameters = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32])

        @test startswith(sprint(show, res), "Result with 3 solutions")
        @test seed(res) isa UInt32
        test_treeviews(res)

        seeded_res = solve(
            F;
            target_parameters = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32],
            seed = seed(res),
        )
        @test seed(seeded_res) == seed(res)
        test_treeviews(res)

        @test length(path_results(res)) == ntracked(res) == 7
        @test length(results(res)) == nresults(res) == 3
        @test length(solutions(res)) == 3
        @test length(findall(is_success, res)) == 3
        @test real_solutions(res) isa Vector{Vector{Float64}}
        @test length(real_solutions(res)) == nreal(res) == 1
        @test length(nonsingular(res)) == nnonsingular(res) == 3
        @test isempty(singular(res))
        @test nsingular(res) == 0
        # @test length(at_infinity(res)) == nat_infinity(res) == 4
        # @test isempty(failed(res))
        # @test nfailed(res) == 0
        @test nexcess_solutions(res) == 0
        @test !isempty(sprint(show, statistics(res)))

        # singular result
        @var x y
        g = System([(29 / 16) * x^3 - 2 * x * y, x^2 - y])
        res = solve(g; start_system = :total_degree)
        @test startswith(sprint(show, res), "Result with 1 solution")
        @test seed(res) isa UInt32
        test_treeviews(res)
        @test !isempty(sprint(show, statistics(res)))
        @test !isempty(sprint(show, res))
    end
end

@testset "ResultIterator" begin

    @testset "Result tests" begin
        d = 2
        @var x y a[1:6]
        F = System(
            [
                (a[1] * x^d + a[2] * y) * (a[3] * x + a[4] * y) + 1,
                (a[1] * x^d + a[2] * y) * (a[5] * x + a[6] * y) + 1,
            ];
            parameters = a,
        )
        res = solve(
            F;
            iterator_only = true,
            target_parameters = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32],
        )

        @test startswith(sprint(show, res), "ResultIterator")
        @test seed(res) isa UInt32

        @test length(path_results(res)) == ntracked(res) == 7
        @test length(results(res)) == nresults(res) == 3
        @test length(solutions(res)) == 3
        @test length(findall(is_success, res)) == 3
        @test real_solutions(res) isa Vector{Vector{Float64}}
        @test length(real_solutions(res)) == nreal(res) == 1
        @test length(nonsingular(res)) == nnonsingular(res) == 3
        @test isempty(singular(res))
        @test nsingular(res) == 0
        @test nexcess_solutions(res) == 0

        # singular result
        @var x y
        g = System([(29 / 16) * x^3 - 2 * x * y, x^2 - y])
        res = solve(g; start_system = :total_degree, iterator_only = true)
        @test startswith(sprint(show, res), "ResultIterator")
        @test seed(res) isa UInt32
        @test !isempty(sprint(show, res))
    end

    @testset "Basic functionality of ResultIterator" begin
        @var x y
        # define the polynomials
        f₁ = y-x^2
        f₂ = y-x^3
        F = [f₁, f₂]
        tsi_polyhedral = solve(F; iterator_only = true, start_system = :polyhedral)
        tsi_total_degree = solve(F; iterator_only = true, start_system = :total_degree)

        @test isa(tsi_polyhedral, ResultIterator)
        @test isa(tsi_total_degree, ResultIterator)

        @test nsolutions(
            tsi_polyhedral;
            only_nonsingular = false,
            multiple_results = true,
        ) == 3
        @test nsolutions(tsi_total_degree; multiple_results = true) == 1
        @test length(tsi_polyhedral) == 3
        @test length(tsi_total_degree) == 6
        @test sum(bitmask(isfinite, tsi_total_degree)) == 3
    end

    @testset "Manual start solutions" begin
        @var x y p
        f₁ = y - x^2 + p
        f₂ = y - x^3 - p
        F = System([f₁, f₂]; variables = [x; y], parameters = [p])
        R = solve(
            F,
            [1, -1];
            iterator_only = true,
            start_parameters = [0],
            target_parameters = [-1],
        )

        @test isa(R, ResultIterator)
        @test nsolutions(R) == 1

        R2 = solve(
            F,
            [[1, -1], [1, -1]];
            iterator_only = true,
            start_parameters = [0],
            target_parameters = [-1],
        )

        @test isa(R2, ResultIterator)
        @test nsolutions(R2) == 2
    end
end
