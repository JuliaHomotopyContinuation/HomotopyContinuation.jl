@testset "ModelKit integration" begin
    @testset "total degree" begin
        @var x y
        f = [2 * x^2 + 3y + 5, x * y + 2 * y^2 + 2x + 3]
        @test_throws ArgumentError solve(f)
        res = solve(f; variable_ordering = [x, y])
        @test nnonsingular(res) == 4

        res = solve(System(f, [x, y]))
        @test nnonsingular(res) == 4

        # overdetermined
        g = [f; sum(f)]
        res = solve(System(g, [x, y]))
        @test nnonsingular(res) == 4

        @var x y z
        @test_throws ArgumentError problem_startsolutions(
            [x - 2y, y^2 + 3 * x^2, x^3 + y^3],
            variable_ordering = [x, y],
        )

        prob, starts = problem_startsolutions(
            [x^2 - 2 * z^2, y^2 + 3 * z^2, z^2 + x^2, y * z + x^2],
            variable_ordering = [x, y, z],
        )
        @test prob isa HC.OverdeterminedProblem{HC.ProjectiveTracking}
        @test starts isa HC.TotalDegreeSolutionIterator
        @test starts.degrees == [2, 2]
    end

    @testset "polyhedral" begin
        # converts into polynomial in the background
        @var x y
        result = solve(
            [2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y];
            variable_ordering = [x, y], seed = 32241, start_system = :polyhedral,
        )
        @test nsolutions(result) == 6

        result = solve(
            System([2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y], [x, y]),
            seed = 32241,
            start_system = :polyhedral,
        )
        @test nsolutions(result) == 6
    end


    @testset "parameter homotopy" begin
        @var x y z a b

        f = System([x^2 - a, y^2 - b], [x, y], [a, b])
        res = solve(
            f,
            [[1, -1], [1, 1], [-1, 1], [-1, -1]];
            start_parameters = [1, 1], target_parameters = [3, 4],
        )
        @test nnonsingular(res) == 4

        f = System([x^2 - z^2 * a, y^2 - z^2 * b], [x, y, z], [a, b])
        res = solve(
            f,
            [[1, -1], [1, 1], [-1, 1], [-1, -1]];
            start_parameters = [1, 1], target_parameters = [3, 4],
        )
        @test nnonsingular(res) == 4


        f = System([x^2 - a, y^2 - b], [x, y], [a, b])
        res = solve(
            f,
            [[1, -1], [1, 1], [-1, 1], [-1, -1]];
            start_parameters = [1, 1],
            target_parameters = [3, 4],
            γ₁ = cis(2.132im),
            γ₀ = cis(0.412im),
        )
        @test nnonsingular(res) == 4
    end

    @testset "Monodromy" begin
        @var x y
        f = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        # define new variables u₁, u₂ and λ₁
        @var u[1:2] λ[1:1]
        # define the jacobian of F
        J = differentiate([f], [x, y])
        # J' defines the transpose of J
        C = [[x, y] - u - J' * λ; f]
        F = System(C, [x; y; λ], u)

        s₀ = [
            0.3785092033935262 + 0.20148856470590054im,
            0.008762115289894847 + 1.0025348274195496im,
            -0.08344446740253277 + 0.01612603621102832im,
        ]
        u₀ = [
            0.3373116901371924 + 0.12535355543080084im,
            0.27964298581644154 + 1.9612473908330865im,
        ]
        @test_throws ArgumentError monodromy_solve(C, s₀, u₀)
        @test_throws ArgumentError monodromy_solve(C, s₀, u₀; variable_ordering = [x; y; λ])

        res = monodromy_solve(
            C,
            s₀,
            u₀;
            variable_ordering = [x; y; λ],
            parameters = u,
            seed = 1312,
            target_solutions_count = 36,
        )
        @test nsolutions(res) == 36
        res = monodromy_solve(
            F,
            s₀,
            u₀;
            threading = false, target_solutions_count = 36, seed = 1312,
        )
        @test nsolutions(res) == 36

        # automatic start pair
        @var x y a b
        f = [x^2 - a, y^2 - b]
        res =
            monodromy_solve(f; parameters = [a, b], variable_ordering = [x, y], seed = 1312)
        @test nsolutions(res) == 4
        F = System(f, [x, y], [a, b])
        res = monodromy_solve(F; seed = 1312)
        @test nsolutions(res) == 4
    end
end
