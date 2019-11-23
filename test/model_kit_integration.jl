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
    end

    @testset "polyhedral" begin
        # converts into polynomial in the background
        @var x y
        result = solve(
            [2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y];
            variable_ordering = [x, y],
            seed = 32241,
            start_system = :polyhedral,
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

        f = System([x^2 - a, y^2 - b], [x, y, z], [a, b])
        res = solve(
            f,
            [[1, -1], [1, 1], [-1, 1], [-1, -1]];
            start_parameters = [1, 1],
            target_parameters = [3, 4],
        )
        @test nonsingular(res) == 4

        f = System([x^2 - z^2 * a, y^2 - z^2 * b], [x, y, z], [a, b])
        res = solve(
            f,
            [[1, -1], [1, 1], [-1, 1], [-1, -1]];
            start_parameters = [1, 1],
            target_parameters = [3, 4],
        )
        @test nonsingular(res) == 4
    end
end
