@testset "ResultIterator" begin
    @testset "Basic functionality" begin
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

    @testset "Parameters" begin
        @var x y p
        f₁ = y - x^2 + p
        f₂ = y - x^3 - p
        F = System([f₁, f₂]; variables = [x; y], parameters = [p])
        R = solve(F; iterator_only = true, target_parameters = [0])

        @test isa(R, ResultIterator)
        @test nsolutions(R; only_nonsingular = false, multiple_results = true) == 3
        @test nsolutions(R; multiple_results = true) == 1
        @test length(R) == 3
        @test sum(bitmask(isfinite, R)) == 3
    end


end
