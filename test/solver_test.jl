@testset "Solver" begin
    @testset "constructor" begin
        @polyvar x y z
        f = [x^2 - 2, x + y - 1]
        solver, starts = solver_startsolutions(f; system = FPSystem)
        @test isa(solver, Solver)
        result = solve!(solver, starts)
        @test all(is_success, result)
    end

    @testset "path jumping" begin
        solver, starts = solver_startsolutions(
            equations(katsura(5));
            system = FPSystem,
            seed = 124232,
            max_corrector_iters = 5,
            accuracy = 1e-3,
        )
        result = solve!(solver, starts, save_all_paths = false)
        @test all(is_nonsingular, result)
        @test all(is_success, result)
    end
end
