@testset "Solver" begin
    @testset "simple systems" begin
        @polyvar x y z
        f = [x^2 - 2, x + y - 1]
        solver, starts = solver_startsolutions(f; system = FPSystem)
        @test isa(solver, Solver)
        result = solve!(solver, starts)
        @test all(is_success, result)
        @test all(!is_failed, result)

        f = equations(griewank_osborne())
        solver, starts = solver_startsolutions(f; seed = 78373, system = FPSystem)
        result = solve!(solver, starts)
        @test nsolutions(result) == 1
        r = first(results(result))
        @test multiplicity(r) == 3
        @test is_singular(r)
        @test !isempty(sprint(show, r))
    end

    @testset "path jumping" begin
        solver, starts = solver_startsolutions(
            equations(katsura(5));
            system = FPSystem,
            seed = 124232,
            max_corrector_iters = 5,
            accuracy = 1e-3,
        )
        result_jumping = solve!(solver, starts; path_jumping_check = false)
        @test nsolutions(result_jumping) < 32

        result = solve!(solver, starts; path_jumping_check = true)
        @test nsolutions(result) == 32
        @test all(is_nonsingular, result)
        @test all(is_success, result)

        # check that path_jumping_check is on by default
        result2 = solve!(solver, starts)
        @test nsolutions(result2) == 32
    end

    @testset "polyhedral" begin
        solver, starts = solver_startsolutions(
            equations(cyclic(5));
            seed = 32241,
            start_system = :polyhedral,
        )
        result = solve!(solver, starts)
        @test nsolutions(result) == 70
        @test ntracked(result) == 70
        @test nreal(result) == 10

        @polyvar x y
        solver, starts = solver_startsolutions(
            [2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y];
            seed = 32241,
            start_system = :polyhedral,
            only_torus = true,
        )
        result = solve!(solver, starts)
        @test nsolutions(result) == 3
        @test ntracked(result) == 3

        solver, starts = solver_startsolutions(
            [2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y];
            seed = 32241,
            start_system = :polyhedral,
            only_torus = false,
        )
        result = solve!(solver, starts)
        @test nsolutions(result) == 6
        @test ntracked(result) == 8
    end

    @testset "solve" begin
        @polyvar x y z
        f = [x^2 - 2, x + y - 1]
        result = solve(f; system = FPSystem)
        @test all(is_success, result)
        @test result isa Result{Vector{ComplexF64}}
    end
end
