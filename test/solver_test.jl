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
        @test nsolutions(result) == 2
        @test result isa Result{Vector{ComplexF64}}

        @test nsolutions(solve(f; system = FPSystem, predictor = Euler())) == 2
        @test nsolutions(solve(f; system = FPSystem, predictor = RK4())) == 2
    end

    @testset "singular" begin
        @polyvar x y
        z = 1
        F = [
            x^2 + 2 * y^2 + 2 * im * y * z,
            (18 + 3 * im) * x * y + 7 * im * y^2 - (3 - 18 * im) * x * z - 14 * y * z -
            7 * im * z^2,
        ]
        result = solve(F; system = FPSystem)
        @test multiplicity(first(nonsingular(result))) == 1
        @test multiplicity(first(singular(result))) == 3
        @test nsingular(result) == 1
        @test nsingular(result; counting_multiplicities = true) == 3
        @test length(singular(result)) == 1
        @test all(r -> multiplicity(r) == 3, singular(result; multiple_results = true))
    end

    @testset "affine and projective solution types" begin
        @polyvar x y
        result = solve([x - 2y, y - 2z])
        @test result isa Result{PVector{ComplexF64,1}}
        @test first(solutions(result)) isa PVector{ComplexF64,1}

        result2 = solve([x - 2, y - 2])
        @test result2 isa Result{Vector{ComplexF64}}
        @test first(solutions(result2)) isa Vector{ComplexF64}

        @polyvar u v
        result3 = solve([x * y - u * v, x - u], variable_groups = [(x, u), (y, v)])
        @test result3 isa Result{PVector{ComplexF64,2}}
        @test first(solutions(result3)) isa PVector{ComplexF64,2}
    end

    @testset "Keep invalid start values" begin
        @polyvar x[1:2] a[1:2]
        F = [x[1]^2 - a[1], x[1] * x[2] - a[1] + a[2]]
        startsolutions = [[100, 100]]
        p₁ = [1, 1]
        p₀ = [3im, 0.5 + 2im]
        res = solve(
            F,
            startsolutions;
            parameters = a,
            start_parameters = p₁,
            target_parameters = p₀,
        )
        @test nfailed(res) == 1
        @test res[1].return_code == :terminated_invalid_startvalue
    end

    @testset "overdetermined" begin
        @polyvar x y z
        f = [x * z - y^2, y - z^2, x - y * z, x + y + z + 1]

        @test nsolutions(solve(f; seed = 213412)) == 3
        @test nsolutions(solve(f; seed = 213412, start_system = :polyhedral)) == 3

        minors = include(joinpath(@__DIR__, "examples", "3_by_5_minors.jl"))
        @test nsolutions(solve(minors; seed = 312321, system = FPSystem)) == 80

        # singular solutions
        @polyvar x y
        f = [
            (x - 2)^4 * (x + y + 1),
            (x^2 + y^2 - 2) * (y - 2),
            (y - 2) * (x^2 + x * y - 3),
        ]
        result = solve(f)
        @test length(nonsingular(result)) == 1
        @test length(singular(result)) == 1
        @test multiplicity(singular(result)[1]) == 4
    end
end
