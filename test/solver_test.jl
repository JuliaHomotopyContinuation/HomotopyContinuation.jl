@testset "Solver" begin
    @testset "simple systems" begin
        @polyvar x y z
        f = [x^2 - 2, x + y - 1]
        solver, starts = solver_startsolutions(f; system = FPSystem)
        @test isa(solver, Solver)
        result = solve(solver, starts)
        @test all(is_success, result)
        @test all(!is_failed, result)

        f = equations(griewank_osborne())
        solver, starts = solver_startsolutions(f; seed = 78373, system = FPSystem)
        result = solve(solver, starts)
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
            max_corrector_iters = 6,
            accuracy = 1e-5,
        )
        result_jumping = solve(solver, starts; path_jumping_check = false)
        @test nsolutions(result_jumping) < 32

        for threading in [false, true]
            result = solve(solver, starts; path_jumping_check = true, threading = threading)
            @test nsolutions(result) == 32
            @test all(is_nonsingular, result)
            @test all(is_success, result)
        end

        # check that path_jumping_check is on by default
        result2 = solve(solver, starts)
        @test nsolutions(result2) == 32
    end

    @testset "polyhedral" begin
        solver, starts = solver_startsolutions(
            equations(cyclic(5));
            seed = 32241,
            start_system = :polyhedral,
        )
        result = solve(solver, starts)
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
        result = solve(solver, starts, threading = false)
        @test nsolutions(result) == 3
        @test ntracked(result) == 3

        solver, starts = solver_startsolutions(
            [2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y];
            seed = 32241,
            start_system = :polyhedral,
            only_torus = false,
        )
        result = solve(solver, starts)
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

    @testset "different input" begin
        F = equations(katsura(5))
        F_hom = homogenize(F)
        res = solve(F_hom, threading = false, patch = RandomPatch())
        @test nfinite(res) == 32
        @test res isa Result{PVector{ComplexF64,1}}
        @test nfinite(solve(F_hom, threading = false, patch = EmbeddingPatch())) == 32
        @test nfinite(solve(F_hom, threading = false, patch = OrthogonalPatch())) == 32

        @polyvar w
        F = equations(cyclic(5))
        result = solve(homogenize(F, w), system = FPSystem, threading = false, homvar = w)
        @test nfinite(result) == 70
        @test result isa Result{Vector{ComplexF64}}
        G = FPSystem(homogenize(F))
        result = solve(G, homvar = 6)
        @test nnonsingular(result) == 70
        @test result isa Result{Vector{ComplexF64}}
    end

    @testset "singular" begin
        # 1 singular solution with multiplicity 3
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

        # 2 solutions with multiplicity 6, projective
        @polyvar x z y
        # This has two roots of multiplicity 6 at the hyperplane z=0
        F = [
            0.75 * x^4 + 1.5 * x^2 * y^2 - 2.5 * x^2 * z^2 + 0.75 * y^4 - 2.5 * y^2 * z^2 +
            0.75 * z^4
            10 * x^2 * z + 10 * y^2 * z - 6 * z^3
        ]
        res = solve(F; system = FPSystem, threading = false, seed = 703127)
        @test nsingular(res) == 2
        @test nsingular(res; counting_multiplicities = true) == 12
        @test multiplicity.(results(res)) == [6, 6]

        # 2 solutions with multiplicity 6, affine
        @polyvar x z
        y = 1
        # This has two roots of multiplicity 6 at the hyperplane z=0
        F = [
            0.75 * x^4 + 1.5 * x^2 * y^2 - 2.5 * x^2 * z^2 + 0.75 * y^4 - 2.5 * y^2 * z^2 +
            0.75 * z^4
            10 * x^2 * z + 10 * y^2 * z - 6 * z^3
        ]
        res = solve(F; system = FPSystem, threading = false, seed = 412312)
        @test nsingular(res) == 2
        @test nsingular(res; counting_multiplicities = true) == 12
        @test multiplicity.(results(res)) == [6, 6]
    end

    @testset "affine and projective solution types" begin
        @polyvar x y z
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
        @test nsolutions(solve(
            f;
            seed = 213412,
            start_system = :polyhedral,
            threading = false,
        )) == 3
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

    @testset "Parameter promotion to Float64" begin
        d = 3
        k = 1

        R = 1.0
        r = 0.5

        @polyvar lam[1:k] p[1:6]
        F = [(R^2 - r^2 + (p[1] + p[4] * lam[1])^2 + (p[2] + p[5] * lam[1])^2 +
              (p[3] + p[6] * lam[1])^2)^2 -
             4 * R^2 * ((p[1] + p[4] * lam[1])^2 + (p[2] + p[5] * lam[1])^2)]
        #Randomly choose a start system
        # with only Complex{Float32 parameters
        p0 = randn(Complex{Float32}, 6)
        F0 = subs(F, p => p0)
        res0 = solve(F0)
        S0 = solutions(res0)

        #Construct the tracker
        tracker = pathtracker(F; parameters = p, generic_parameters = p0)
        @test HC.basehomotopy(tracker.core_tracker.homotopy).p isa NTuple{
            2,
            Vector{ComplexF64},
        }

        #Let's solve another system using PathTracking
        p_current = randn(6)
        result = track(tracker, S0[1]; target_parameters = p_current)
        sol = solution(result)
        F_SP = SPSystem(F; parameters = p)
        @test residual(result) ≈ euclidean_norm(F_SP(sol, p_current)) atol = 1e-13
    end

    @testset "early stop callback" begin
        @polyvar x
        first_result = nothing
        results = solve(
            [(x - 3) * (x + 6) * (x + 2)];
            system = FPSystem,
            stop_early_cb = r -> begin
                first_result = r
                true
            end,
            threading = false,
        )
        @test length(results) == 1
        @test first(results) == first_result

        nresults = 0
        result = solve(
            [(x - 3) * (x + 6) * (x + 2)];
            system = FPSystem,
            stop_early_cb = r -> (nresults += 1) == 2,
            threading = false,
        )
        @test length(result) == 2


        # threading
        @polyvar x y z
        first_result = nothing
        # this has 5^3 = 125 solutions, so we should definitely stop early if we
        # have less than 64 threads
        results = solve(
            [
             (x - 3) * (x + 6) * (x + 2) * (x - 2) * (x + 2.5),
             (y + 2) * (y - 2) * (y + 3) * (y + 5) * (y - 1),
             (z + 2) * (z - 2) * (z + 3) * (z + 5) * (z - 2.1),
            ];
            system = FPSystem,
            stop_early_cb = r -> begin
                first_result = r
                true
            end,
            threading = true,
        )
        @test length(results) < 125
    end

    @testset "variable ordering" begin
        # total degree
        @polyvar x y z
        f = [x + 1, y + 2, z + 3]
        @test solution(first(solve(f))) ≈ [-1, -2, -3]
        @test solution(first(solve(f; variable_ordering = [z, y, x]))) ≈ [-3, -2, -1]
        @test_throws ArgumentError solve(f; variable_ordering = [z, y])

        # parameter
        @polyvar x a y b
        F = [x^2 - a, x * y - a + b]
        p = [a, b]

        r1 = solve(F, [√2, √2]; parameters = p, p₁ = [2, 0], p₀ = [2, 4])
        @test is_success(first(r1))
        s1 = solution(first(r1))

        r2 = solve(
            F,
            [√2, √2];
            variable_ordering = [y, x],
            parameters = p,
            p₁ = [2, 0],
            p₀ = [2, 4],
        )
        s2 = solution(first(r2))

        @test s1 ≈ s2[[2, 1]]
        @test_throws ArgumentError solve(
            F,
            [√2, √2];
            variable_ordering = [y],
            parameters = p,
            p₁ = [2, 0],
            p₀ = [2, 4],
        )

        # start target
        @polyvar x y
        G = [x^2 - 2, x * y - 2]
        F = [x^2 - 2, x * y + 2]
        r3 = solve(G, F, [√2, √2]; variable_ordering = [y, x])
        s3 = solution(first(r3))
        @test s1 ≈ s3[[2, 1]]

        @test_throws ArgumentError solve(G, F, [√2, √2]; variable_ordering = [x])
    end

    @testset "parameter change" begin
        @polyvar x y z p[1:3]
        F = [
            x + 3 + 2y + 2 * y^2 - p[1],
            (x - 2 + 5y) * z + 4 - p[2] * z,
            (x + 2 + 4y) * z + 5 - p[3] * z,
        ]
        q = randn(ComplexF64, 3)
        S = solutions(solve(subs(F, p => q)))
        # create some fake parameter values
        params = [
            [1.3440639843905677, -2.5074017514201343, -0.39441026816083785],
            [-0.5905470689391262, 0.662942297773609, -0.4383381251082712],
            [-0.07181848550172094, -0.9892961521599329, 0.3648603261873276],
            [2.37916274514681, -0.4796116819359556, 1.4850511902292287],
            [1.0080243329063356, 0.1156613479279384, -0.4496252569645542],
            [0.5658877469626858, 0.2383098725793224, -0.13845976523959966],
            [-0.22823657335108488, -0.7111958812641874, -0.7478056668942166],
            [-0.255059563703982, -0.9004822204456397, -0.3389293738788412],
            [0.3392911355486216, 0.05817061767333153, 0.33923705766735207],
            [1.0861653811709309, -0.31406883615041653, -0.17861101511219216],
        ]
        # create a `Solver` to reuse for the path tracking
        F_solver = solver(F; parameters = p, generic_parameters = q)
        # solve the system F for all paramaters p in params
        params_solutions = map(params) do p
            solutions(solve(F_solver, S; target_parameters = p))
        end
        @test sum(length, params_solutions) == 20
    end

    @testset "Many parameters solver" begin
        # Setup
        @polyvar x y
        f = x^2 + y^2 - 1

        @polyvar a b c
        l = a * x + b * y + c
        F = [f, l]

        # Compute start solutions S₀ for given start parameters p₀
        p₀ = randn(ComplexF64, 3)
        S₀ = solutions(solve(subs(F, [a, b, c] => p₀)))
        # The parameters we are intersted in
        params = [rand(3) for i = 1:100]

        result1 = solve(
            F,
            S₀,
            ;
            start_parameters = p₀,
            target_parameters = params,
            parameters = [a, b, c],
            threading = true,
        )
        @test typeof(result1) == Vector{Tuple{Result{Vector{ComplexF64}},Vector{Float64}}}
        result1 = solve(
            F,
            S₀,
            ;
            start_parameters = p₀,
            target_parameters = params,
            parameters = [a, b, c],
            threading = false,
        )
        @test typeof(result1) == Vector{Tuple{Result{Vector{ComplexF64}},Vector{Float64}}}

        # Only keep real solutions
        result2 = solve(
            F,
            S₀,
            ;
            start_parameters = p₀,
            target_parameters = params,
            parameters = [a, b, c],
            transform_result = (r, p) -> real_solutions(r),
            threading = true,
        )
        @test typeof(result2) == Vector{Vector{Vector{Float64}}}
        @test !isempty(result2)

        # Now instead of an Array{Array{Array{Float64,1},1},1} we want to have an
        # Array{Array{Float64,1},1}
        result3 = solve(
            F,
            S₀,
            ;
            start_parameters = p₀,
            target_parameters = params,
            parameters = [a, b, c],
            transform_result = (r, p) -> real_solutions(r),
            flatten = true,
            threading = false,
        )
        @test typeof(result3) == Vector{Vector{Float64}}
        @test !isempty(result3)

        # The passed `params` do not directly need to be the target parameters.
        # Instead they can be some more concrete informations (e.g. an index)
        # and we can them by using the `transform_parameters` method
        result4 = solve(
            F,
            S₀,
            ;
            start_parameters = p₀,
            target_parameters = 1:100,
            parameters = [a, b, c],
            transform_result = (r, p) -> (real_solutions(r), p),
            transform_parameters = _ -> rand(3),
        )
        @test typeof(result4) == Vector{Tuple{Vector{Vector{Float64}},Int64}}
    end
end
