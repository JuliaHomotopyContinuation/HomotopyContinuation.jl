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
        @test length(real_solutions(res; atol = 1e-16)) == 1
        @test length(real_solutions(res; atol = 0.0)) == 0
        @test length(real_solutions(res; atol = 0.0, rtol = 1e-16)) == 1
        @test length(real_solutions(res; atol = Inf, rtol = Inf)) == 3
        @test length(real_solutions(res; atol = 1.0, rtol = 1e-16)) == 1
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
        param = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32]
        res = solve(F; iterator_only = true, target_parameters = param)

        @test startswith(sprint(show, res), "ResultIterator")
        @test seed(res) isa UInt32

        @test length(path_results(res)) == ntracked(res) == 7
        @test nresults(res) == 3
        @test nsolutions(res) == 3
        real_sols = collect(real_solutions(res))
        @test length(real_sols) == nreal(res) == 1
        @test nnonsingular(res) == 3
        @test isempty(singular(res))
        @test nsingular(res) == 0
        @test nexcess_solutions(res) == 0

        # pass bit-vector 
        B = BitVector([1, 0, 0, 0, 0, 0, 0])
        res_B = solve(F; iterator_only = true, bitmask = B, target_parameters = param)
        @test length(res_B) == 1

        # singular result
        @var x y
        g = System([(29 / 16) * x^3 - 2 * x * y, x^2 - y])
        res = solve(g; start_system = :total_degree, iterator_only = true)
        @test startswith(sprint(show, res), "ResultIterator")
        @test seed(res) isa UInt32
        @test !isempty(sprint(show, res))


        # pass iterator as start solutions
        @var x y p
        f₁ = y - x^2 + p
        f₂ = y - x^3 - p
        F = System([f₁, f₂]; variables = [x; y], parameters = [p])
        R = solve(
            F,
            [[1, 1], [-1, 1]];
            iterator_only = true,
            start_parameters = [0],
            target_parameters = [-1],
        )
        RR = solve(
            F,
            R;
            iterator_only = true,
            start_parameters = [-1],
            target_parameters = [-2],
        )
        @test length(collect(RR)) == 2
    end

    @testset "Basic functionality of ResultIterator" begin
        @var x y
        # define the polynomials
        f₁ = y - x^2
        f₂ = y - x^3
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

        BM = bitmask_filter(isfinite, tsi_total_degree)
        @test length(BM) == sum(bitmask(isfinite, tsi_total_degree)) == 3

        t = trace(BM)
        @test norm([1.0 + 0.0im, 1.0 + 0.0im] - t) < 1e-12
    end

    @testset "Manual start solutions" begin
        @var x y p
        f₁ = y - x^2 + p
        f₂ = y - x^3 - p
        F = System([f₁, f₂]; variables = [x; y], parameters = [p])
        R = solve(
            F,
            [1, 1];
            iterator_only = true,
            start_parameters = [0],
            target_parameters = [-1],
        )

        @test isa(R, ResultIterator)
        @test nsolutions(R) == 1
        nsolutions(R)

        R2 = solve(
            F,
            [[1, 1], [1, 1]];
            iterator_only = true,
            start_parameters = [0],
            target_parameters = [-1],
        )

        @test isa(R2, ResultIterator)
        @test nsolutions(R2) == 2
    end

    @testset "Many parameters" begin
        ## these tests are the same in in solve_test.jl
        @var x y
        f = x^2 + y^2 - 1

        @var a b c
        l = a * x + b * y + c
        F = [f, l]

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
            iterator_only = true,
        )
        r1 = collect.(result1)
        @test !isempty(r1)

        result1 = solve(
            F,
            S₀,
            ;
            start_parameters = p₀,
            target_parameters = params,
            parameters = [a, b, c],
            show_progress = false,
            threading = false,
            iterator_only = true,
        )
        r1 = collect.(result1)
        @test !isempty(r1)

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
            iterator_only = true,
        )
        r2 = collect.(result2)
        @test !isempty(r2)

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
            iterator_only = true,
        )
        r3 = collect.(result3)
        @test !isempty(r3)

        # The passed `params` do not directly need to be the target parameters.
        # Instead they can be some more concrete informations (e.g. an index)
        result4 = solve(
            F,
            S₀,
            ;
            start_parameters = p₀,
            target_parameters = 1:100,
            parameters = [a, b, c],
            transform_result = (r, p) -> (real_solutions(r), p),
            transform_parameters = _ -> rand(3),
            iterator_only = true,
        )
        r4 = collect.(result4)
        @test !isempty(r4)
    end

    @testset "Target subspaces" begin
        @var x y
        F = System([x^2 + y^2 - 5], [x, y])
        l1 = rand_subspace(2; codim = 1)
        l2 = rand_subspace(2; codim = 1)

        r1 = solve(F; target_subspace = l1, intrinsic = true, iterator_only = true)

        w1 = collect(r1)
        @test length(w1) == 2
        @test w1 isa Vector{PathResult}

        r2 = solve(F, r1; start_subspace = l1, target_subspace = l2, iterator_only = true)

        w2 = collect(r2)
        @test length(w2) == 2
        pt = solution(w2[1])
        @test norm(l2(pt)) < 1e-12
    end

    @testset "Compression" begin
        F = cyclic(5)
        R = solve(F)
        # this function encodes the zeros R of F as an iterator with total degree start system 
        function compress(F, R; gamma = cis(rand() * 2pi))
            d = degrees(F)
            v = variables(F)
            k = length(R)
            S = total_degree_start_solutions(d)
            G = System(gamma .* [vi^di - 1 for (vi, di) in zip(v, d)], variables = v)
            T = solutions(solve(F, G, R))

            U = UniquePoints(first(T), 0)
            for (i, t) in enumerate(T)
                add!(U, t, i)
            end
            B = BitVector()
            for (i, s) in enumerate(S)
                _, j = add!(U, s, i + k)
                push!(B, !j)
            end

            solve(G, F, S; iterator_only = true, bitmask = B)
        end

        C = compress(F, R)

        @test C isa ResultIterator
        R_test = collect(C)
        @test length(R_test) == 70
        @test all(is_success, R_test)
    end
end
