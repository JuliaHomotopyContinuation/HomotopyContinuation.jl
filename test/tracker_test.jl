@testset "Tracker" begin
    @testset "tracking $mode - AD: $AD" for mode in [InterpretedSystem, CompiledSystem],
        AD = 0:3

        @var x a y b
        F = System([x^2 - a, x * y - a + b], [x, y], [a, b])

        tracker = Tracker(
            ParameterHomotopy(mode(F), [1, 0], [2, 4]),
            options = TrackerOptions(automatic_differentiation = AD),
        )
        s = [1, 1]
        res = track(tracker, s, 1, 0)
        @test is_success(res)
        @test steps(res) ≤ 5
        @test isa(solution(res), Vector{ComplexF64})
        @test solution(res) ≈ [sqrt(2), -sqrt(2)]
        @test length(solution(res)) == 2
        @test is_success(track(tracker, res, 0, 1))

        code = track!(tracker, s, 1, 0)
        @test is_success(code)
        s0 = copy(tracker.state.x)
        @unpack μ, ω = tracker.state
        @test is_success(track(tracker, s0, 0, 1))
        @test is_success(track(tracker, s0, 0, 1, μ = μ, ω = ω))
    end

    @testset "iterator" begin
        @var x a y b
        F = System([x^2 - a, x * y - a + b], [x, y], [a, b])
        tracker = Tracker(ParameterHomotopy(F, [1, 0], [2, 4]))
        s = [1, 1]

        # path iterator
        typeof(first(iterator(tracker, s, 1.0, 0.0))) == Tuple{Vector{ComplexF64},Float64}

        tracker.options.max_step_size = 0.01
        @test length(collect(iterator(tracker, s, 1.0, 0.0))) ≥ 101

        F = System([x - a], [x], [a])
        ct = Tracker(
            ParameterHomotopy(F, [1], [2]),
            options = TrackerOptions(max_step_size = 0.015625),
        )
        Xs = Vector{ComplexF64}[]
        for (x, t) in iterator(ct, [1.0], 1.0, 0.0)
            push!(Xs, x)
        end

        @test length(Xs) ≥ length(1:0.015625:2)
    end

    @testset "Change parameters" begin
        @var x a y b
        F = System([x^2 - a, x * y - a + b]; parameters = [a, b])
        s = [1.0, 1.0 + 0im]
        tracker = Tracker(ParameterHomotopy(F, [2.2, 3.2], [2.2, 3.2]))
        start_parameters!(tracker, [1, 0])
        target_parameters!(tracker, [2, 4])
        res = track(tracker, s, 1.0, 0.0)
        @test is_success(res)
    end

    @testset "Straight Line Homotopy" begin
        @var x y
        F = System([x^2 + y^2 - 2.3, 2x + 3y - 4], [x, y])
        G = System((0.2 + 0.4im) .* [x^2 - 1, y - 1], [x, y])
        H = StraightLineHomotopy(G, F)
        S = [[1, 1], [-1, 1]]
        tracker = Tracker(H, options = (automatic_differentiation = 3,))

        @test is_success(track(tracker, S[1], 1, 0))
        @test is_success(track(tracker, S[2], 1, 0))

        @test is_invalid_startvalue(track(tracker, [100, -100], 1, 0))
    end

    @testset "Compiled/InterpretedHomotopy" begin
        @var x y t γ
        H = Homotopy(
            t * γ * [x^2 - 1, y^2 - 1] + (1 - t) .* [2 * x^2 - 4 + 2x * y, y^2 - 5],
            [x, y],
            t,
            [γ],
        )
        S = total_degree_start_solutions([2, 2])
        for make in [CompiledHomotopy, InterpretedHomotopy]
            tracker = Tracker(
                fix_parameters(make(H), [randn(ComplexF64)]);
                options = (automatic_differentiation = 3,),
            )
            @test all(is_success, track.(tracker, S))
        end

        # also test without parameters
        γ = randn(ComplexF64)
        H = Homotopy(
            t * γ * [x^2 - 1, y^2 - 1] + (1 - t) .* [2 * x^2 - 4 + 2x * y, y^2 - 5],
            [x, y],
            t,
        )
        for make in [CompiledHomotopy, InterpretedHomotopy]
            tracker = Tracker(make(H); options = (automatic_differentiation = 1,))
            @test all(is_success, track.(tracker, S))
        end
    end

    @testset "invalid_startvalue_singular_jacobian" begin
        # https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl/issues/454
        p0 = [
            0 0
            1 0
            1 1
            2 1
            3/2 1/2
            3.0 1.0
            4 1
            5 1
            9/2 1/2
            5 0
            6 0
        ]
        E = [
            1 2
            2 3
            2 5
            3 4
            3 5
            4 5
            4 6
            5 9
            6 7
            7 8
            7 9
            8 9
            8 10
            9 10
            10 11
        ]
        pinnedvertices = [1; 6; 11]
        freevertices = [2; 3; 4; 5; 7; 8; 9; 10]

        @var x[1:size(p0)[1], 1:size(p0)[2]] # x[1:11, 1:2] # create all variables
        xvarz_moving_frame = [Variable(:x, i, k) for i in freevertices for k = 1:2]

        # create random, real-valued, linear equation in the moving frame variables
        # created so that it passes through the initial configuration p0
        a = randn(1, length(xvarz_moving_frame)) # random coefficients
        a = [0 1.0 0 1.0 0 1.0 0 1.0 0 1.0 0 1.0 0 1.0 0 1.0]
        b0 = evaluate(
            a * xvarz_moving_frame,
            [Variable(:x, i, k) => p0[i, k] for i in freevertices for k = 1:2]...,
        )[1]
        # the parameters to move linear "slice" later in a parameter homotopy, just move
        # the constant term "b"
        bvarz = [Variable(:b)]

        # the linear equation with parameters
        L = (a*xvarz_moving_frame)[1] .- Variable(:b)

        ε = 1e-3 #1e-10
        p1 = p0 + ε * randn(size(p0))

        b1 = subs(
            a * xvarz_moving_frame,
            [Variable(:x, i, k) => p1[i, k] for i in freevertices for k = 1:2]...,
        )[1]
        b1 = to_number(b1) # Float64(b1) throws an error

        function edge_equation(i, j)
            eqn = sum([(x[i, k] - x[j, k])^2 for k = 1:2])
            eqn += -sum([(p0[i, k] - p0[j, k])^2 for k = 1:2])
        end
        fs = [edge_equation(E[m, 1], E[m, 2]) for m = 1:size(E)[1]]

        # pin the vertices by substitution
        gs = [
            subs(
                fij,
                [Variable(:x, i, k) => p0[i, k] for i in pinnedvertices for k = 1:2]...,
            ) for fij in fs
        ]

        G = System(vcat(gs, L); variables = xvarz_moving_frame, parameters = bvarz)
        startsolutions0 = [p0[i, k] for i in freevertices for k = 1:2]
        tracker = Tracker(ParameterHomotopy(G, [b0], [b1]; compile = false))
        result = track(tracker, startsolutions0, 1, 0)
        @test is_invalid_startvalue(result)
        @test result.return_code == :terminated_invalid_startvalue_singular_jacobian
    end

    include("test_cases/steiner_higher_prec.jl")
    include("test_cases/four_bar.jl")
end
