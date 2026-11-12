struct CorrectorCountingHomotopy{H} <: HC.AbstractHomotopy
    homotopy::H
    evaluations::Base.RefValue{Int}
    jacobians::Base.RefValue{Int}
    extended_jacobians::Base.RefValue{Int}
end
CorrectorCountingHomotopy(H) = CorrectorCountingHomotopy(H, Ref(0), Ref(0), Ref(0))
Base.size(H::CorrectorCountingHomotopy) = size(H.homotopy)
HC.ModelKit.taylor!(u, v::Val, H::CorrectorCountingHomotopy, x, t, incremental::Bool) =
    taylor!(u, v, H.homotopy, x, t, incremental)
HC.extended_jacobian_supported(H::CorrectorCountingHomotopy) =
    HC.extended_jacobian_supported(H.homotopy)
function HC.ModelKit.evaluate!(u, H::CorrectorCountingHomotopy, x, t)
    H.evaluations[] += 1
    evaluate!(u, H.homotopy, x, t)
end
function HC.ModelKit.evaluate_and_jacobian!(u, U, H::CorrectorCountingHomotopy, x, t)
    H.jacobians[] += 1
    H.extended_jacobians[] += eltype(x) === HC.DoubleDouble.ComplexDF64
    evaluate_and_jacobian!(u, U, H.homotopy, x, t)
end

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

    @testset "projective tracking $mode - AD: $AD" for mode in
                                                       [InterpretedSystem, CompiledSystem],
        AD = 0:3

        @var x a y b z
        F = System([x^2 - a * z^2, x * y + (b - a) * z^2], [x, y, z], [a, b])
        H = ParameterHomotopy(mode(F), [1, 0], [2, 4])
        tracker = Tracker(
            on_affine_chart(H, (2,)),
            options = TrackerOptions(automatic_differentiation = AD),
        )

        s = [1, 1, 1]
        res = track(tracker, s, 1, 0)
        @test is_success(res)
        @test isa(solution(res), Vector{ComplexF64})
        x₀ = abs(solution(res)[end])
        @test affine_chart(PVector(solution(res))) ≈ [sqrt(2), -sqrt(2)] rtol = 1e-12 / x₀
    end

    @testset "Paths through zero" begin
        @var y p
        cases = (
            (y * (p - y), 1 + im, 5, 0, 0),
            ((y - 50) * (p - y), -1, -2, 50, 50),
            (y - p, 1, -2, 1, -2),
            (y - p^3, 1, -2, 1, -8),
            # At the start, the first three derivatives vanish, but the path
            # is not constant.
            (y - p^4, 0, 1, 0, 1),
        )
        for (f, p_start, p_target, y_start, y_target) in cases
            F = System([f], [y], [p])
            tracker = Tracker(ParameterHomotopy(F, [p_start], [p_target]))
            result = track(tracker, [y_start], 1, 0)
            @test is_success(result)
            @test norm(solution(result) - [y_target]) < 1e-10
        end
    end

    @testset "Parameter solve with a zero solution" begin
        @var y p
        F = System([y * (p - y)], [y], [p])
        S1 = solve(F; target_parameters = [1 + im], show_progress = false)
        S2 = solve(
            F,
            S1;
            start_parameters = [1 + im],
            target_parameters = [5.0],
            show_progress = false,
        )
        @test nsolutions(S2) == 2
        @test all(is_success, S2)
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

    @testset "path info" begin
        @var x a y b
        F = System([x^2 - a, x * y - a + b], [x, y], [a, b])
        tracker = Tracker(ParameterHomotopy(F, [1, 0], [2, 4]))
        info = path_info(tracker, [1, 1], 1, 0)
        @test !isempty(sprint(show, info))
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

    @testset "Newton correction at an unchanged point" begin
        @var x y t
        symbolic = Homotopy([x - 1 - t, 2 * (y - 1)], [x, y], t)
        for mode in (InterpretedHomotopy, CompiledHomotopy, MixedHomotopy),
            extended in (false, true),
            weighted in (false, true),
            δ in (0.0, 2.0^-60)

            H = CorrectorCountingHomotopy(mode(symbolic))
            start = ComplexF64[1, 1]
            out = similar(start)
            nc = HC.NewtonCorrector(0.125, start, 2)
            J = HC.Jacobian(zeros(ComplexF64, 2, 2))
            norm = weighted ? WeightedNorm(ones(2), InfNorm()) : InfNorm()
            r = HC.newton!(
                out,
                nc,
                H,
                start,
                δ,
                J,
                norm;
                μ = eps(),
                ω = 1.0,
                extended_precision = extended,
                first_correction = true,
            )
            @test HC.is_converged(r)
            @test out == start
            @test r.accuracy ≈ δ
            @test r.iters == 2
            @test H.jacobians[] == 1
            @test J.factorizations[] == 1
            @test H.evaluations[] == Int(extended)
            extended && @test r.μ_low ≈ δ
        end
    end

    @testset "Newton correction with changing iterates" begin
        @var x y t
        symbolic = Homotopy([x^2 - 2, y^2 - 3], [x, y], t)
        for extended in (false, true), accurate in (false, true)
            H = CorrectorCountingHomotopy(InterpretedHomotopy(symbolic))
            start = ComplexF64[1.42, 1.74]
            out = similar(start)
            nc = HC.NewtonCorrector(0.125, start, 2)
            J = HC.Jacobian(zeros(ComplexF64, 2, 2))
            r = HC.newton!(
                out,
                nc,
                H,
                start,
                0.0,
                J,
                InfNorm();
                μ = eps(),
                ω = 1.0,
                extended_precision = extended,
                accurate_μ = accurate,
                first_correction = true,
            )
            @test HC.is_converged(r)
            @test out ≈ [sqrt(2), sqrt(3)] atol = 1e-14
            @test r.accuracy <= 1e-14
            @test H.jacobians[] >= 2
            @test all(isfinite, HC.matrix(J))
        end

        # Rejection and singularity must not become successful stagnation exits.
        for H in (
            InterpretedHomotopy(Homotopy([x^2 + 1, y^2 + 1], [x, y], t)),
            InterpretedHomotopy(Homotopy([x + y - 1, x + y - 1], [x, y], t)),
        )
            start = ComplexF64[0.1, 0.1]
            nc = HC.NewtonCorrector(0.125, start, 2)
            r = HC.newton!(
                similar(start),
                nc,
                H,
                start,
                0.0,
                HC.Jacobian(zeros(ComplexF64, 2, 2)),
                InfNorm();
                μ = eps(),
                ω = 1.0,
            )
            @test !HC.is_converged(r)
        end
    end

    @testset "Persistent path precision" begin
        @var x y p q
        F = System(
            [x + y - 1, x + y - 1 + p * ((1 + q) * y - 1)];
            variables = [x, y],
            parameters = [p, q],
        )
        for compile in (false, true, :mixed)
            H = CorrectorCountingHomotopy(
                ParameterHomotopy(F, [2.0^-60, 1], [2.0^-60, 0]; compile),
            )
            T = Tracker(H)
            @test HC.init!(T, [0.5, 0.5])
            @test T.corrector.full_extended_precision[]
            @test T.predictor.used_extended_precision
            low = H.jacobians[] - H.extended_jacobians[]
            high = H.extended_jacobians[]
            @test HC.step!(T)
            @test H.extended_jacobians[] > high
            @test H.jacobians[] - H.extended_jacobians[] == low
            @test T.corrector.full_extended_precision[]
            @test T.state.used_extended_prec
            HC.init!(T, 0.0)
            @test T.corrector.full_extended_precision[]
            W = T.corrector.extended_workspace[]
            P = T.predictor.extended_workspace

            # Supplying initial estimates bypasses Newton initialization, so
            # here the predictor discovers the need for escalation first.
            T2 = Tracker(H)
            @test HC.init!(T2, [0.5, 0.5]; ω = 1.0, μ = eps())
            @test T2.predictor.used_extended_precision
            @test T2.corrector.full_extended_precision[]
            @test T2.corrector.extended_workspace[] === nothing
            low = H.jacobians[] - H.extended_jacobians[]
            @test HC.step!(T2)
            @test H.jacobians[] - H.extended_jacobians[] == low
            @test T2.corrector.extended_workspace[] !== nothing

            HC.parameters!(H.homotopy, [1, 1], [1, 0])
            @test HC.init!(T, [0.5, 0.5])
            @test !T.corrector.full_extended_precision[]
            @test !T.predictor.used_extended_precision
            @test T.corrector.extended_workspace[] === W
            @test T.predictor.extended_workspace === P

            HC.parameters!(H.homotopy, [2.0^-60, 1], [2.0^-60, 0])
            T.options.extended_precision = false
            high = H.extended_jacobians[]
            @test !HC.init!(T, [0.5, 0.5])
            @test !T.corrector.full_extended_precision[]
            @test H.extended_jacobians[] == high
        end
    end

    include("test_cases/steiner_higher_prec.jl")
    include("test_cases/four_bar.jl")
end
