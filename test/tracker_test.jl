@testset "Tracker" begin

    @testset "tracking" begin
        @var x a y b
        F = System([x^2 - a, x * y - a + b], [x, y], [a, b])

        tracker = Tracker(ParameterHomotopy(F, [1, 0], [2, 4]))
        s = [1, 1]
        res = track(tracker, s, 1, 0)
        @test is_success(res)
        @test isa(solution(res), Vector{ComplexF64})
        @test length(solution(res)) == 2
        @test is_success(track(tracker, res, 0, 1))

        code, μ, ω = track!(tracker, s, 1, 0)
        @test is_success(code)
        @test !isnan(μ)
        @test !isnan(ω)
        s0 = copy(tracker.state.x)
        @test is_success(track(tracker, s0, 0, 1))
        @test is_success(track(tracker, s0, 0, 1, μ = μ, ω = ω))

        s = @SVector [1, 1]
        @test is_success(track(tracker, s, 1, 0))
    end

    @testset "iterator" begin
        @var x a y b
        F = System([x^2 - a, x * y - a + b], [x, y], [a, b])
        tracker = Tracker(ParameterHomotopy(F, [1, 0], [2, 4]))
        s = [1, 1]

        # path iterator
        typeof(first(iterator(tracker, s, 1.0, 0.0))) == Tuple{Vector{ComplexF64},Float64}

        tracker.options.max_step_size = 0.01
        @test length(collect(iterator(tracker, s, 1.0, 0.0))) == 101

        F = System([x - a], [x], [a])
        ct = Tracker(ParameterHomotopy(F, [1], [2]), max_step_size = 0.1)
        Xs = Vector{ComplexF64}[]
        for (x, t) in iterator(ct, [1.0], 1.0, 0.0)
            push!(Xs, x)
        end
        @test round.(Int, real.(first.(Xs)) .* 10) == collect(10:20)
    end

    @testset "path info" begin
        @var x a y b
        F = System([x^2 - a, x * y - a + b], [x, y], [a, b])
        tracker = Tracker(ParameterHomotopy(F, [1, 0], [2, 4]))
        info = path_info(tracker, [1,1], 1, 0)
        @test !isempty(sprint(show, info))
    end

    include("test_cases/steiner_higher_prec.jl")
    include("test_cases/four_bar.jl")
end
