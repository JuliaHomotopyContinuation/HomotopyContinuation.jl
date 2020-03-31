@testset "PathTracker" begin
    @var x y
    f = [2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3, 2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5]
    H, starts = total_degree_homotopy(f, [x, y]; gamma = 1.3im + 0.4)
    S = collect(starts)
    tracker = PathTracker(Tracker(H))
    @test !isempty(sprint(show, tracker.options))

    res = track.(tracker, S)

    @test is_success(res[1])
    @test all(isfinite, solution(res[1]))
    @test is_finite(res[1])
    @test winding_number(res[1]) == nothing
    @test !is_failed(res[1])
    @test accuracy(res[1]) < 1e-12
    @test steps(res[1]) < 10
    @test accepted_steps(res[1]) < 10
    @test rejected_steps(res[1]) == 0
    @test last_path_point(res[1]) == nothing
    @test is_success(res[2])
    @test is_at_infinity(res[3])
    @test is_at_infinity(res[4])
    @test !isempty(sprint(show, res[1]))
    @test valuation(res[3]) ≈ [-1,-1] rtol = 1e-3

    # singular stuff
    @var x
    f = [(x - 10)^2]
    H, starts = total_degree_homotopy(f, [x])
    tracker = PathTracker(Tracker(H))
    res = track.(tracker, starts)

    @test winding_number(res[1]) == 2
    @test last_path_point(res[1]) isa Tuple{Vector{ComplexF64}, Float64}
    @test 0 < last(last_path_point(res[1])) < 0.1
end
