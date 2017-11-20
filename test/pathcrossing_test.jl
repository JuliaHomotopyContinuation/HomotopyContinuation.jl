@testset "Pathcrossing" begin
    H, startvalues = randomhomotopy(StraightLineHomotopy, 4)

    solver = Solver(H; path_precision=1e-3, corrector_maxiters=5, apply_gammatrick=false)

    tracked_paths = map(startvalues) do startvalue
        track!(solver.pathtracker, startvalue, 1.0, 0.1)
        PathtrackerResult(solver.pathtracker, false)
    end
    crossed_indices = HomotopyContinuation.check_crossed_paths(tracked_paths, 1e-8, true)

    before = length(crossed_indices)
    @test before ≥ 0
    HomotopyContinuation.pathcrossing_check!(tracked_paths, solver)
    crossed_indices = HomotopyContinuation.check_crossed_paths(tracked_paths, 1e-8, true)
    @test before == 0 || before > length(crossed_indices)

    res = solve(H, startvalues, parallel_type=:none, path_precision=1e-3, corrector_maxiters=5, apply_gammatrick=false)
end
