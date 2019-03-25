using HomotopyContinuation
using TestSystems
using BenchmarkTools
p1 = TotalDegreeProblem(equations(katsura(7)()))
P = ProjectiveProblem(p1)

sols = embed.(P, totaldegree_solutions(p1.system) |> collect)

tracker = CoreTracker(P.homotopy, first(sols), 1.0, 0.0)

track(tracker, sols, 1.0, 0.0)

Profile.clear_malloc_data()

@btime track($tracker, $sols, 1.0, 0.0)


H, x, t = tracker.cache.homotopy, tracker.state.x, HomotopyContinuation.CoreTrackers.current_t(tracker.state)
y = (randn(9).*0.1 .- 0.1)
z = x + y

pc = tracker.method.predictor_corrector
pcc = tracker.cache.predictor_corrector

@btime HomotopyContinuation.PredictionCorrection.refine!($z,
     $pc, $pcc,
    $H, $t, 1e-7, 20)
