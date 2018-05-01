using HomotopyContinuation
using TestSystems
using BenchmarkTools
p1 = Problems.TotalDegreeProblem(equations(katsura7()))
P = Problems.ProjectiveStartTargetProblem(p1)

sols = Problems.embed.(P, Utilities.totaldegree_solutions(p1.system) |> collect)

tracker = PathTracking.PathTracker(P.homotopy, first(sols), 1.0, 0.0)

PathTracking.track(tracker, sols, 1.0, 0.0)

Profile.clear_malloc_data()

@btime PathTracking.track($tracker, $sols, 1.0, 0.0)


H, x, t = tracker.cache.homotopy, tracker.state.x, HomotopyContinuation.PathTrackers.current_t(tracker.state)
y = (randn(9).*0.1 .- 0.1)
z = x + y

pc = tracker.method.predictor_corrector
pcc = tracker.cache.predictor_corrector

@btime HomotopyContinuation.PredictionCorrection.refine!($z,
     $pc, $pcc,
    $H, $t, 1e-7, 20)
