using HomotopyContinuation
using Compat.Test
import DynamicPolynomials
const PolyImpl = DynamicPolynomials
using TestSystems
using BenchmarkTools
using Bertini

PolyImpl.@polyvar x

katsuras = [katsura5(), katsura6(), katsura7(), katsura8(), katsura9(), katsura10()]
cyclics = [katsura5(), katsura6(), katsura7(), katsura8(), katsura9(), katsura10()]

F = equations(cyclic7())

@time R = solve(F);
count(r -> r.returncode == :success, R)


@time R = solve(F, options=PathTracking.Options(corrector_maxiters=2));
mean(r -> r.iterations, R)
mytimings = map(katsuras) do K
    F = equations(K)
    solve(F)

    mean(1:10) do i
        t = @elapsed(solve(F, options=PathTracking.Options(corrector_maxiters=2)))
        println(t)
        t
    end

end


bertini_timings = map(katsuras) do K
    F = equations(K)
    bertini(F)
    mean(_ -> first(bertini(F, ENDGAMENUM=2)), 1:10)
end

bertini(equations(katsura9()), MAXNEWTONITS=3)

BLAS.set_num_threads(1)
@time solve(equations(katsura10()))


@time R = solve(F)
count(r -> r.returncode != :at_infinity, R)


P1 = Problems.TotalDegreeProblem(F)
P = Problems.ProjectiveStartTargetProblem(P1)
sols = Utilities.totaldegree_solutions(F) |> collect
start_sols = Problems.embed.(P, sols)
s = start_sols[1]
tracker = PathTracking.PathTracker(P.homotopy, s, rand(), rand(), options=PathTracking.Options(tol=1e-7, refinement_tol=1e-13))
endgamer = Endgame.Endgamer(Endgame.Cauchy(), tracker)
for k = 1:120
    r = PathTracking.track(tracker, start_sols[k], 1.0, 0.1)
    er = Endgame.play(endgamer, copy(r.x.data), 0.1)
    if er.returncode != :success
        @show k
        break
    end
end

r = PathTracking.track(tracker, start_sols[3], 1.0, 0.1)
er = Endgame.play(endgamer, copy(r.x.data), 0.1)


endgamer = Endgame.Endgamer(Endgame.Cauchy(), tracker)


Utilities.infinity_norm(er.x, 1)


@btime Endgame.play!($endgamer, $r.x.data, $0.1)


y = copy(tracker.state.x)
w = y[2:end] / y[1]

[g(PolyImpl.variables(G)=> w) for g in G]

P.homotopy(tracker.state.x, 0.0)





solve(G)









solve(F)


@time solve(F)
@time solve(F)
@time solve(F)
@time solve(F)
@time solve(F)

@time R = solve(G)
count(r -> r.residual < 1e-4, R)
mean(r -> r.iterations, R)


A = rand(3, 3)
b = rand(3)
P1 = Problems.TotalDegreeProblem(F)
P = Problems.ProjectiveStartTargetProblem(P1)

sols = Utilities.totaldegree_solutions(F) |> collect

start_sols = Problems.embed.(P, sols)

s = start_sols[1]


const PathTrackers = HomotopyContinuation.PathTrackers
tracker = PathTracking.PathTracker(P.homotopy, s, 1.0, 0.0)

PathTracking.track(tracker, s, 1.0, 0.0)
