# Hand-rolled monodromy loops through `solve`: p₀ → p_a → p_b → p₀ with random complex
# p_a, p_b, starting from the full solution set at p₀ (benchmarks/cache). Every leg should
# return all N solutions, distinct and certified; after a loop the set must equal the one
# at p₀ (a bijection), and the loop's permutation of it is reported.
#   julia --project=benchmarks/harness benchmarks/harness/monodromy.jl <cache name> [loops = 10] [seed = 1]
using HomotopyContinuation, Serialization, Printf, Random, LinearAlgebra
include(joinpath(@__DIR__, "load_system.jl"))
const ROOT = joinpath(@__DIR__, "..", "..")
name = ARGS[1]
loops = parse(Int, get(ARGS, 2, "10"))
rng = MersenneTwister(parse(Int, get(ARGS, 3, "1")))
d = deserialize(joinpath(ROOT, "benchmarks/cache", name * ".jls"))
F = fixed(build_system(load_system(d.system), Variable, System); compile = false)
p₀ = ComplexF64.(d.p₀)
S = [ComplexF64.(x) for x in d.start]
N = length(S)
reld(a, b) = maximum(abs, a .- b) / max(1.0, maximum(abs, b))
leg_solve(X, a, b) = solve(F, X; start_parameters = a, target_parameters = b,
    threading = false, show_progress = false)

# index of the solution in Y nearest to x (bucketed by the rounded first coordinates)
function matcher(Y)
    key(x) = (round(Int, 1e3 * real(x[1]) / max(1, abs(x[1]))), round(Int, 1e3 * imag(x[1]) / max(1, abs(x[1]))))
    B = Dict{Tuple{Int,Int},Vector{Int}}()
    for (j, y) in enumerate(Y)
        push!(get!(B, key(y), Int[]), j)
    end
    x -> begin
        best, bj = Inf, 0
        k = key(x)
        for a = -1:1, b = -1:1, j in get(B, (k[1] + a, k[2] + b), Int[])
            e = reld(x, Y[j]); e < best && ((best, bj) = (e, j))
        end
        (bj, best)
    end
end

@printf("%s: N = %d, %d parameters, %d variables, %d loops\n", name, N, length(p₀), length(S[1]), loops)
leg_solve(S[1:2], p₀, p₀ .+ 0.01) # compile
certify(F, S[1:2], p₀; show_progress = false)
X = copy(S)
total_lost = 0
for l = 1:loops
    pa = randn(rng, ComplexF64, length(p₀)) .* (norm(p₀) / sqrt(length(p₀)))
    pb = randn(rng, ComplexF64, length(p₀)) .* (norm(p₀) / sqrt(length(p₀)))
    for (leg, (a, b)) in enumerate(((p₀, pa), (pa, pb), (pb, p₀)))
        alive = findall(x -> all(isfinite, x), X) # a lost path stays lost
        t = @elapsed R = leg_solve(X[alive], a, b)
        sol = solutions(R; only_nonsingular = false)
        C = certify(F, sol, b; show_progress = false)
        nd = ndistinct_certified(C)
        global total_lost += N - nd
        @printf("  loop %d leg %d: success %5d distinct %5d certified distinct %5d singular %d failed %d %6.2f s\n",
            l, leg, count(is_success, path_results(R)), length(sol), nd, nsingular(R),
            count(r -> !is_success(r), path_results(R)), t)
        flush(stdout)
        # keep the start order: X[i] is where start i went (a failed path gets NaN)
        Xn = [fill(NaN + 0im, length(S[1])) for _ = 1:N]
        for r in path_results(R)
            is_success(r) && (Xn[alive[r.path_number]] = solution(r))
        end
        global X = Xn
    end
    # back at p₀: the same set, permuted
    m = matcher(S)
    hit = zeros(Int, N)
    worst = 0.0
    perm = zeros(Int, N)
    for i = 1:N
        all(isfinite, X[i]) || continue
        j, e = m(X[i])
        (j > 0 && e < 1e-8) || continue
        hit[j] += 1; perm[i] = j; worst = max(worst, e)
    end
    @printf("  loop %d back at p₀: matched %d of %d, each original solution hit once: %s, worst rel. distance %.1e, moved by the loop %d\n",
        l, count(>(0), perm), N, all(==(1), hit), worst, count(i -> perm[i] != 0 && perm[i] != i, 1:N))
    flush(stdout)
end
println("total missing certified solutions over all legs: ", total_lost)
