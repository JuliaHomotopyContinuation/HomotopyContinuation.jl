# Full solution sets for systems where HC.jl's `solve` is too slow or unreliable, computed
# with HC.jl monodromy. Writes benchmarks/cache/solsets/<name>.jls = (; p, solutions, ...);
# benchmarks/harness/reference.jl uses these as start data (p₀ = p) when present.
#
#   julia --threads=auto --project=benchmarks/harness benchmarks/harness/solution_sets.jl [name ...]
using HomotopyContinuation, Random, Serialization, TOML
include(joinpath(@__DIR__, "load_system.jl"))

const OUT = joinpath(@__DIR__, "..", "cache", "solsets")
mkpath(OUT)

# path => expected number of solutions; four-bar starts from HC.jl's test start pair
const SETS = [
    "hcjl/tritangents.toml" => 720,
    "hcjl/fano_quintic_lines.toml" => 2875,
    "tier_a/alt_nine_point_four_bar.toml" => 8652,
    "tier_a/steiner_conics.toml" => 3264,
]

function start_pair(path, F)
    path == "tier_a/alt_nine_point_four_bar.toml" || return nothing
    d = TOML.parsefile(joinpath(SYSROOT, "hcjl/four_bar_tracker_case.toml"))
    p = parse.(ComplexF64, d["start_parameters"])
    x = parse.(ComplexF64, d["start_solution"])
    r = newton(F, x, p)
    is_success(r) || error("four-bar start solution did not refine")
    (solution(r), p)
end

function solution_set(path, expected)
    sys = load_system(path)
    F = build_system(sys, Variable, System)
    sp = start_pair(path, F)
    t = @elapsed R = isnothing(sp) ?
        monodromy_solve(F; target_solutions_count = expected, max_loops_no_progress = 10,
            show_progress = false) :
        monodromy_solve(F, [sp[1]], sp[2]; target_solutions_count = expected,
            max_loops_no_progress = 10, show_progress = false)
    S = solutions(R)
    data = (; system = path, p = parameters(R), solutions = S, expected,
        method = "HC.jl $(pkgversion(HomotopyContinuation)) monodromy_solve", seconds = t, threads = Threads.nthreads())
    serialize(joinpath(OUT, sys.name * ".jls"), data)
    println(rpad(sys.name, 30), length(S), " / $expected solutions in ", round(t, digits = 1), " s")
    flush(stdout)
end

for (path, n) in SETS
    name = replace(basename(path), ".toml" => "")
    (isempty(ARGS) || name in ARGS) && solution_set(path, n)
end
