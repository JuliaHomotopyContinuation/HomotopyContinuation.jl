# Scenario caches: start solutions and a reference tracking run of every scenario.
#
#   julia --project=benchmarks/harness -e 'using Pkg; Pkg.instantiate()'   # once
#   julia --project=benchmarks/harness benchmarks/harness/reference.jl [system.toml | scenario ...]
#
# For each scenario: the start data (p₀, p₁, all solutions at p₀) from the artifact
# `benchmark_starts` (starts.jl); then every start solution tracked p₀ → p₁ with the
# checkout's `Tracker` (compile = :all, single thread). Writes benchmarks/cache/<name>.jls
# with the start data and per-path results (status, endpoint, accepted/rejected steps,
# seconds), the reference that run.jl compares against. Existing caches are not rebuilt,
# so the reference stays fixed (`hc_version` records which HC.jl built it).
#
# A scenario that is not in the artifact yet (a new system or scenario) gets new start
# data: seeded random complex parameters p₀, p₁ and all solutions at p₀ from `solve`, or
# from benchmarks/cache/solsets/<name>.jls if present (then p₀ is the parameter of that
# set, see solution_sets.jl). export_starts.jl then puts it into a new artifact.
using HomotopyContinuation, Random, Serialization
include(joinpath(@__DIR__, "load_system.jl"))
include(joinpath(@__DIR__, "starts.jl"))

const CACHE = joinpath(@__DIR__, "..", "cache")
mkpath(CACHE)
const STARTS = scenario_starts()
const HC_VERSION = string(pkgversion(HomotopyContinuation), " @ ",
    readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`))

# Extra scenarios: name => (system, fixed target parameters). The start data is the same
# as for the generic scenario of that system; only p₁ differs.
# steiner_conics_real: generic -> real, the five real conics of HC.jl
# test/certification_test.jl ("3264"), which have 3264 real tangent conics.
const STEINER_REAL = Float64[
    10124547 // 662488724, 8554609 // 755781377, 5860508 // 2798943247,
    -251402893 // 1016797750, -25443962 // 277938473, 1,
    520811 // 1788018449, 2183697 // 542440933, 9030222 // 652429049,
    -12680955 // 370629407, -24872323 // 105706890, 1,
    6537193 // 241535591, -7424602 // 363844915, 6264373 // 1630169777,
    13097677 // 39806827, -29825861 // 240478169, 1,
    13173269 // 2284890206, 4510030 // 483147459, 2224435 // 588965799,
    33318719 // 219393000, 92891037 // 755709662, 1,
    8275097 // 452566634, -19174153 // 408565940, 5184916 // 172253855,
    -23713234 // 87670601, 28246737 // 81404569, 1,
]
const SCENARIOS = Dict{String, Tuple{String, Any}}(
    "steiner_conics_real" => ("tier_a/steiner_conics.toml", STEINER_REAL),
)
# <name>_instance: generic → the stored instance of the systems from HC.jl issues
# (target_parameters, or coefficients 1 for coefficient families)
for path in HARNESS_SYSTEMS
    occursin("_issue", path) || continue
    SCENARIOS[replace(basename(path), ".toml" => "") * "_instance"] = (path, instance_parameters(path))
end
# <name>_g1.._g3: generic -> generic, three more random complex targets per hard system
# (an Integer target is the seed of p₁ = randn(MersenneTwister(seed), ComplexF64, m)).
const HARD = ["hcjl/tritangents.toml", "hcjl/fano_quintic_lines.toml",
    "tier_a/alt_nine_point_four_bar.toml", "tier_a/steiner_conics.toml"]
for path in HARD, k in 1:3
    SCENARIOS[replace(basename(path), ".toml" => "") * "_g$k"] = (path, k)
end
# steiner_conics_real_g<k>: generic -> real again, but starting at the p₁ of
# steiner_conics_g<k>, from HC.jl's endpoints there (all paths must succeed, distinct).
const CHAINED = Dict("steiner_conics_real_g$k" => ("steiner_conics_g$k", STEINER_REAL) for k in 1:3)

function chained_start(from)
    d = deserialize(joinpath(CACHE, from * ".jls"))
    all(p -> p.status === :success, d.hc) || error("$from: not all HC.jl paths succeeded")
    S = [p.x for p in d.hc]
    n = length(S)
    length(unique(s -> round.(s, sigdigits = 8), S)) == n || error("$from: endpoints not distinct")
    d.p₁, S
end

# Coefficient families whose stored coefficients span many orders of magnitude: random
# parameters are divided by the magnitude of the stored coefficient of their term, so the
# random members are well scaled and only the stored instance (c = 1) is badly scaled.
const NORMALIZED = ["rc_circuit_badly_scaled_issue414"]
function coefficient_scales(sys)
    env = Dict{Symbol,Any}(v => 1.0 for v in sys.vars)
    [abs(_ev(term.args[3], env)) for e in sys.eqs for term in e.args[2:end]]
end

function reference(path; name = nothing, target = nothing, from = nothing)
    sys = load_system(path)
    name = something(name, sys.name)
    F = build_system(sys, Variable, System)
    rng = MersenneTwister(0x6e8f)
    p₀ = randn(rng, ComplexF64, length(sys.params))
    p₁ = randn(rng, ComplexF64, length(sys.params))
    if sys.name in NORMALIZED
        σ = coefficient_scales(sys)
        p₀ ./= σ
        p₁ ./= σ
    end
    if target isa Integer
        p₁ = randn(MersenneTwister(target), ComplexF64, length(sys.params))
    elseif !isnothing(target)
        p₁ = ComplexF64.(target)
    end
    solset = joinpath(CACHE, "solsets", sys.name * ".jls")
    if haskey(STARTS, name)
        st = STARTS[name]
        st.system == path || error("$name: the artifact has system $(st.system), not $path")
        p₀, p₁, S, t_solve = st.p₀, st.p₁, st.start, st.seconds
    elseif !isnothing(from)
        (p₀, S), t_solve = chained_start(from), 0.0
    elseif isfile(solset) # full set from solution_sets.jl (monodromy)
        set = deserialize(solset)
        p₀, S, t_solve = set.p, set.solutions, set.seconds
    else
        t_solve = @elapsed res = solve(F; target_parameters = p₀, threading = false, show_progress = false)
        S = solutions(res; only_nonsingular = true)
    end
    T = Tracker(ParameterHomotopy(F, p₀, p₁; compile = :all))
    track(T, S[1], 1.0, 0.0) # compile
    paths = map(S) do s
        t = @elapsed r = track(T, s, 1.0, 0.0)
        (; status = Symbol(r.return_code), x = copy(solution(r)),
            accepted = r.accepted_steps, rejected = r.rejected_steps, seconds = t)
    end
    data = (; name, system = path, vars = sys.vars, params = sys.params, p₀, p₁,
        start = S, hc = paths, hc_version = HC_VERSION, t_start_solve = t_solve)
    serialize(joinpath(CACHE, name * ".jls"), data)
    ok = count(p -> p.status === :success, paths)
    println(rpad(name, 34), length(S), " paths, HC.jl success $ok, ",
        round(sum(p -> p.seconds, paths) * 1e3, digits = 1), " ms tracking")
end

# args: system paths or scenario names; without args, all missing caches are built: the
# artifact's scenarios, then the ones defined here but not in the artifact yet (chained
# scenarios last, after the caches they start from)
defined = [HARNESS_SYSTEMS; sort!(collect(keys(SCENARIOS))); sort!(collect(keys(CHAINED)))]
todo = isempty(ARGS) ? [sort!(collect(keys(STARTS))); defined] : ARGS
for arg in todo
    name = haskey(SCENARIOS, arg) || haskey(CHAINED, arg) || haskey(STARTS, arg) ? arg :
        replace(basename(arg), ".toml" => "")
    isempty(ARGS) && isfile(joinpath(CACHE, name * ".jls")) && continue
    if haskey(STARTS, name)
        reference(STARTS[name].system; name)
    elseif haskey(SCENARIOS, arg)
        path, target = SCENARIOS[arg]
        reference(path; name = arg, target)
    elseif haskey(CHAINED, arg)
        from, target = CHAINED[arg]
        reference(SCENARIOS[from][1]; name = arg, target, from)
    else
        reference(arg)
    end
end
