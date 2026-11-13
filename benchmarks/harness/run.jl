# Track every cached scenario's start solutions p₀ → p₁ with the checkout's `Tracker` and
# compare with the reference run stored in the cache.
#
#   julia --project=benchmarks/harness benchmarks/harness/run.jl <label> [interpreted|compiled] [scaled] [option=value ...]
#
# Trailing `option=value` arguments set `TrackerOptions` fields, e.g. `parameters=fast`
# or `extended_precision=false`. `HC_BENCH_ONLY=<regex>` restricts the scenarios.
#
# `scaled` tracks every system after a fixed random rescaling G(y) = c .* F(D y), with
# row factors cᵢ and column factors D = diag(dⱼ) drawn from 10^U(-3, 3) (seeded), from
# y₀ = D⁻¹ x₀. The scaling wraps the unchanged program (`Scaled` in scaled.jl): inputs
# are multiplied by d, outputs by c. Endpoints are mapped back (x = D y) before they are
# compared with the reference and certified, and system names get the suffix `~scaled`.
# A scaling-invariant tracker gives the same steps and failures with and without it.
#
# Needs benchmarks/cache/*.jls from reference.jl.
#
# Outputs (see README.md):
#   benchmarks/results/<commit>_<label>_<mode>.tsv   one row per path
#   benchmarks/log/runs.tsv                          one row per system, appended
#   benchmarks/cache/endpoints/<commit>_<label>_<mode>.jls   endpoints, for certification
# `<commit>` is the short HEAD hash, suffixed `-dirty` if the worktree has changes to
# src/ or benchmarks/harness/ (such runs are for exploration, not for the log).
using HomotopyContinuation, Serialization, Statistics, Printf, Dates, Random
include(joinpath(@__DIR__, "load_system.jl"))
include(joinpath(@__DIR__, "scaled.jl"))

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CACHE = joinpath(ROOT, "benchmarks", "cache")
const RESULTS = joinpath(ROOT, "benchmarks", "results")
# `HC_BENCH_RUNS=<path>` writes the per-system rows elsewhere (concurrent sweeps, merged afterwards)
const RUNS = get(ENV, "HC_BENCH_RUNS", joinpath(ROOT, "benchmarks", "log", "runs.tsv"))

const RUN_COLUMNS = [
    "date", "commit", "label", "mode", "julia", "cpu", "threads", "system", "paths",
    # robustness
    "success", "step_size_too_small", "max_steps", "accuracy_limit", "ill_conditioned",
    "other_failures", "duplicates", "mismatch_ref", "ref_success",
    # cost
    "accepted", "rejected", "extended_precision", "total_s", "median_us", "p90_us",
    "ref_accepted", "ref_median_us",
    "results_file",
]

const USAGE = "usage: run.jl <label> [interpreted|compiled] [scaled] [option=value ...]"
pos = filter(a -> !occursin('=', a), ARGS)
label = length(pos) ≥ 1 ? pos[1] : error(USAGE)
mode = length(pos) ≥ 2 ? pos[2] : "interpreted"
mode in ("interpreted", "compiled") || error(USAGE)
scaled = length(pos) ≥ 3 && pos[3] == "scaled"
tracker_options = map(filter(a -> occursin('=', a), ARGS)) do a
    k, v = split(a, '='; limit = 2)
    F = fieldtype(TrackerOptions, Symbol(k))
    Symbol(k) => (F <: Number ? parse(F, v) : Symbol(v))
end
only = Regex(get(ENV, "HC_BENCH_ONLY", ""))

function commit_id()
    sha = readchomp(`git -C $ROOT rev-parse --short HEAD`)
    dirty = !isempty(readchomp(`git -C $ROOT status --porcelain -- src benchmarks/harness`))
    dirty ? sha * "-dirty" : sha
end
cpu() = Sys.cpu_info()[1].model

reldist(a, b) = maximum(abs, a .- b) / max(1.0, maximum(abs, b))

# number of successful endpoints that coincide with an earlier successful endpoint
function duplicates(xs)
    d = 0
    for i in eachindex(xs), j = 1:i-1
        reldist(xs[i], xs[j]) < 1e-6 && (d += 1; break)
    end
    d
end

commit = commit_id()
mkpath(RESULTS)
mkpath(dirname(RUNS))
results_file = "$(commit)_$(label)_$(mode)$(scaled ? "_scaled" : "").tsv"
endpoints = Dict{String,Any}() # name => (; status, x), certified by certify.jl
new_runs = !isfile(RUNS)
runs = open(RUNS, "a")
new_runs && println(runs, join(RUN_COLUMNS, '\t'))
date = Dates.format(now(), "yyyy-mm-ddTHH:MM")

open(joinpath(RESULTS, results_file), "w") do io
    println(io, join(["system", "path", "status", "accepted", "rejected",
        "extended_precision", "seconds", "t", "accuracy", "ref_status", "ref_accepted",
        "ref_seconds", "matches_ref"], '\t'))
    @printf("%-34s %6s %8s %6s %6s %9s %9s %10s %10s\n", "system", "paths", "success",
        "dups", "≠ref", "steps", "ref steps", "µs/path", "ref µs/path")
    for f in sort(filter(f -> endswith(f, ".jls") && occursin(only, f), readdir(CACHE)))
        data = deserialize(joinpath(CACHE, f))
        sys = load_system(data.system)
        name = hasproperty(data, :name) ? data.name : sys.name
        G = fixed(build_system(sys, Variable, System); compile = mode == "compiled" ? :all : false)
        d = ones(length(sys.vars))
        if scaled
            c, d = harness_scaling(sys)
            G = Scaled(G, c, d)
        end
        T = Tracker(ParameterHomotopy(G, data.p₀, data.p₁);
            options = TrackerOptions(; tracker_options...))
        track(T, data.start[1] ./ d, 1.0, 0.0) # compile
        results = map(data.start) do s
            t = @elapsed r = track(T, s ./ d, 1.0, 0.0)
            # back to the original coordinates
            (; status = r.return_code, x = d .* solution(r), t = real(r.t),
                r.accuracy, accepted = r.accepted_steps, rejected = r.rejected_steps,
                extended = r.extended_precision_used, seconds = t)
        end
        endpoints[name] = (; status = [r.status for r in results], x = [r.x for r in results])
        ok(r) = r.status === :success
        mism = 0
        for (i, r) in enumerate(results)
            ref = data.hc[i]
            both = ok(r) && ref.status === :success
            m = both && reldist(r.x, ref.x) < 1e-6
            mism += both && !m
            println(io, join(Any[name, i, r.status, r.accepted, r.rejected, r.extended,
                r.seconds, r.t, r.accuracy, ref.status, ref.accepted, ref.seconds, m], '\t'))
        end
        ts = [r.seconds for r in results]
        nstat(s) = count(r -> r.status === s, results)
        named = [:success, :terminated_step_size_too_small, :terminated_max_steps,
            :terminated_accuracy_limit, :terminated_ill_conditioned]
        row = (;
            date, commit, label, mode, julia = string(VERSION), cpu = cpu(),
            threads = Threads.nthreads(), system = scaled ? name * "~scaled" : name,
            paths = length(results), success = nstat(:success),
            step_size_too_small = nstat(:terminated_step_size_too_small),
            max_steps = nstat(:terminated_max_steps),
            accuracy_limit = nstat(:terminated_accuracy_limit),
            ill_conditioned = nstat(:terminated_ill_conditioned),
            other_failures = count(r -> r.status ∉ named, results),
            duplicates = duplicates([r.x for r in results if ok(r)]),
            mismatch_ref = mism, ref_success = count(p -> p.status === :success, data.hc),
            accepted = sum(r -> r.accepted, results), rejected = sum(r -> r.rejected, results),
            extended_precision = count(r -> r.extended, results),
            total_s = round(sum(ts), sigdigits = 4),
            median_us = round(1e6 * median(ts), digits = 2),
            p90_us = round(1e6 * quantile(ts, 0.9), digits = 2),
            ref_accepted = sum(p -> p.accepted, data.hc),
            ref_median_us = round(1e6 * median(p -> p.seconds, data.hc), digits = 2),
            results_file,
        )
        @assert collect(string.(keys(row))) == RUN_COLUMNS
        println(runs, join(values(row), '\t'))
        flush(runs)
        @printf("%-34s %6d %8d %6d %6d %9d %9d %10.1f %10.1f\n", row.system, row.paths,
            row.success, row.duplicates, mism, row.accepted, row.ref_accepted,
            row.median_us, row.ref_median_us)
    end
end
close(runs)
mkpath(joinpath(CACHE, "endpoints"))
serialize(joinpath(CACHE, "endpoints", replace(results_file, ".tsv" => ".jls")), endpoints)
println("\nlogged to $RUNS, per-path results in benchmarks/results/$results_file")
