# Run `solve` (path tracking with the endgame) on every cached scenario: the start
# solutions at p₀ to the target parameters p₁, and write one row per scenario.
#
#   julia --project=benchmarks/harness benchmarks/harness/solve.jl <label> [interpreted|compiled] [scaled] [option=value ...]
#
# Without options this runs `solve`'s defaults (single thread). `option=value` arguments
# set `TrackerOptions` fields as in run.jl. `HC_BENCH_ONLY=<regex>` restricts the
# scenarios. Rows go to benchmarks/log/solve_runs.tsv (or `HC_BENCH_RUNS`), the endpoints
# to benchmarks/cache/endpoints/<commit>_solve-<label>_<mode>[_scaled].jls in the format of
# run.jl, so certify.jl and compare.jl work on them unchanged.
using HomotopyContinuation, Serialization, Printf, Dates
include(joinpath(@__DIR__, "load_system.jl"))
include(joinpath(@__DIR__, "scaled.jl"))

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CACHE = joinpath(ROOT, "benchmarks", "cache")
const RUNS = get(ENV, "HC_BENCH_RUNS", joinpath(ROOT, "benchmarks", "log", "solve_runs.tsv"))
const COLUMNS = ["date", "commit", "label", "mode", "system", "paths", "success",
    "nonsingular", "singular", "distinct", "at_infinity", "excess", "failed", "seconds"]

const USAGE = "usage: solve.jl <label> [interpreted|compiled] [scaled] [option=value ...]"
pos = filter(a -> !occursin('=', a), ARGS)
label = isempty(pos) ? error(USAGE) : pos[1]
mode = length(pos) ≥ 2 ? pos[2] : "interpreted"
mode in ("interpreted", "compiled") || error(USAGE)
scaled = length(pos) ≥ 3 && pos[3] == "scaled"
tracker_options = map(filter(a -> occursin('=', a), ARGS)) do a
    k, v = split(a, '='; limit = 2)
    F = fieldtype(TrackerOptions, Symbol(k))
    Symbol(k) => (F <: Number ? parse(F, v) : Symbol(v))
end
only = Regex(get(ENV, "HC_BENCH_ONLY", ""))

sha = readchomp(`git -C $ROOT rev-parse --short HEAD`)
dirty = !isempty(readchomp(`git -C $ROOT status --porcelain -- src benchmarks/harness`))
commit = dirty ? sha * "-dirty" : sha
date = Dates.format(now(), "yyyy-mm-ddTHH:MM")
mkpath(dirname(RUNS))
new = !isfile(RUNS)
io = open(RUNS, "a")
new && println(io, join(COLUMNS, '\t'))
endpoints = Dict{String,Any}()
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
    starts = [s ./ d for s in data.start]
    kw = (; start_parameters = data.p₀, target_parameters = data.p₁, threading = false,
        show_progress = false, tracker_options = TrackerOptions(; tracker_options...))
    solve(G, starts[1:1]; kw...) # compile
    t = @elapsed R = solve(G, starts; kw...)
    # endpoints in start order and original coordinates
    status = fill(:missing, length(starts))
    xs = [fill(complex(NaN), length(d)) for _ in starts]
    for r in path_results(R)
        status[r.path_number] = r.return_code
        xs[r.path_number] = d .* solution(r)
    end
    endpoints[name] = (; status, x = xs)
    row = (; date, commit, label, mode, system = scaled ? name * "~scaled" : name,
        paths = length(starts), success = count(==(:success), status),
        nonsingular = nnonsingular(R), singular = nsingular(R),
        distinct = nsolutions(R; only_nonsingular = false), at_infinity = nat_infinity(R),
        excess = nexcess_solutions(R), failed = nfailed(R), seconds = round(t, sigdigits = 4))
    @assert collect(string.(keys(row))) == COLUMNS
    println(io, join(values(row), '\t'))
    flush(io)
    @printf("%-40s %6d success %6d nonsingular %6d singular %5d distinct %6d at infinity %5d failed %5d  %.2f s\n",
        row.system, row.paths, row.success, row.nonsingular, row.singular, row.distinct,
        row.at_infinity, row.failed, t)
end
close(io)
mkpath(joinpath(CACHE, "endpoints"))
serialize(joinpath(CACHE, "endpoints",
    "$(commit)_solve-$(label)_$(mode)$(scaled ? "_scaled" : "").jls"), endpoints)
