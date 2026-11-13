# Certify endpoints with HC.jl's interval-arithmetic (Krawczyk) certification
# (Breiding, Rose, Timme 2023): which successful endpoints are certified solutions, how
# many of them are distinct, and how many are certified real.
#
#   julia --threads=auto --project=benchmarks/harness benchmarks/harness/certify.jl <commit>_<label>_<mode>
#   julia --threads=auto --project=benchmarks/harness benchmarks/harness/certify.jl --starts
#
# The first form certifies a run's endpoints from benchmarks/cache/endpoints/<run>.jls
# (tracker `run`) and the cached reference endpoints of the same paths (tracker
# `reference`), at the target parameters p₁ of each scenario. The
# second certifies the start sets of the artifact `benchmark_starts` (starts.jl) at their
# parameters p₀: each should be certified complete and distinct. Appends one
# row per system and tracker to benchmarks/log/certification.tsv.
using HomotopyContinuation, Serialization, Dates, LinearAlgebra
include(joinpath(@__DIR__, "load_system.jl"))
include(joinpath(@__DIR__, "starts.jl"))

const CACHE = joinpath(@__DIR__, "..", "cache")
# `HC_BENCH_CERTIFICATION=<path>` writes the rows elsewhere
const LOG = get(ENV, "HC_BENCH_CERTIFICATION", joinpath(@__DIR__, "..", "log", "certification.tsv"))
mkpath(dirname(LOG))
const COLUMNS = ["date", "run", "system", "tracker", "paths", "success", "certified",
    "distinct_certified", "real_certified", "distinct_real_certified", "seconds"]

# HC.jl's certify throws (LAPACKException) when an endpoint's Jacobian is exactly singular
# in Float64. Such endpoints are not certifiable: drop them and certify the rest together.
function certify_robust(F, xs, p)
    kw = (; show_progress = false, threading = Threads.nthreads() > 1)
    try
        return certify(F, xs, p; kw...)
    catch e
        e isa Union{CompositeException,TaskFailedException,LinearAlgebra.LAPACKException} || rethrow()
        ok = filter(xs) do x
            try
                certify(F, [x], p; show_progress = false, threading = false)
                true
            catch e
                e isa LinearAlgebra.LAPACKException || rethrow()
                false
            end
        end
        return isempty(ok) ? nothing : certify(F, ok, p; kw...)
    end
end

function certify_row(F, xs, p, run, name, tracker, npaths)
    t = @elapsed C = isempty(xs) ? nothing : certify_robust(F, xs, p)
    row = (; date = Dates.format(now(), "yyyy-mm-ddTHH:MM"), run, system = name, tracker,
        paths = npaths, success = length(xs),
        certified = isnothing(C) ? 0 : ncertified(C),
        distinct_certified = isnothing(C) ? 0 : ndistinct_certified(C),
        real_certified = isnothing(C) ? 0 : nreal_certified(C),
        distinct_real_certified = isnothing(C) ? 0 : ndistinct_real_certified(C),
        seconds = round(t, digits = 2))
    @assert collect(string.(keys(row))) == COLUMNS
    row
end

function log_rows(rows)
    new = !isfile(LOG)
    open(LOG, "a") do io
        new && println(io, join(COLUMNS, '\t'))
        for r in rows
            println(io, join(values(r), '\t'))
            println(rpad(r.system, 46), rpad(r.tracker, 10), r.success, " ok, ", r.certified,
                " certified, ", r.distinct_certified, " distinct, ",
                r.distinct_real_certified, " distinct real  (", r.seconds, " s)")
        end
    end
end

system_of(path) = build_system(load_system(path), Variable, System)

if ARGS == ["--starts"]
    dir = starts_dir()
    rows = map(TOML.parsefile(joinpath(dir, "index.toml"))["sets"]) do set
        certify_row(system_of(set["system"]), read_start_set(dir, set),
            pairs_to_complex(set["p"]), "starts", set["name"], "start set", set["nsolutions"])
    end
    log_rows(rows)
else
    run = only(ARGS)
    endpoints = deserialize(joinpath(CACHE, "endpoints", run * ".jls"))
    rows = NamedTuple[]
    for name in sort(collect(keys(endpoints)))
        data = deserialize(joinpath(CACHE, name * ".jls"))
        F = system_of(data.system)
        e = endpoints[name]
        wx = [x for (s, x) in zip(e.status, e.x) if s === :success]
        hx = [p.x for p in data.hc if p.status === :success]
        push!(rows, certify_row(F, wx, data.p₁, run, name, "run", length(e.x)))
        push!(rows, certify_row(F, hx, data.p₁, run, name, "reference", length(data.hc)))
    end
    log_rows(rows)
end
