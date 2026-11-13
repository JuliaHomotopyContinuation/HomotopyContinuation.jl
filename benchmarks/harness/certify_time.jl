# Single-threaded time of `certify` on a run's endpoints, after a warm-up call
# (the `seconds` column of benchmarks/log/certification.tsv includes compilation for the
# first system of a run and uses all threads).
#
#   julia --project=benchmarks/harness benchmarks/harness/certify_time.jl <commit>_<label>_<mode> [system ...]
using HomotopyContinuation, Serialization, Printf
include(joinpath(@__DIR__, "load_system.jl"))

const CACHE = joinpath(@__DIR__, "..", "cache")
run = ARGS[1]
E = deserialize(joinpath(CACHE, "endpoints", run * ".jls"))
for name in (length(ARGS) > 1 ? ARGS[2:end] : sort(collect(keys(E))))
    data = deserialize(joinpath(CACHE, name * ".jls"))
    F = build_system(load_system(data.system), Variable, System)
    e = E[name]
    xs = [x for (s, x) in zip(e.status, e.x) if s === :success]
    isempty(xs) && continue
    certify(F, xs[1:1], data.p₁; show_progress = false, threading = false)
    t = @elapsed C = certify(F, xs, data.p₁; show_progress = false, threading = false)
    @printf("%-36s %6d endpoints  %d/%d/%d certified/distinct/real  %.1f µs/endpoint\n", name,
        length(xs), ncertified(C), ndistinct_certified(C), nreal_certified(C), 1e6t / length(xs))
end
