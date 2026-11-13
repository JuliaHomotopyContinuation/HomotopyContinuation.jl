# Export the start data of the scenario caches (benchmarks/cache/*.jls) as the artifact
# `benchmark_starts` (format in starts.jl) and bind it in Artifacts.toml.
#
#   julia --project=benchmarks/harness benchmarks/harness/export_starts.jl <url> <tarball>
#
# Writes the artifact's gzipped tarball to <tarball>; upload it so that it downloads from
# <url> (e.g. a GitHub release asset), then commit Artifacts.toml. Run this only when the
# start data changes (new systems or scenarios); give the new tarball a new <url>.
import Pkg
using Pkg.Artifacts, Serialization, TOML
include(joinpath(@__DIR__, "load_system.jl"))
include(joinpath(@__DIR__, "starts.jl"))

const CACHE = joinpath(@__DIR__, "..", "cache")
length(ARGS) == 2 || error("usage: export_starts.jl <url> <tarball>")
url, tarball = ARGS

caches = Dict(replace(f, ".jls" => "") => deserialize(joinpath(CACHE, f))
              for f in readdir(CACHE) if endswith(f, ".jls"))

# Scenarios with the same system, p₀ and start solutions share a start set, named after its
# first scenario. A set whose p₀ is the target of another scenario of the same system (the
# chained scenarios of reference.jl) is that scenario's reference endpoints.
groups = Dict{Any,Vector{String}}()
for name in sort(collect(keys(caches)))
    d = caches[name]
    push!(get!(groups, (d.system, d.p₀, d.start), String[]), name)
end
function set_name_and_method(names)
    d = caches[first(names)]
    for (name, e) in caches
        e.system == d.system && e.p₁ == d.p₀ && name ∉ names &&
            return name * "_endpoints", "reference endpoints of scenario $name (Tracker)"
    end
    solset = joinpath(CACHE, "solsets", load_system(d.system).name * ".jls")
    isfile(solset) && deserialize(solset).p == d.p₀ && return first(names), deserialize(solset).method
    first(names), "solve"
end

sets = Dict{String,Any}[]
scenarios = Dict{String,Any}[]
hash = create_artifact() do dir
    for names in values(groups)
        d = caches[first(names)]
        name, method = set_name_and_method(names)
        A = reduce(hcat, d.start)::Matrix{ComplexF64}
        open(io -> write(io, htol.(reinterpret(Float64, A))), joinpath(dir, name * ".bin"), "w")
        push!(sets, Dict("name" => name, "system" => system, "nvariables" => size(A, 1),
            "nsolutions" => size(A, 2), "p" => complex_to_pairs(p₀), "method" => method,
            "hc_version" => string(d.hc_version), "seconds" => d.t_start_solve))
        for n in names
            push!(scenarios, Dict("name" => n, "set" => name,
                "target" => complex_to_pairs(caches[n].p₁)))
        end
    end
    sort!(sets; by = s -> s["name"])
    sort!(scenarios; by = s -> s["name"])
    open(joinpath(dir, "index.toml"), "w") do io
        TOML.print(io, Dict("format" => 1, "sets" => sets, "scenarios" => scenarios);
            sorted = true)
    end
    # the artifact must reproduce the caches exactly
    S = scenario_starts(dir)
    keys(S) == keys(caches) || error("scenarios differ")
    for (n, d) in caches
        s = S[n]
        (s.system, s.p₀, s.p₁, s.start) == (d.system, d.p₀, d.p₁, d.start) || error("$n differs")
    end
end
sha256 = archive_artifact(hash, tarball)
bind_artifact!(joinpath(@__DIR__, "Artifacts.toml"), "benchmark_starts", hash;
    download_info = [(url, sha256)], lazy = true, force = true)
println(length(sets), " start sets, ", length(scenarios), " scenarios; git-tree-sha1 ",
    bytes2hex(hash.bytes), ", tarball ", tarball, " (", round(filesize(tarball) / 2^20, digits = 1),
    " MB, sha256 ", sha256, ")")
