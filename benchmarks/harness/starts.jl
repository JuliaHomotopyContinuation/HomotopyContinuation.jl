# The scenarios' start data, from the artifact `benchmark_starts` (Artifacts.toml; written
# by export_starts.jl, downloaded on first use). The artifact holds
#
#   index.toml   format, and per start set and scenario:
#     [[sets]]       name, system (path relative to benchmarks/systems), nvariables,
#                    nsolutions, p (the parameters of the set, as [re, im] pairs), method,
#                    hc_version, seconds (time to compute the set)
#     [[scenarios]]  name, set, target (p₁, as [re, im] pairs)
#   <set>.bin    the set's solutions: an nvariables × nsolutions ComplexF64 matrix,
#                column-major, little-endian (one solution per column)
#
# A scenario tracks the solutions of its set from p₀ = the set's p to p₁ = its target.
# Changing the start data means a new artifact (export_starts.jl), never editing one.
import Pkg, TOML

starts_dir() = Pkg.Artifacts.ensure_artifact_installed("benchmark_starts",
    joinpath(@__DIR__, "Artifacts.toml"))

pairs_to_complex(a) = [complex(Float64(z[1]), Float64(z[2])) for z in a]
complex_to_pairs(v) = [[real(z), imag(z)] for z in v]

function read_start_set(dir, set)
    A = Matrix{ComplexF64}(undef, set["nvariables"], set["nsolutions"])
    path = joinpath(dir, set["name"] * ".bin")
    filesize(path) == sizeof(A) || error("$path: expected $(sizeof(A)) bytes")
    open(io -> read!(io, A), path)
    F = reinterpret(Float64, A)
    F .= ltoh.(F)
    [A[:, j] for j in axes(A, 2)]
end

"""
    scenario_starts(dir = starts_dir()) -> Dict(name => (; system, p₀, p₁, start, method, hc_version, seconds))
"""
function scenario_starts(dir = starts_dir())
    index = TOML.parsefile(joinpath(dir, "index.toml"))
    index["format"] == 1 || error("unsupported start data format $(index["format"])")
    sets = Dict(s["name"] => s for s in index["sets"])
    starts = Dict{String,Vector{Vector{ComplexF64}}}()
    Dict(map(index["scenarios"]) do sc
        set = sets[sc["set"]]
        start = get!(() -> read_start_set(dir, set), starts, set["name"])
        sc["name"] => (; system = set["system"], p₀ = pairs_to_complex(set["p"]),
            p₁ = pairs_to_complex(sc["target"]), start, method = set["method"],
            hc_version = set["hc_version"], seconds = Float64(set["seconds"]))
    end)
end
