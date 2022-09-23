import BenchmarkTools
using HomotopyContinuation
using LinearAlgebra
using PrettyTables

include("../test/test_systems.jl")


function benchsystem!(suite, T, f, mode)
    I = Ref(mode(f))

    u = zeros(T, length(f))
    U = zeros(T, length(f), nvariables(f))
    x = randn(T, nvariables(f))
    p = randn(T, nparameters(f))

    v1 = TaylorVector{2}(randn(T, 2, length(f)))
    v2 = TaylorVector{3}(randn(T, 3, length(f)))
    v3 = TaylorVector{4}(randn(T, 4, length(f)))
    tx0 = TaylorVector{1}(randn(T, 1, nvariables(f)))
    tx1 = TaylorVector{2}(randn(T, 2, nvariables(f)))
    tx2 = TaylorVector{3}(randn(T, 3, nvariables(f)))

    tp1 = TaylorVector{2}(randn(T, 2, nparameters(f)))

    suite["evaluate"] = BenchmarkTools.@benchmarkable evaluate!($u, ($I)[], $x, $p)
    suite["evaluate_and_jacobian"] =
        BenchmarkTools.@benchmarkable evaluate_and_jacobian!($u, $U, ($I)[], $x, $p)
    suite["taylor_1"] =
        BenchmarkTools.@benchmarkable taylor!($v1, $(Val(1)), ($I)[], $tx0, $p)
    suite["taylor_2"] =
        BenchmarkTools.@benchmarkable taylor!($v2, $(Val(2)), ($I)[], $tx1, $tp1)
end

suite = BenchmarkGroup()
for mode in [InterpretedSystem, CompiledSystem]
    mode_suite = BenchmarkTools.BenchmarkGroup()
    suite[string(mode)] = mode_suite
    for (name, system) in TEST_SYSTEM_COLLECTION
        mode_suite[name] = BenchmarkGroup()
        benchsystem!(mode_suite[name], ComplexF64, system, mode)
    end
end

function run_benchmark(file_name)

    loadparams!(suite, BenchmarkTools.load("$(@__DIR__)/params.json")[1], :evals, :samples)

    res = BenchmarkTools.run(suite, verbose = true, seconds = 3)

    BenchmarkTools.save(file_name, res)
    return res
end