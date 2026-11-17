# julia --project=. benchmarks/modelkit_extended.jl [--unfused] [--compiled]
# --unfused removes only the compound double-double specializations in this
# process, leaving precision dispatch unchanged for an equal-precision baseline.
# Run the two modes separately on an otherwise idle machine.
using HomotopyContinuation, Random, Printf
using HomotopyContinuation.DoubleDouble: ComplexDF64
const MK = HomotopyContinuation.ModelKit

if "--unfused" in ARGS
    for (f, n) in (
        (MK.op_muladd, 3),
        (MK.op_mulsub, 3),
        (MK.op_submul, 3),
        (MK.op_mulmuladd, 4),
        (MK.op_mulmulsub, 4),
    )
        Base.delete_method(which(f, Tuple{ntuple(_ -> ComplexDF64, n)...}))
    end
end

include(joinpath(@__DIR__, "../test/test_systems.jl"))

function evaluation_time(f)
    f()
    best = Inf
    for _ = 1:7
        start = time_ns()
        for _ = 1:1000
            f()
        end
        best = min(best, (time_ns() - start) / 1000)
    end
    best
end

function benchmark_extended()
    rng = MersenneTwister(0x6271)
    constructor = "--compiled" in ARGS ? CompiledSystem : InterpretedSystem
    println(constructor, ", ", "--unfused" in ARGS ? "unfused" : "fused")
    println("system          evaluate_ns   eval_jac_ns   taylor3_ns")
    for (name, F) in TEST_SYSTEM_COLLECTION
        name in ("fano_quintic", "four_bar", "steiner", "cyclo") || continue
        G = constructor(F)
        x = ComplexDF64.(randn(rng, ComplexF64, nvariables(F)))
        p = ComplexDF64.(randn(rng, ComplexF64, nparameters(F)))
        u = zeros(ComplexDF64, length(F))
        U = zeros(ComplexDF64, length(F), nvariables(F))
        tx = TaylorVector{3}(ComplexDF64.(randn(rng, ComplexF64, 3, nvariables(F))))
        ev = evaluation_time(() -> evaluate!(u, G, x, p))
        jac = evaluation_time(() -> evaluate_and_jacobian!(u, U, G, x, p))
        tay = evaluation_time(() -> taylor!(u, Val(3), G, tx, p))
        @printf("%-14s %11.1f %13.1f %12.1f\n", name, ev, jac, tay)
    end
end

benchmark_extended()
