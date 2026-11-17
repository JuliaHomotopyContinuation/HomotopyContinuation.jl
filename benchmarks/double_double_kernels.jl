# julia --project=. benchmarks/double_double_kernels.jl
# Compare ModelKit's fused operations with ordinary complex double-double
# arithmetic. Inputs vary and every result is stored to prevent constant folding.
using HomotopyContinuation, Random, Printf
using HomotopyContinuation.DoubleDouble: DoubleF64, ComplexDF64
const MK = HomotopyContinuation.ModelKit

function arithmetic_batch!(out, f, x, y, z, w)
    @inbounds for i in eachindex(out)
        out[i] = f(x[i], y[i], z[i], w[i])
    end
end

function arithmetic_time(out, f, x, y, z, w)
    arithmetic_batch!(out, f, x, y, z, w)
    best = Inf
    for _ = 1:7
        start = time_ns()
        for _ = 1:1000
            arithmetic_batch!(out, f, x, y, z, w)
        end
        best = min(best, (time_ns() - start) / (1000length(out)))
    end
    best
end

function benchmark_double_double()
    rng = MersenneTwister(0xf053)
    inputs = ntuple(4) do _
        [
            ComplexDF64(
                DoubleF64(randn(rng), randn(rng) * 1e-17),
                DoubleF64(randn(rng), randn(rng) * 1e-17),
            ) for _ = 1:1024
        ]
    end
    out = similar(first(inputs))
    for (name, baseline, fused) in (
        ("multiply", (x, y, z, w) -> x * y, (x, y, z, w) -> MK.op_mul(x, y)),
        ("multiply-add", (x, y, z, w) -> x * y + z, (x, y, z, w) -> MK.op_muladd(x, y, z)),
        (
            "two products",
            (x, y, z, w) -> x * y + z * w,
            (x, y, z, w) -> MK.op_mulmuladd(x, y, z, w),
        ),
    )
        before = arithmetic_time(out, baseline, inputs...)
        after = arithmetic_time(out, fused, inputs...)
        @printf(
            "%-14s generic %.2f ns, fused %.2f ns (%.2fx)\n",
            name,
            before,
            after,
            before / after
        )
    end
end

benchmark_double_double()
