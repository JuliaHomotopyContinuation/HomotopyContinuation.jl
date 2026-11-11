# Compare the existing single-column algorithm with paired-column updates.
using HomotopyContinuation, Random, Printf
const HC = HomotopyContinuation

function timing(f, A, B, piv, count)
    start = time_ns()
    for _ = 1:count
        copyto!(B, A)
        f(B, piv)
    end
    (time_ns() - start) / count
end
function run()
    rng = MersenneTwister(42)
    for split in (false, true), n in (2, 3, 4, 16, 24, 32, 48, 100)
        A = randn(rng, ComplexF64, n, n)
        split && (A = HC.StructArrays.StructArray(A))
        B, piv = copy(A), zeros(Int, n)
        timing(HC.lu_rank1!, A, B, piv, 10)
        timing(HC.lu!, A, B, piv, 10)
        old, new = Inf, Inf
        for _ = 1:5
            old = min(old, timing(HC.lu_rank1!, A, B, piv, 1000))
            new = min(new, timing(HC.lu!, A, B, piv, 1000))
        end
        @printf(
            "%-11s n=%3d scalar/optimized: %.2fx\n",
            split ? "StructArray" : "Matrix",
            n,
            old / new
        )
    end
end
run()
