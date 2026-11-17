# Run with julia --project=. benchmarks/linear_algebra_kernels.jl
# Compare paired-column solves with the previous scalar-column algorithm. Both
# include the same RHS copy; compilation is excluded. Times are minima of seven
# interleaved batches, in ns per solve. These are kernel, not solver, timings.
using HomotopyContinuation, LinearAlgebra, Random, Printf
const HC = HomotopyContinuation

function scalar_upper!(A, b)
    @inbounds for j = size(A, 2):-1:1
        xj = b[j] = (@fastmath A[j, j] \ b[j])
        for i = 1:j-1
            b[i] -= A[i, j] * xj
        end
    end
    b
end

function scalar_lower!(A, b)
    @inbounds for j = 1:size(A, 2)
        xj = b[j]
        for i = j+1:size(A, 2)
            b[i] -= A[i, j] * xj
        end
    end
    b
end

function kernel_time(f, A, rhs, x, repetitions)
    start = time_ns()
    for _ = 1:repetitions
        copyto!(x, rhs)
        f(A, x)
    end
    (time_ns() - start) / repetitions
end

function benchmark_triangular_solves()
    rng = MersenneTwister(0x6271)
    println("storage       n   upper scalar/paired (ns)   lower scalar/paired (ns)")
    for split in (false, true), n in (3, 6, 12, 24, 32, 48, 100)
        dense = randn(rng, ComplexF64, n, n) / n + 2I
        A = split ? HC.StructArrays.StructArray(dense) : dense
        rhs = randn(rng, ComplexF64, n)
        x = similar(rhs)
        functions = (scalar_upper!, HC.ldiv_upper!, scalar_lower!, HC.ldiv_unit_lower!)
        for f in functions
            kernel_time(f, A, rhs, x, 10)
        end
        times = fill(Inf, 4)
        for _ = 1:7
            for (i, f) in enumerate(functions)
                times[i] = min(times[i], kernel_time(f, A, rhs, x, 20_000))
            end
        end
        @printf(
            "%-11s %3d   %9.1f / %9.1f   %9.1f / %9.1f\n",
            split ? "StructArray" : "Matrix",
            n,
            times...
        )
    end
end

benchmark_triangular_solves()
