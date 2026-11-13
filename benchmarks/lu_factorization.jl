# julia --project=. benchmarks/lu_factorization.jl
# Times the LU factorization and solve of square MatrixWorkspaces (split real and
# imaginary parts), as used for the Jacobians during path tracking. The factorization
# time includes copying the matrix into the workspace; compilation is excluded. Times
# are minima of five batches, in ns.
using HomotopyContinuation, LinearAlgebra, Random, Printf
const HC = HomotopyContinuation

function batch_time(f, count)
    start = time_ns()
    for _ = 1:count
        f()
    end
    (time_ns() - start) / count
end

function benchmark_workspace_lu()
    rng = MersenneTwister(42)
    println("  n   factorize (ns)   solve (ns)")
    for n in (2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 100)
        A = randn(rng, ComplexF64, n, n)
        b = randn(rng, ComplexF64, n)
        x = similar(b)
        WS = HC.MatrixWorkspace(A)
        factorize() = (copyto!(WS.A, A); HC.updated!(WS); HC.factorize!(WS))
        solve() = ldiv!(x, WS, b)
        factorize()
        solve()
        lu, sv = Inf, Inf
        for _ = 1:5
            lu = min(lu, batch_time(factorize, 1000))
            sv = min(sv, batch_time(solve, 1000))
        end
        @printf("%3d   %14.1f   %10.1f\n", n, lu, sv)
    end
end

benchmark_workspace_lu()
