# julia --project=. benchmarks/qr_factorization.jl
# Times the QR factorization and the least-squares solve of overdetermined
# MatrixWorkspaces (split real and imaginary parts), as used for the Jacobians of
# overdetermined systems. The factorization time includes copying the matrix into the
# workspace; compilation is excluded. Times are minima of five batches, in ns.
using HomotopyContinuation, LinearAlgebra, Random, Printf
const HC = HomotopyContinuation

function batch_time(f, count)
    start = time_ns()
    for _ = 1:count
        f()
    end
    (time_ns() - start) / count
end

function benchmark_workspace_qr()
    rng = MersenneTwister(42)
    println("  m    n   factorize (ns)   solve (ns)")
    for (m, n) in ((3, 2), (5, 4), (6, 5), (11, 10), (21, 20), (41, 40), (60, 20))
        A = randn(rng, ComplexF64, m, n)
        b = randn(rng, ComplexF64, m)
        x = zeros(ComplexF64, n)
        WS = HC.MatrixWorkspace(A)
        factorize() = (copyto!(WS.A, A); HC.updated!(WS); HC.factorize!(WS))
        solve() = ldiv!(x, WS, b)
        factorize()
        solve()
        qr, sv = Inf, Inf
        for _ = 1:5
            qr = min(qr, batch_time(factorize, 1000))
            sv = min(sv, batch_time(solve, 1000))
        end
        @printf("%3d  %3d   %14.1f   %10.1f\n", m, n, qr, sv)
    end
end

benchmark_workspace_qr()
