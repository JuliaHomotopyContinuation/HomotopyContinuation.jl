# Shared extended-precision storage, checked linear solves, and capabilities.
# Predictor and Newton algorithms decide when to use this machinery.

# Shared, lazily allocated storage for extended predictor and Newton solves.
struct ExtendedLinearWorkspace
    A::Matrix{ComplexDF64}
    factors::Matrix{ComplexDF64}
    residual::Vector{ComplexDF64}
    correction::Vector{ComplexDF64}
    pivots::Vector{Int}
end
ExtendedLinearWorkspace(n) = ExtendedLinearWorkspace(
    zeros(ComplexDF64, n, n),
    zeros(ComplexDF64, n, n),
    zeros(ComplexDF64, n),
    zeros(ComplexDF64, n),
    zeros(Int, n),
)

function extended_factorize!(W::ExtendedLinearWorkspace)
    copyto!(W.factors, W.A)
    lu_rank2!(W.factors, W.pivots)
    W
end

function extended_ldiv!(W::ExtendedLinearWorkspace)
    n = length(W.correction)
    copyto!(W.correction, W.residual)
    @inbounds for k = 1:n
        p = W.pivots[k]
        W.correction[k], W.correction[p] = W.correction[p], W.correction[k]
    end
    @inbounds for k = 1:n, i = k+1:n
        W.correction[i] -= W.factors[i, k] * W.correction[k]
    end
    @inbounds for k = n:-1:1
        W.correction[k] /= W.factors[k, k]
        for i = 1:k-1
            W.correction[i] -= W.factors[i, k] * W.correction[k]
        end
    end
    all(isfinite, W.correction) || return false
    # Termwise backward error permits cancellation in the right-hand side.
    @inbounds for i = 1:n
        residual = -W.residual[i]
        scale = abs(W.residual[i])
        for j = 1:n
            term = W.A[i, j] * W.correction[j]
            residual += term
            scale += abs(term)
        end
        abs(residual) <= 1e-25 * scale || return false
    end
    true
end

struct ExtendedPredictorWorkspace
    linear::ExtendedLinearWorkspace
    tx::TaylorVector{4,ComplexDF64}
end
ExtendedPredictorWorkspace(n) =
    ExtendedPredictorWorkspace(ExtendedLinearWorkspace(n), TaylorVector{4}(ComplexDF64, n))

# Escalate the Jacobian and factorization only when residual-only refinement
# stagnates. User-defined systems keep their existing evaluation contract.
extended_jacobian_supported(::AbstractHomotopy) = false
extended_jacobian_supported(::Union{CompiledHomotopy,InterpretedHomotopy,MixedHomotopy}) =
    true
extended_jacobian_supported(H::ParameterHomotopy) =
    H.F isa Union{CompiledSystem,InterpretedSystem,MixedSystem}
extended_jacobian_supported(H::FixedParameterHomotopy) =
    extended_jacobian_supported(H.homotopy)
extended_jacobian_supported(H::AffineChartHomotopy) =
    extended_jacobian_supported(H.homotopy)
