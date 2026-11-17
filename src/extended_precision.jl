# Extended-precision cache.
# Predictor and Newton method decide when to use this.

# Dense double-double storage for one n×n solve A·c = r.
# `A` and `residual` are filled by evaluate_and_jacobian!,
# `factors`/`pivots` hold the LU factorization
# the solution is stored in `correction` .
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

# Solve using the factors from `extended_factorize!`, writing to `correction`.
# Returns `false` if the solution is non-finite or fails a termwise backward-error
# test |A c − r|ᵢ ≤ 1e-25 · (|r|ᵢ + Σⱼ|Aᵢⱼ cⱼ|). The check is termwise.
# The output `false` means "no usable correction".
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

# Helper function that decides whether extended precision for the Jacobian is available for this homotopy.
# By default, this returns false for user-defined AbstractHomotopys, since their evaluate_and_jacobian! need not accept ComplexDF64. 
extended_jacobian_supported(::AbstractHomotopy) = false
extended_jacobian_supported(::Union{CompiledHomotopy,InterpretedHomotopy,MixedHomotopy}) =
    true
extended_jacobian_supported(H::ParameterHomotopy) =
    H.F isa Union{CompiledSystem,InterpretedSystem,MixedSystem}
extended_jacobian_supported(H::FixedParameterHomotopy) =
    extended_jacobian_supported(H.homotopy)
extended_jacobian_supported(H::AffineChartHomotopy) =
    extended_jacobian_supported(H.homotopy)
