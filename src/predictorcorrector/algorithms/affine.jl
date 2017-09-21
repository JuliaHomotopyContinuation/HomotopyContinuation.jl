export AffinePredictorCorrector

struct AffinePredictorCorrector <: AbstractPredictorCorrectorAlgorithm{Val{false}} end

# nothing to do here
precondition!(x::Vector, H::AbstractHomotopy, ::AffinePredictorCorrector) = x

"""
    predict(alg::Affine, H, J_H, ∂H∂t, x, t, Δt)

Tangent predictor. `dx_xt` and `dt_xt` are the evaluations of partial derivatives from our homotopy at x and t, i.e.
``D_xH(x,t)`` and ``D_tH(x,t)``.
It is based on the explicit euler method:
```math
x(t)−∆t·D_xH^{-1}(x,t)·D_tH(x,t)
```

We have ``v = D_xH^{-1}(x,t)·D_tH(x,t)``,
therefore ``D_xH(x,t)·v = D_tH(x,t)``
and thus ``v = \\(D_xH(x,t), D_tH(x,t)``
"""
function predict!(
    u::Vector{T}, A::Matrix{T}, b::Vector{T},
    H::AbstractHomotopy{T}, J_H!::Function, Hdt!::Function,
    x::Vector{T}, t::Number, Δt::Number, alg::AffinePredictorCorrector) where T

    # fill A and b
    J_H!(A, x, t)
    Hdt!(b, x, t)
    # compute A \ b and store the result in b
    LU = lufact!(A)
    A_ldiv_B!(LU, b)

    x .- Δt .* b
end

function correct!(
    u::Vector{T},
    A::Matrix{T},
    b::Vector{T},
    H::AbstractHomotopy,
    J_H!::Function,
    x::Vector{T},
    t::Number,
    tol::Float64,
    max_iterations::Int,
    alg::AffinePredictorCorrector
) where {T}
    u .= x
    k = 0
    while true
        k += 1
        Homotopy.evaluate!(b, H, u, t)
        if norm(b, Inf) < tol
            return true
        elseif k > max_iterations
            return false
        end
        J_H!(A, u, t)
        # compute A \ b and store the result in b
        LU = lufact!(A)
        A_ldiv_B!(LU, b)
        u .= u .- b
    end
end
