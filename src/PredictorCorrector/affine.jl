struct Affine <: AbstractPredictorCorrectorAlgorithm{false} end

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
function predict(alg::Affine, H::AbstractHomotopy, J_H, ∂H∂t, x, t, Δt)
    x .- Δt .* \(J_H(x,t), ∂H∂t(x,t))
end

function correct!(
    u::Vector{T},
    alg::Affine,
    H::AbstractHomotopy,
    J_H,
    x::Vector{T},
    t,
    tol::Float64,
    max_iterations::Int
) where {T}
    N = length(x)
    res = zeros(x)
    dx_xt = zeros(T, N, N)
    u .= x
    for _ in 1:max_iterations
        res .= evaluate(H, u, t)
        if norm(res) < tol
            return true
        end
        dx_xt .= J_H(u, t)
        u .= u .- \(dx_xt, res)
    end

    return norm(res) < tol
end
