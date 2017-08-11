"""
    Spherical(homogenization_variable)

Spherical is a projective algorithm, i.e. you have to provide a variable with which the system will be homogenized.
"""
struct Spherical{V<:MP.AbstractVariable} <:AbstractPredictorCorrectorHomConAlgorithm
    homogenization_variable::V
end

"""
Spherical tangent predictor following the formulation of Chen[^1]. Denote by ``\hat{H}`` the lift of the homotopy ``H``
as a projective mapping.
Euler's method yields
```math
\begin{bmatrix}
D_x\hat{H}^(x,t) \\
x^{H}
\end{bmatrix}
\dot{x} =
\begin{bmatrix}
-D_t\hat{H}^(x,t) \\
0
\end{bmatrix}
```
and we can then obtain the spherical projection:
```math
\cos(\norm{\dot{x}}_2 Δt) x + \sin(\norm{\dot{x}}_2 Δt) \frac{\dot{x}}{\norm{\dot{x}}_2}
```

[^1]: Tianran Chen and Tien-Yien Li.
    “Spherical projective path tracking for homotopy continuation methods”.
    Communications in Information and Systems 12(3):195-220 (2012)
"""
function predict(alg::Spherical, H::AbstractHomotopy{T}, J_H::F1, ∂H∂t::F2, x::Vector{T}, t::Float64, Δt::Float64) where {T<:Complex,F1,F2}
    dot_x = \([J_H(x,t); x'], [-∂H∂t(x,t); 0])
    norm_dot_x = norm(dot_x)
    # TODO: conditining? i.e. check singular values
    cos(-norm_dot_x*Δt) .* x .+ sin(-norm_dot_x*Δt) .* dot_x ./ norm_dot_x
end

function correct!(
    u::Vector{T},
    alg::Spherical,
    H::AbstractHomotopy{T},
    J_H::F,
    x::Vector{T},
    t::Float64,
    tol::Float64,
    max_iterations::Int
) where {T,F}
    N = length(x)
    res = zeros(x)
    dx_xt = zeros(T, N, N)
    Δx = zeros(x)
    u .= x

    for k = 1:max_iterations
        res[1:N-1] = -evaluate(H, u, t)
        # println("newton iteration: $(k), res: $(norm(res))")
        if norm(res) < tol
            return true
        end

        dx_xt[1:(N-1), :] .= J_H(u, t)
        dx_xt[N, :] = u'

        Δx .= \(dx_xt, res)
        norm_Δx = norm(Δx)
        u .= cos(norm_Δx) .* u .+ sin(norm_Δx) .* Δx ./ norm_Δx
    end

    norm(res) < tol
end

function hom_var_index(H::AbstractHomotopy, alg::Spherical)
    index = findfirst(vars(H), alg.homogenization_variable)
    Nullable(index, index > 0)
end
prepare_homotopy(H::AbstractHomotopy, alg::Spherical) = homogenize(H, alg.homogenization_variable)

function prepare_start_value(start_value, H::AbstractHomotopy, alg::Spherical{<:MP.AbstractVariable})
    if length(start_value) == nvars(H)
        return start_value
    end
    projective(start_value, get(hom_var_index(H,alg)))
end