export SphericalPredictorCorrector

struct SphericalPredictorCorrector <: AbstractPredictorCorrectorAlgorithm{Val{true}} end

function precondition!(x::Vector, H::AbstractHomotopy, ::SphericalPredictorCorrector)
    # we work on the unit sphere
    normalize!(x)
end

"""
    predict(alg::Affine, H, J_H, Hdt, x, t, Δt)

Spherical tangent predictor following the formulation of Chen[^1]. Denote by ``\\hat{H}`` the lift of the homotopy ``H``
as a projective mapping.
Euler's method yields
```math
\\begin{bmatrix}
D_x\\hat{H}^(x,t) \\\\
x^{H}
\\end{bmatrix}
\\dot{x} =
\\begin{bmatrix}
-D_t\\hat{H}^(x,t) \\\\
0
\\end{bmatrix}
```
and we can then obtain the spherical projection:
```math
\\cos(\\norm{\\dot{x}}_2 Δt) x + \\sin(\\norm{\\dot{x}}_2 Δt) \\frac{\\dot{x}}{\\norm{\\dot{x}}_2}
```

[^1]: Tianran Chen and Tien-Yien Li.
    “Spherical projective path tracking for homotopy continuation methods”.
    Communications in Information and Systems 12(3):195-220 (2012)
"""
function predict!(
    u::Vector{T},
    A::Matrix{T},
    b::Vector{T},
    H::AbstractHomotopy{T},
    J_H!::Function,
    Hdt!::Function,
    x::Vector{T},
    t::Number,
    Δt::Number,
    alg::SphericalPredictorCorrector) where T
    # fill A
    filljacobi!(A, J_H!, x, t)
    # fill b
    v = @view b[1:end-1]
    Hdt!(v, x, t)
    scale!(v, -one(T))
    b[end] = zero(T)

    # this computes A x = b and stores the result x in b
    LU = lufact!(A)
    A_ldiv_B!(LU, b)
    # lets rename it to make the rest clearer
    dot_x = b

    norm_dot_x = norm(dot_x)
    θ = -norm_dot_x*Δt
    u .= cos(θ) .* x .+ (sin(θ) / norm_dot_x) .* dot_x
    u
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
    alg::SphericalPredictorCorrector) where {T}
    u .= x

    k = 0
    while true
        k += 1
        v = @view b[1:end-1]
        Homotopy.evaluate!(v, H, u, t)
        scale!(v, -one(T))
        b[end] = zero(T)


        # println("newton iteration: $(k), res: $(norm(res))")
        if norm(b, Inf) < tol
            return true
        elseif k > max_iterations
            return false
        end

        filljacobi!(A, J_H!, u, t)

        # this computes A x = b and stores the result x in b
        LU = lufact!(A)
        A_ldiv_B!(LU, b)
        # we rename the result for further computations to Δx
        Δx = b
        norm_Δx = norm(Δx)
        u .= cos(norm_Δx) .* u .+ (sin(norm_Δx) / norm_Δx) .* Δx
    end
end

function filljacobi!(A::Matrix{T}, J_H!::Function, u::Vector{T}, t::Number) where T
    U = @view A[1:end-1, :]
    J_H!(U, u, t)
    for j=1:size(A,2)
        A[end, j] = conj(u[j])
    end
end
