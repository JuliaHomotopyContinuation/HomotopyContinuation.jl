export FixedStepSize

struct FixedStepSize <: AbstractPredictorCorrectorAlgorithm{Val{false}} end

# nothing to do here
precondition!(x::Vector, H::AbstractHomotopy, ::FixedStepSize) = x

"""
    predict(alg::Affine, H, J_H, ∂H∂t, x, t, Δt)

Uses the current iterate x as input for the Newton method corrector

"""
function predict!(
    u::Vector{T}, A::Matrix{T}, b::Vector{T},
    H::AbstractHomotopy{T}, J_H!::Function, Hdt!::Function,
    x::Vector{T}, t::Number, Δt::Number, alg::FixedStepSize) where T

    x
    
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
    alg::FixedStepSize
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
        # this computes A x = b and stores the result x in b
        LU = lufact!(A)
        A_ldiv_B!(LU, b)
        u .= u .- b
    end
end
