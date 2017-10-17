export SphericalPredictorCorrector

# struct SphericalPredictorCorrector <: PathtrackingAlgorithm end

function precondition!(x::Vector, H::AbstractHomotopy, ::SphericalPredictorCorrector)
    # we work on the unit sphere
    normalize!(x)
end

"""
    predict(alg::Affine, H, J_H, Hdt, x, t, Δt)

Spherical tangent predictor.
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
    #scale!(v, -one(T))
    b[end] = zero(T)

    # this computes A x = b and stores the result x in b
    LU = lufact!(A)
    A_ldiv_B!(LU, b)
    @. u = x - Δt * b
    normalize!(u)
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
    k = 0
    while true
        k += 1
        v = @view b[1:end-1]
        Homotopy.evaluate!(v, H, u, t)
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
        u .= u .- b
        normalize(u)
    end
end

function filljacobi!(A::Matrix{T}, J_H!::Function, u::Vector{T}, t::Number) where T
    U = @view A[1:end-1, :]
    J_H!(U, u, t)
    for j=1:size(A,2)
        A[end, j] = conj(u[j])
    end
end
