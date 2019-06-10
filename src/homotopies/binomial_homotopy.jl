export BinomialHomotopy

"""
    BinomialHomotopy

Construct the homotopy ``x^A - t γ + (1-t) b``
where each column in A denotes a monomial.
"""
struct BinomialHomotopy <: AbstractHomotopy
    A::Matrix{Int}
    b::Vector{ComplexF64}
    γ::Vector{Float64}

end

function BinomialHomotopy(A::Matrix{Int}, b::Vector{ComplexF64})
    γ = sign.(real.(b))
    BinomialHomotopy(A, b, γ)
end

cache(H::BinomialHomotopy, x, t) = HomotopyNullCache()
Base.size(H::BinomialHomotopy) = size(H.A)


function evaluate!(u, H::BinomialHomotopy, x, t, c::HomotopyNullCache)
    @unpack A, b, γ = H
    n = length(b)
    @inbounds for i in 1:n
        xAᵢ = one(x[1])
        for j in 1:n
            if A[j, i] != 0
                xAᵢ *= x[j]^A[j, i]
            end
        end
        u[i] = xAᵢ - (t * γ[i] + (1.0 - t) * b[i])
    end
    u
end

function dt!(u, H::BinomialHomotopy, x, t, c::HomotopyNullCache)
    @unpack b, γ = H
    n = length(b)
    @inbounds for i in 1:n
        u[i] = b[i] - γ[i]
    end
    u
end

function jacobian!(U, H::BinomialHomotopy, x, t, c::HomotopyNullCache)
    @unpack A, b, γ = H
    n = length(b)
    @inbounds for i in 1:n
        xAᵢ = one(x[1])
        for j in 1:n
            if A[j, i] != 0
                xAᵢ *= x[j]^A[j, i]
            end
        end
        for j in 1:n
            if A[j,i] == 0
                U[i,j] = 0
            else
                U[i,j] = A[j,i] * (@fastmath xAᵢ / x[j])
            end
        end
    end
    U
end

function evaluate_and_jacobian!(u, U, H::BinomialHomotopy, x, t, c::HomotopyNullCache)
    @unpack A, b, γ = H
    n = length(b)
    @inbounds for i in 1:n
        xAᵢ = one(x[1])
        for j in 1:n
            if A[j, i] != 0
                xAᵢ *= x[j]^A[j, i]
            end
        end
        u[i] = xAᵢ - (t * γ[i] + (1.0 - t) * b[i])

        for j in 1:n
            if A[j,i] == 0
                U[i,j] = 0
            else
                U[i,j] = A[j,i] * (@fastmath xAᵢ / x[j])
            end
        end
    end
    nothing
end


function evaluate(H::BinomialHomotopy, x, t, c::HomotopyNullCache)
    T = promote_type(eltype(x), ComplexF64)
    evaluate!(zeros(T, size(H.A, 1)), H, x, t, c)
end
(H::BinomialHomotopy)(x, t, c=cache(H, x, t)) = evaluate(H, x, t, c)

function jacobian(H::BinomialHomotopy, x, t, c::HomotopyNullCache)
    T = promote_type(eltype(x), ComplexF64)
    jacobian!(zeros(T, size(H.A)), H, x, t, c)
end
function dt(H::BinomialHomotopy, x, t, c::HomotopyNullCache)
    dt!(zeros(ComplexF64, size(H.A, 1)), H, x, t, c)
end
