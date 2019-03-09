"""
    Jacobian{T, <:Factorization}

Data structure holding a Jacobian matrix `J` with it's factorization.
Additional contains a vector `D` of scalings which row equilibarate `J`.

    Jacobian(J)

Setup a Jacobian with Jacobian matrices of the form `J`.
"""
mutable struct Jacobian{T, F<:LinearAlgebra.Factorization}
    J::Matrix{T} # Jacobian matrix
    D::Vector{Float64} # Scaling factors D
    fac::F # Factorization of D * J
    b::Vector{T} # Vector for overdetermined systems
    r::Vector{T} # Vector for iterative refinement
end

function Jacobian(A::AbstractMatrix)
    m, n = size(A)
    if m == n
        fac = LinearAlgebra.lu(A)
    else
        fac = LinearAlgebra.qr(A)
    end
    J = copy(A)
    D = ones(m)
    r = similar(J, size(J, 1))
    b = copy(r)
    Jacobian(J, D, fac, b, r)
end

"""
    update!(jacobian, J)

Update the `jacobian` struct with the new Matrix `J`.
"""
function update!(Jac::Jacobian, J::AbstractMatrix)
    copyto!(Jac.J, J)
    Jac.fac = factorize!(Jac.fac, Jac.J)
    Jac
end

"""
    update_jacobian!(jacobian)

This computes a factorization of the stored Jacobian matrix.
Call this instead of `update!` if `jacobian.J` got updated.
"""
function updated_jacobian!(Jac::Jacobian)
    Jac.fac = factorize!(Jac.fac, Jac.J)
    Jac
end



"""
    solve!([x,] A, b)

Solve ``Ax=b`` inplace. This overwrites `A` and `b`
and stores the result in `x`. If `x` is not provided result is stored in `b`
"""
function solve!(x, A::StridedMatrix, b::StridedVecOrMat)
    m, n = size(A)
    if m == n
        if x !== b
            copyto!(x, b)
        end
        # solve using an LU factorization
        lu_factorization!(A, x, nothing)
        # now forward and backward substitution
        ldiv_unit_lower!(A, x)
        ldiv_upper!(A, x)
    else
        LinearAlgebra.ldiv!(x, LinearAlgebra.qr!(A), b)
    end
    x
end
solve!(A, b) = solve!(b, A, b)

"""
    solve!([x, ] factorization, b)

Solve ``Ax=b`` inplace where already a factorization of `A` is provided.
This stores the result in `b`.
"""
function solve!(x, Jac::Jacobian, b::AbstractVector)
    if x !== b
        if length(b) == length(x)
            rhs = x
        else
            # overdetermined, we need an intermediate structure
            rhs = Jac.b
        end
        copyto!(rhs, b)
    else
        rhs = x
    end

    solve!(x, Jac.fac, rhs)
end
function solve!(x, LU::LinearAlgebra.LU, b::AbstractVector)
    if x !== b
        copyto!(x, b)
    end
     _ipiv!(LU, x)
     ldiv_unit_lower!(LU.factors, x)
     ldiv_upper!(LU.factors, x)
     x
 end
solve!(x, fact::LinearAlgebra.Factorization, b::AbstractVector) = LinearAlgebra.ldiv!(x, fact, b)


"""
    factorization(A::AbstractMatrix)

Compute a factorization of `A`. This depends on the shape of `A`.
"""
function factorization(A::AbstractMatrix)
    m, n = size(A)
    if m == n
        LinearAlgebra.lu(A)
    else
        LinearAlgebra.qr(A)
    end
end

"""
    factorize!(factorization, A)

Compute a factorization of `A` with the same type as `factorization`.
"""
function factorize!(LU::LinearAlgebra.LU, A::AbstractMatrix)
    copyto!(LU.factors, A)
    lu_factorization!(LU.factors, nothing, LU.ipiv)
    LU
end
function factorize!(QR::LinearAlgebra.QRCompactWY, A::AbstractMatrix)
    copyto!(QR.factors, A)
    LinearAlgebra.qr!(QR.factors)
end


# This is an adoption of LinearAlgebra.generic_lufact!
# with 3 changes:
# 1) For choosing the pivot we use abs2 instead of abs
# 2) Instead of using the robust complex division we
#    just use the naive division. This shouldn't
#    lead to problems since numerical diffuculties only
#    arise for very large or small exponents
# 3) We fold lu! and ldiv! into one routine
#    this has the effect that we do not need to allocate
#    the pivot vector anymore and also avoid the allocations
#    coming from the LU wrapper
@inline function lu_factorization!(A::AbstractMatrix{T},
                           b::Union{AbstractVector, Nothing}=nothing,
                           ipiv::Union{Vector{I}, Nothing}=nothing,
                           ::Val{Pivot} = Val(true)) where {T,I,Pivot}
    m, n = size(A)
    minmn = min(m,n)
    # LU Factorization
    @inbounds begin
        for k = 1:minmn
            # find index max
            kp = k
            if Pivot
                amax = zero(real(T))
                for i = k:m
                    if T <: Complex
                        absi = abs2(A[i,k])
                    else
                        absi = abs(A[i,k])
                    end
                    if absi > amax
                        kp = i
                        amax = absi
                    end
                end
            end
            if !(ipiv isa Nothing)
                ipiv[k] = kp
            end
            if !iszero(A[kp,k])
                if k != kp
                    # Interchange
                    for i = 1:n
                        tmp = A[k,i]
                        A[k,i] = A[kp,i]
                        A[kp,i] = tmp
                    end
                    if !(b isa Nothing)
                        b[k], b[kp] = b[kp], b[k]
                    end
                end
                # Scale first column
                Akkinv = @fastmath inv(A[k,k])
                for i = k+1:m
                    A[i,k] *= Akkinv
                end
            end
            # Update the rest
            for j = k+1:n
                for i = k+1:m
                    A[i,j] -= A[i,k]*A[k,j]
                end
            end
        end
    end
    A
end


function _ipiv!(A::LinearAlgebra.LU, b::AbstractVector)
    for i = 1:length(A.ipiv)
        if i != A.ipiv[i]
            _swap_rows!(b, i, A.ipiv[i])
        end
    end
    b
end

function _swap_rows!(B::StridedVector, i::Integer, j::Integer)
    B[i], B[j] = B[j], B[i]
    B
end

@inline function ldiv_upper!(A::AbstractMatrix, b::AbstractVector, x::AbstractVector = b)
    n = size(A, 2)
    for j in n:-1:1
        @inbounds iszero(A[j,j]) && throw(LinearAlgebra.SingularException(j))
        @inbounds xj = x[j] = (@fastmath A[j,j] \ b[j])
        for i in 1:(j-1)
            @inbounds b[i] -= A[i,j] * xj
        end
    end
    b
end
@inline function ldiv_unit_lower!(A::AbstractMatrix, b::AbstractVector, x::AbstractVector = b)
    n = size(A, 2)
    @inbounds for j in 1:n
        xj = x[j] = b[j]
        for i in j+1:n
            b[i] -= A[i,j] * xj
        end
    end
    x
end


"""
    iterative_refinement_step!(x, A, b, fac, T, r)

Apply one step of iterative refinement where the residual is computed with precision `T`.
Returns the euclidean norm of the update.
"""
function iterative_refinement_step!(x, A, b, fac, r, ::Type{T}=eltype(x)) where T
    residual!(r, A, x, b, T)
    δx = solve!(fac, r)
    norm_δx = euclidean_norm(δx)
    @inbounds for i in eachindex(x)
        x[i] = convert(T, x[i]) - convert(T, δx[i])
    end

    return norm_δx
end

"""
    residual!(u, A, x, b, ::Type{T}=eltype(u))

Compute the residual `Ax-b` in precision `T` and store in `u`.
"""
function residual!(u::AbstractVector, A, x, b, ::Type{T}=eltype(u)) where {T}
    @boundscheck size(A, 1) == length(b) && size(A,2) == length(x)
    m, n = size(A)
    @inbounds for i in 1:m
        dot = zero(T)
        for j in 1:n
            dot = muladd(convert(T, A[i,j]), convert(T, x[j]), dot)
        end
        u[i] = dot - convert(T, b[i])
    end
    u
end


"""
    solve_with_digits_lost!(x, Jac::Jac, b)::Float64

Solve `Jac.J * x = b` and apply one step of iterative refinment to get an estimate of the
(relative) lost digits.
Let ``δx`` be the update of the iterative refinement step.
Then, the estimate is computed by ``log₁₀(||δx||₂/ ϵ(||x||₂))`` where ``ϵ`` is the machine precision (`eps` in Julia).
This is a lowerbound of the logarithm of the condition number, i.e., ``log₁₀(κ(J))``.
The estimate is returned.
"""
function solve_with_digits_lost!(x::AbstractVector, Jac::Jacobian, b::AbstractVector)
    solve!(x, Jac.fac, b)
    norm_x = euclidean_norm(x)
    norm_δx = iterative_refinement_step!(x, Jac.J, b, Jac.fac, Jac.r)
    log₁₀(norm_δx / eps(norm_x))
end
