issquare(A::AbstractMatrix) = isequal(size(A)...)

"""
    row_scaling!(A, tol=maximum(size(A))^2 * eps(real(eltype(A))))

Divide each row of `A` by its ∞-norm. This minimizes the ∞ condition number
for all possible scalings. In order to avoid scaling near zero rows back to
1 we only scale a row if its ∞-norm is larger than `tol`.
"""
function row_scaling!(A::AbstractMatrix{T}, tol=max(size(A)...)^2 * eps(real(T))) where {T}
    @inbounds for i=1:size(A, 1)
        rᵢ = abs(A[i, 1])
        for j=2:size(A, 2)
            rᵢ = max(rᵢ, abs(A[i, j]))
        end
        rᵢ > tol || continue
        for j=1:size(A, 2)
            A[i, j] /= rᵢ
        end
    end
    A
end

function row_scaling!(A::AbstractMatrix{T}, D::AbstractVector, tol=max(size(A)...)^2 * eps(real(T))) where {T}
    @inbounds for i=1:size(A, 1)
        rᵢ = abs(A[i, 1])
        for j=2:size(A, 2)
            rᵢ = max(rᵢ, abs(A[i, j]))
        end
        if rᵢ < tol
            D[i] = one(rᵢ)
            continue
        end
        D[i] = rᵢ

        for j=1:size(A, 2)
            A[i, j] /= rᵢ
        end
    end
    A
end

@enum ActiveFactorization begin
    QR_FACTORIZATION
    LU_FACTORIZATION
end

"""
    Jacobian{T, <:Factorization}

Data structure holding a Jacobian matrix `J` with it's factorization.
Additional contains a vector `D` of scalings which row equilibarate `J`.

    Jacobian(J)

Setup a Jacobian with Jacobian matrices of the form `J`.
"""
mutable struct Jacobian{T}
    J::Matrix{T} # Jacobian matrix
    D::Vector{Float64} # Inverse of scaling factors D
    # lu is only there for square systems
    lu::Union{Nothing, LinearAlgebra.LU{T,Matrix{T}}} # LU Factorization of D * J
    qr::LinearAlgebra.QRPivoted{T,Matrix{T}} # Pivoted QR Factorization of D * J
    active_factorization::ActiveFactorization
    b::Vector{T} # Vector for overdetermined systems
    r::Vector{T} # Vector for iterative refinement
    perm::Vector{Int} # work copy of the permutation in the qr
    qr_work::Vector{ComplexF64} # preallocated memory for pivoted qr
    qr_rwork::Vector{Float64} # preallocated memory for pivoted qr
    ormqr_work::Vector{ComplexF64}
    # Numerical Informations
    corank::Int # The numerical corank
    cond::Float64 # Estimate of the Jacobian
    # The relative number of digits lost during the solution of the linear systems
    # in Newton's method. See `solve_with_digits_lost!` in utilities/linear_algebra.jl
    # for how this is computed.
    digits_lost::Union{Nothing,Float64}
end

function Jacobian(A::AbstractMatrix, corank::Int=0)
    m, n = size(A)
    lu = m ≠ n ? nothing : LinearAlgebra.lu(A)

    work = Vector{ComplexF64}(undef, 1)
    ormqr_work = Vector{ComplexF64}(undef, 1)
    rwork = Vector{Float64}(undef, 2n)
    jpvt = zeros(BlasInt, n)
    tau = similar(A, min(m, n))
    qr = LinearAlgebra.QRPivoted(geqp3!(copy(A), jpvt, tau, work, rwork)...)
    active_factorization = lu === nothing ? QR_FACTORIZATION : LU_FACTORIZATION
    J = copy(A)
    D = ones(m)
    r = similar(J, size(J, 1))
    b = copy(r)
    perm = zeros(Int, length(qr.p))
    cond = 1.0
    digits_lost = nothing
    Jacobian(J, D, lu, qr, active_factorization, b, r, perm, work, rwork, ormqr_work,
        corank, cond, digits_lost)
end

"""
    update!(jacobian, J, corank::Int=jacobian.corank)

Update the `jacobian` struct with the new Matrix `J`.
"""
function update!(Jac::Jacobian, J::AbstractMatrix; update_infos::Bool=false)
    copyto!(Jac.J, J)
    updated_jacobian!(Jac; update_infos=update_infos)
    Jac
end

"""
    update_jacobian!(jacobian)

This computes a factorization of the stored Jacobian matrix.
Call this instead of `update!` if `jacobian.J` got updated.
"""
function updated_jacobian!(Jac::Jacobian{T}; update_infos::Bool=false) where {T}
    if !update_infos && issquare(Jac.J) && Jac.corank == 0
        Jac.lu !== nothing || return Jac
        copyto!(Jac.lu.factors, Jac.J)
        lu_factorization!(Jac.lu.factors, nothing, Jac.lu.ipiv)
        Jac.active_factorization = LU_FACTORIZATION
    else
        copyto!(Jac.qr.factors, Jac.J)

        m, n = size(Jac.J)
        if m == n
            # only apply row scaling for square matrices since
            # otherwise we change the problem
            row_scaling!(Jac.qr.factors, Jac.D)
        end
        # this computes a pivoted qr factorization, i.e.,
        # qr!(Jac.qr.factors, Val(true)) but without allocating new memory
        geqp3!(Jac.qr.factors, Jac.qr.jpvt, Jac.qr.τ, Jac.qr_work, Jac.qr_rwork)

        ε = max(n,m) * eps(real(T)) * 100
        # check rank 0
        rnm = min(n,m)
        r₁ = abs(real(Jac.qr.factors[1,1]))
        if r₁ < ε
            Jac.corank = rnm
        else
            for i in 2:rnm
                if ε * r₁ > abs(real(Jac.qr.factors[i,i]))
                    Jac.corank = rnm - i + 1
                    break
                end
            end
        end
        # compute subcondition number
        Jac.cond = r₁ / abs(real(Jac.qr.factors[rnm, rnm]))
        Jac.active_factorization = QR_FACTORIZATION
    end
    Jac
end

"""
    solve!([x, ], Jac::Jacobian, b)

Solve `Jac.J * x = b` inplace. Assumes `Jac` contains already the correct factorization.
This stores the result in `x`.
"""
function solve!(x, Jac::Jacobian, b::AbstractVector; update_digits_lost=false)
    if Jac.active_factorization == LU_FACTORIZATION
        lu = Jac.lu
        if lu !== nothing
            custom_ldiv!(x, lu, b)
        end
    elseif Jac.active_factorization == QR_FACTORIZATION
        if issquare(Jac.J)
            # we applied row scaling for square matrices
            @inbounds for i in eachindex(Jac.b)
                Jac.b[i] = b[i] / Jac.D[i]
            end
            custom_ldiv!(x, Jac.qr, Jac.b, Jac.corank, Jac.perm, Jac.ormqr_work )
        else
            Jac.b .= b
            custom_ldiv!(x, Jac.qr, Jac.b, Jac.corank, Jac.perm, Jac.ormqr_work )
        end
    end
    if issquare(Jac.J) && Jac.corank == 0 && update_digits_lost
        norm_x = euclidean_norm(x)
        norm_δx = iterative_refinement_step!(x, Jac, b)
        Jac.digits_lost = log₁₀(norm_δx / eps(norm_x))
    elseif !issquare(Jac.J)
        Jac.digits_lost = nothing
    end
    x
end
solve!(A::Jacobian, b) = solve!(b, A::Jacobian, b)

function custom_ldiv!(x, LU::LinearAlgebra.LU, b::AbstractVector)
    if x !== b
        copyto!(x, b)
    end
     _ipiv!(LU, x)
     ldiv_unit_lower!(LU.factors, x)
     ldiv_upper!(LU.factors, x)
     x
end

const LA = LinearAlgebra

# This is an adoption of the Julia implementation of ldiv! for QRPivoted
function custom_ldiv!(x, A::LinearAlgebra.QRPivoted{ComplexF64},
                      b::Vector{ComplexF64}, corank::Int,
                      perm::Vector{Int}, ormqr_work::Vector{ComplexF64})
    mA, nA = size(A.factors)
    nr = min(mA,nA)
    nrhs = length(b)
    if nr == 0
        return x
    end
    if corank == nr
        fill!(x, zero(eltype(x)))
        return x
    end
    rank = nr - corank

    if rank == nA == nr
        # The following is equivalent to LA.lmul!(LA.adjoint(A.Q), b) but
        # we can also pass the preallocated work vector
        ormqr!('L', 'C', A.factors, A.τ, b, ormqr_work)
        ldiv_upper!(A.factors, b)
        @inbounds for i in 1:nr
            x[i] = b[i]
            perm[i] = A.p[i]
        end
        Base.invpermute!!(x, perm)
    else
        C, τ = LA.LAPACK.tzrzf!(A.factors[1:rank,:])
        LA.ldiv!(LA.UpperTriangular(C[1:rank,1:rank]),view(LA.lmul!(LA.adjoint(A.Q), b), 1:rank))
        for i in 1:rank
            x[i] = b[i]
        end
        for i in rank+1:nA
            x[i] = zero(eltype(x))
        end
        LA.LAPACK.ormrz!('L', eltype(x)<:Complex ? 'C' : 'T', C, τ, reshape(x, length(x), 1))
        invpermute!(x, A.p)
    end
    return x
end



"""
    iterative_refinement_step!(x, A, b, fac, T, r)

Apply one step of iterative refinement where the residual is computed with precision `T`.
Returns the euclidean norm of the update.
"""
function iterative_refinement_step!(x, Jac::Jacobian, b, ::Type{T}=eltype(x)) where T
    residual!(Jac.r, Jac.J, x, b, T)
    δx = solve!(Jac, Jac.r)
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


# """
#     solve_with_digits_lost!(x, Jac::Jac, b)::Float64
#
# Solve `Jac.J * x = b` and apply one step of iterative refinment to get an estimate of the
# (relative) lost digits.
# Let ``δx`` be the update of the iterative refinement step.
# Then, the estimate is computed by ``log₁₀(||δx||₂/ ϵ(||x||₂))`` where ``ϵ`` is the machine precision (`eps` in Julia).
# This is a lower bound of the logarithm of the condition number, i.e., ``log₁₀(κ(J))``.
# The estimate is returned.
# """
# function solve_with_digits_lost!(x::AbstractVector, Jac::Jacobian, b::AbstractVector)
#     solve!(x, Jac, b)
#     norm_x = euclidean_norm(x)
#     norm_δx = iterative_refinement_step!(x, Jac, b)
#     log₁₀(norm_δx / eps(norm_x))
# end


###########################
# CUSTOM LU Factorization #
###########################

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


_ipiv!(A::LinearAlgebra.LU, b::AbstractVector) = apply_ipiv!(b, A.ipiv)

function apply_ipiv!(b::AbstractVector, ipiv)
    @inbounds for i = 1:length(ipiv)
        if i != ipiv[i]
            _swap_rows!(b, i, ipiv[i])
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

import LinearAlgebra: BlasInt

"""
    geqp3!(A::AbstractMatrix{ComplexF64},
        jpvt::AbstractVector{BlasInt},
        tau::AbstractVector{ComplexF64},
        work::Vector{ComplexF64},
        rwork=Vector{Float64}(undef, 2*size(A,2)))

Version of `LAPACK.geqp3!` which also accepts the work buffers.
"""
function geqp3!(A::AbstractMatrix{ComplexF64},
		jpvt::AbstractVector{BlasInt},
		tau::AbstractVector{ComplexF64},
		work::Vector{ComplexF64},
		rwork=Vector{Float64}(undef, 2*size(A,2)))
	m,n = size(A)
	lda = stride(A,2)
	if lda == 0
		return return A, tau, jpvt
	end # Early exit
	lwork = BlasInt(-1)
	info = Ref{BlasInt}()
	for i = 1:2
		ccall((LinearAlgebra.LAPACK.@blasfunc(zgeqp3_), LinearAlgebra.LAPACK.liblapack), Cvoid,
			  (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
			   Ptr{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
			   Ptr{Float64}, Ptr{BlasInt}),
			  m, n, A, lda,
			  jpvt, tau, work, lwork,
			  rwork, info)
		LinearAlgebra.LAPACK.chklapackerror(info[])
		if i == 1
			lwork = BlasInt(real(work[1]))
			resize!(work, lwork)
		end
	end
	return A, tau, jpvt
end

"""
    ormqr!(side, trans, A, tau, C, work)

Version of `LAPACK.ormqr!` which also accepts the work buffer.
Computes `Q * C` (`trans = N`), `transpose(Q) * C` (`trans = T`), `adjoint(Q) * C`
(`trans = C`) for `side = L` or the equivalent right-sided multiplication
for `side = R` using `Q` from a `QR` factorization of `A` computed using
`geqrf!`. `C` is overwritten.
"""
function ormqr!(side::AbstractChar, trans::AbstractChar, A::Matrix{ComplexF64},
                tau::Vector{ComplexF64}, C::Vector{ComplexF64},
                work::Vector{ComplexF64})
    m,n = (size(C, 1), 1)
    mA  = size(A, 1)
    k   = length(tau)
    if side == 'L' && m != mA
        throw(DimensionMismatch("for a left-sided multiplication, the first dimension of C, $m, must equal the second dimension of A, $mA"))
    end
    if side == 'R' && n != mA
        throw(DimensionMismatch("for a right-sided multiplication, the second dimension of C, $m, must equal the second dimension of A, $mA"))
    end
    if side == 'L' && k > m
        throw(DimensionMismatch("invalid number of reflectors: k = $k should be <= m = $m"))
    end
    if side == 'R' && k > n
        throw(DimensionMismatch("invalid number of reflectors: k = $k should be <= n = $n"))
    end
    lwork = BlasInt(-1)
    info  = Ref{BlasInt}()
    for i = 1:2  # first call returns lwork as work[1]
        ccall((LA.LAPACK.@blasfunc(zunmqr_), LA.LAPACK.liblapack), Cvoid,
              (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
               Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
               Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
               Ptr{BlasInt}),
              side, trans, m, n,
              k, A, max(1,stride(A,2)), tau,
              C, max(1, stride(C,2)), work, lwork,
              info)
        LA.LAPACK.chklapackerror(info[])
        if i == 1
            lwork = BlasInt(real(work[1]))
            resize!(work, lwork)
        end
    end
    C
end
