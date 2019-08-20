import DoubleFloats: Double64


issquare(A::AbstractMatrix) = isequal(size(A)...)

"""
    row_scaling!(A, D::AbstractVector, tol=maximum(size(A))^2 * eps(real(eltype(A))))

Divide each row of `A` by its ∞-norm. This minimizes the ∞ condition number
for all possible scalings. In order to avoid scaling near zero rows back to
1 we only scale a row if its ∞-norm is larger than `tol`.
"""
function row_scaling!(A::AbstractMatrix{T}, D::AbstractVector, tol=max(size(A)...)^2 * eps(real(T))) where {T}
    D .= 0.0
    d_max = convert(real(T), -Inf)
    @inbounds for j=1:size(A, 2), i=1:size(A, 1)
        D[i] = @fastmath max(D[i], abs(A[i,j]))
        d_max = @fastmath max(d_max, D[i])
    end

    # In the case that the row norms are all less than 1 we have a lower bound
    # the maximal row norm
    # To avoid that we scale a zero matrix, we need a tol for wich we don't scale any more
    bound = clamp(d_max, tol, 1.0)
    D .= max.(D, bound)

    @inbounds for j=1:size(A, 2), i=1:size(A, 1)
        A[i, j] /= D[i]
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
    corank_proposal::Int
    cond::Float64 # Estimate of the Jacobian
    # The relative number of digits lost during the solution of the linear systems
    # in Newton's method. See `solve_with_digits_lost!` in utilities/linear_algebra.jl
    # for how this is computed.
    digits_lost::Union{Nothing,Float64}
end

function Jacobian(A::AbstractMatrix, corank::Int=0)
    m, n = size(A)
    lu = m ≠ n ? nothing : LinearAlgebra.lu!(Random.randn!(complex.(float.(A))))

    work = Vector{ComplexF64}(undef, 1)
    ormqr_work = Vector{ComplexF64}(undef, 1)
    rwork = Vector{Float64}(undef, 2n)
    jpvt = zeros(BlasInt, n)
    B = ComplexF64.(A)
    tau = similar(B, min(m, n))
    geqp3!(B, jpvt, tau, work, rwork)
    qr = LinearAlgebra.QRPivoted(B, tau, jpvt)
    active_factorization = lu === nothing ? QR_FACTORIZATION : LU_FACTORIZATION
    J = copy(B)
    D = ones(m)
    r = similar(J, size(J, 1))
    b = copy(r)
    perm = zeros(Int, length(qr.p))
    cond = 1.0
    digits_lost = nothing
    Jacobian(J, D, lu, qr, active_factorization, b, r, perm, work, rwork, ormqr_work,
        corank, 0, cond, digits_lost)
end

function update_rank!(Jac::Jacobian)
    Jac.corank = Jac.corank_proposal
    Jac
end

function reset!(jacobian::Jacobian)
    jacobian.cond = 1.0
    jacobian.corank = 0
    jacobian.corank_proposal = 0
    jacobian.digits_lost = nothing
    jacobian
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
function updated_jacobian!(Jac::Jacobian{ComplexF64}; update_infos::Bool=false) where {T}
    if !update_infos && issquare(Jac.J) && Jac.corank == 0
        Jac.lu !== nothing || return Jac
        copyto!(Jac.lu.factors, Jac.J)
        lu!(Jac.lu.factors, nothing, Jac.lu.ipiv)
        Jac.active_factorization = LU_FACTORIZATION
    else
        copyto!(Jac.qr.factors, Jac.J)

        m, n = size(Jac.J)
        if m == n
            # only apply row scaling for square matrices since
            # otherwise we change the problem
            row_scaling!(Jac.qr.factors, Jac.D, 1)
        end
        # this computes a pivoted qr factorization, i.e.,
        # qr!(Jac.qr.factors, Val(true)) but without allocating new memory
        geqp3!(Jac.qr.factors, Jac.qr.jpvt, Jac.qr.τ, Jac.qr_work, Jac.qr_rwork)

        ε = max(n,m) * eps() * 10
        # check rank 0
        rnm = min(n,m)
        r₁ = abs(real(Jac.qr.factors[1,1]))
        corank = 0
        if r₁ < ε
            corank = rnm
            Jac.cond = inv(ε)
        else
            for i in 2:rnm
                rᵢ = abs(real(Jac.qr.factors[i,i]))
                if ε * r₁ > rᵢ
                    corank = rnm - i + 1
                    break
                end
            end
            Jac.cond = max(1.0, r₁) / abs(real(Jac.qr.factors[rnm, rnm]))
        end
        Jac.corank_proposal = corank

        Jac.active_factorization = QR_FACTORIZATION
    end
    Jac
end

"""
    solve!([x, ], Jac::Jacobian, b)

Solve `Jac.J * x = b` inplace. Assumes `Jac` contains already the correct factorization.
This stores the result in `x`.
"""
function solve!(x, Jac::Jacobian, b::AbstractVector; update_digits_lost::Bool=false, refinement_step::Bool=false)
    if Jac.active_factorization == LU_FACTORIZATION
        lu = Jac.lu
        if lu !== nothing
            lu_ldiv!(x, lu, b)
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
    if issquare(Jac.J) && Jac.corank == 0 && update_digits_lost && !refinement_step
        norm_x = euclidean_norm(x)
        norm_δx = iterative_refinement_step!(x, Jac, b)
        Jac.digits_lost = log₁₀(norm_δx / eps(norm_x))
    elseif !issquare(Jac.J)
        Jac.digits_lost = nothing
    end
    x
end
solve!(A::Jacobian, b) = solve!(b, A::Jacobian, b)

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
        ldiv_upper!(A.factors, b; singular_exception=false)
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
    δx = solve!(Jac.r, Jac, Jac.r; refinement_step=true)
    norm_δx = euclidean_norm(δx)
    @inbounds for i in eachindex(x)
        x[i] = convert(T, x[i]) - convert(T, δx[i])
    end

    return norm_δx
end
