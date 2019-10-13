import LinearAlgebra: BlasInt
const LA = LinearAlgebra

@enum Factorization begin
    QR_FACT
    LU_FACT
end
struct MatrixWorkspace{T} <: AbstractMatrix{Complex{T}}
    A::Matrix{Complex{T}} # Matrix
    d::Vector{T} # Inverse of scaling factors
    # factorizations
    fact::Base.RefValue{Factorization}
    factorized::Base.RefValue{Bool}
    lu::Union{Nothing,LA.LU{Complex{T},Matrix{Complex{T}}}} # LU Factorization of D * J
    qr::LA.QRPivoted{Complex{T},Matrix{Complex{T}}} # Pivoted QR Factorization of D * J
    # qr worksapce
    qr_perm::Vector{Int} # work copy of the permutation in the qr
    qr_work::Vector{Complex{T}} # preallocated memory for pivoted qr
    qr_rwork::Vector{T} # preallocated memory for pivoted qr
    qr_ormqr_work::Vector{Complex{T}}
    # overdetermined
    qr_rhs::Vector{Complex{T}}
    # iterative refinement
    ir_r::Vector{Complex{T}}
    ir_δx::Vector{Complex{T}}
    residual_abs::Vector{T}
    # cond estimators
    lu_cond_work::Vector{Complex{T}}
    lu_cond_rwork::Vector{T}
    qr_cond_work::Vector{ComplexF64}
    qr_cond_rwork::Vector{T}
    # inf_norm_est
    inf_norm_est_work::Vector{Complex{T}}
    inf_norm_est_rwork::Vector{T}
    row_scaling::Vector{T}
end

MatrixWorkspace(Â::AbstractMatrix) =
    MatrixWorkspace(Â, promote_type(real(eltype(Â)), Float64))

function MatrixWorkspace(Â::AbstractMatrix, ::Type{T}) where {T}
    A = Matrix{Complex{T}}(Â)
    m, n = size(A)
    d = ones(m)
    fact = m == n ? LU_FACT : QR_FACT
    # lu
    lu = m == n ? LA.lu(A; check = false) : nothing
    # qr
    qr = LA.qr(A, Val(true))
    qr_perm = zeros(Int, length(qr.p))
    qr_work = Vector{Complex{T}}(undef, 1)
    qr_ormqr_work = Vector{Complex{T}}(undef, 1)
    qr_rwork = Vector{T}(undef, 2n)

    ir_r = Vector{Complex{T}}(undef, m)
    ir_δx = Vector{Complex{T}}(undef, n)
    residual_abs = Vector{T}(undef, m)
    qr_rhs = copy(ir_r)
    lu_cond_work = Vector{Complex{T}}(undef, 2n)
    lu_cond_rwork = Vector{T}(undef, 2n)
    qr_cond_work = Vector{Complex{T}}(undef, 2n)
    qr_cond_rwork = Vector{T}(undef, n)
    inf_norm_est_work = Vector{Complex{T}}(undef, n)
    inf_norm_est_rwork = Vector{T}(undef, n)
    row_scaling = ones(m)
    MatrixWorkspace(
        A,
        d,
        Ref(fact),
        Ref(true),
        lu,
        qr,
        qr_perm,
        qr_work,
        qr_rwork,
        qr_ormqr_work,
        qr_rhs,
        ir_r,
        ir_δx,
        residual_abs,
        lu_cond_work,
        lu_cond_rwork,
        qr_cond_work,
        qr_cond_rwork,
        inf_norm_est_work,
        inf_norm_est_rwork,
        row_scaling,
    )
end

Base.size(MW::MatrixWorkspace) = size(MW.A)

import Base: @propagate_inbounds
@propagate_inbounds Base.getindex(MW::MatrixWorkspace, i::Integer) = getindex(MW.A, i)
@propagate_inbounds Base.getindex(MW::MatrixWorkspace, i::Integer, j::Integer) =
    getindex(MW.A, i, j)
@propagate_inbounds Base.setindex!(MW::MatrixWorkspace, x, i::Integer) =
    setindex!(MW.A, x, i)
@propagate_inbounds Base.setindex!(MW::MatrixWorkspace, x, i::Integer, j::Integer) =
    setindex!(MW.A, x, i, j)
@propagate_inbounds Base.copyto!(MW::MatrixWorkspace, A::AbstractArray) = copyto!(MW.A, A)

"""
    updated!(MW::MatrixWorkspace)

Indicate that the matrix `MW` got updated.
"""
updated!(MW::MatrixWorkspace) = (MW.factorized[] = false; MW)

"""
    update!(MW::MatrixWorkspace, A::Matrix)

Update the matrix in `MW` with `A`.
"""
function update!(MW::MatrixWorkspace, A::Matrix)
    @boundscheck (size(MW.A) == size(A) || throw(ArgumentError("Matrix of invalid size.")))
    @inbounds copyto!(MW.A, A)
    updated!(MW)
end

"""
    factorization!(MW::MatrixWorkspace, fact::Factorization)

Update the used factorization in `MW` to fact.
"""
function factorization!(MW::MatrixWorkspace, fact::Factorization)
    if fact == LU_FACT && MW.lu === nothing
        throw(ArgumentError("Cannot set `LU_FACT` for non-square matrix."))
    end
    MW.factorized[] = MW.fact[] == fact
    MW.fact[] = fact
    MW
end

##################
# Factorizations #
##################

###########################
# CUSTOM LU Factorization #
###########################
# This is an adoption of LA.generic_lufact!
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
function lu!(
    A::AbstractMatrix{T},
    b::Union{AbstractVector,Nothing} = nothing,
    ipiv::Union{Vector{I},Nothing} = nothing,
    ::Val{Pivot} = Val(true),
) where {T,I<:Integer,Pivot}
    m, n = size(A)
    minmn = min(m, n)
    # LU Factorization
    @inbounds begin
        for k = 1:minmn
            # find index max
            kp = k
            if Pivot
                amax = zero(real(T))
                for i = k:m
                    if T <: Complex
                        absi = abs2(A[i, k])
                    else
                        absi = abs(A[i, k])
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
            if !iszero(A[kp, k])
                if k != kp
                    # Interchange
                    for i = 1:n
                        tmp = A[k, i]
                        A[k, i] = A[kp, i]
                        A[kp, i] = tmp
                    end
                    if !(b isa Nothing)
                        b[k], b[kp] = b[kp], b[k]
                    end
                end
                # Scale first column
                # Akkinv = @fastmath inv(A[k,k])
                Akkinv = inv(A[k, k])
                for i = k+1:m
                    A[i, k] *= Akkinv
                end
            end
            # Update the rest
            for j = k+1:n
                A_kj = A[k, j]
                for i = k+1:m
                    A[i, j] -= A[i, k] * A_kj
                end
            end
        end
    end
    A
end

#############
# Custom QR #
#############
"""
    geqp3!(A::AbstractMatrix{ComplexF64},
        jpvt::AbstractVector{BlasInt},
        tau::AbstractVector{ComplexF64},
        work::Vector{ComplexF64},
        rwork=Vector{Float64}(undef, 2*size(A,2)))

Version of `LAPACK.geqp3!` which also accepts the work buffers.
"""
function geqp3!(
    A::AbstractMatrix{ComplexF64},
    jpvt::AbstractVector{BlasInt},
    tau::AbstractVector{ComplexF64},
    work::Vector{ComplexF64},
    rwork = Vector{Float64}(undef, 2 * size(A, 2)),
)
    m, n = size(A)
    lda = stride(A, 2)
    lda == 0 && return nothing
    jpvt .= BlasInt(0)
    lwork = BlasInt(-1)
    info = Ref{BlasInt}()
    for i = 1:2
        ccall(
            (LA.LAPACK.@blasfunc(zgeqp3_), LA.LAPACK.liblapack),
            Cvoid,
            (
             Ref{BlasInt},
             Ref{BlasInt},
             Ptr{ComplexF64},
             Ref{BlasInt},
             Ptr{BlasInt},
             Ptr{ComplexF64},
             Ptr{ComplexF64},
             Ref{BlasInt},
             Ptr{Float64},
             Ptr{BlasInt},
            ),
            m,
            n,
            A,
            lda,
            jpvt,
            tau,
            work,
            lwork,
            rwork,
            info,
        )
        LA.LAPACK.chklapackerror(info[])
        if i == 1
            lwork = BlasInt(real(work[1]))
            resize!(work, lwork)
        end
    end
    nothing
end

function factorize!(WS::MatrixWorkspace)
    if WS.fact[] == LU_FACT
        @inbounds copyto!(WS.lu.factors, WS.A)
        lu!(WS.lu.factors, nothing, WS.lu.ipiv)
    elseif WS.fact[] == QR_FACT
        @inbounds copyto!(WS.qr.factors, WS.A)
        geqp3!(WS.qr.factors, WS.qr.jpvt, WS.qr.τ, WS.qr_work, WS.qr_rwork)
    end
    WS.factorized[] = true
    WS
end


##########
## ldiv ##
##########
@inline _ipiv!(A::LA.LU, b::AbstractVector) = apply_ipiv!(b, 1:length(A.ipiv), A.ipiv)
@inline _inverse_ipiv!(A::LA.LU, b::StridedVecOrMat) =
    apply_ipiv!(b, length(A.ipiv):-1:1, A.ipiv)
@inline function apply_ipiv!(b::AbstractVector, range::OrdinalRange, ipiv)
    @inbounds for i in range
        if i != ipiv[i]
            _swap_rows!(b, i, ipiv[i])
        end
    end
    b
end
@inline function _swap_rows!(B::StridedVector, i::Integer, j::Integer)
    B[i], B[j] = B[j], B[i]
    B
end


@inline function ldiv_upper!(
    A::AbstractMatrix,
    b::AbstractVector,
    x::AbstractVector = b;
    singular_exception::Bool = false,
)
    n = size(A, 2)
    @inbounds for j = n:-1:1
        singular_exception && iszero(A[j, j]) && throw(LA.SingularException(j))
        # xj = x[j] = (@fastmath A[j,j] \ b[j])
        xj = x[j] = A[j, j] \ b[j]
        for i = 1:(j-1)
            b[i] -= A[i, j] * xj
        end
    end
    b
end
@inline function ldiv_unit_lower!(
    A::AbstractMatrix,
    b::AbstractVector,
    x::AbstractVector = b,
)
    n = size(A, 2)
    @inbounds for j = 1:n
        xj = x[j] = b[j]
        for i = j+1:n
            b[i] -= A[i, j] * xj
        end
    end
    x
end

function lu_ldiv!(x, LU::LA.LU, b::AbstractVector; check::Bool = false)
    x === b || copyto!(x, b)
    _ipiv!(LU, x)
    ldiv_unit_lower!(LU.factors, x)
    ldiv_upper!(LU.factors, x; singular_exception = check)
    x
end


function ldiv_adj_unit_lower!(A::AbstractMatrix, b::AbstractVector, x::AbstractVector = b)
    n = size(A, 1)
    @inbounds for j = n:-1:1
        z = b[j]
        for i = n:-1:j+1
            z -= conj(A[i, j]) * x[i]
        end
        x[j] = z
    end
    x
end

function ldiv_adj_upper!(
    A::AbstractMatrix,
    b::AbstractVector,
    x::AbstractVector = b;
    singular_exception::Bool = false,
)
    n = size(A, 1)
    @inbounds for j = 1:n
        z = b[j]
        for i = 1:j-1
            z -= conj(A[i, j]) * x[i]
        end
        iszero(A[j, j]) && singular_exception && throw(SingularException(j))
        x[j] = @fastmath conj(A[j, j]) \ z
    end
    x
end

function lu_ldiv_adj!(x, LU::LA.LU, b::AbstractVector; check::Bool = false)
    x === b || copyto!(x, b)
    ldiv_adj_upper!(LU.factors, x; singular_exception = check)
    ldiv_adj_unit_lower!(LU.factors, x)
    _inverse_ipiv!(LU, x)
    x
end

"""
    ormqr!(side, trans, A, tau, C, work)

Version of `LAPACK.ormqr!` which also accepts the work buffer.
Computes `Q * C` (`trans = N`), `transpose(Q) * C` (`trans = T`), `adjoint(Q) * C`
(`trans = C`) for `side = L` or the equivalent right-sided multiplication
for `side = R` using `Q` from a `QR` factorization of `A` computed using
`geqrf!`. `C` is overwritten.
"""
function ormqr!(
    side::AbstractChar,
    trans::AbstractChar,
    A::Matrix{ComplexF64},
    tau::Vector{ComplexF64},
    C::Vector{ComplexF64},
    work::Vector{ComplexF64},
)
    m, n = (size(C, 1), 1)
    mA = size(A, 1)
    k = length(tau)
    lwork = BlasInt(-1)
    info = Ref{BlasInt}()
    for i = 1:2  # first call returns lwork as work[1]
        ccall(
            (LA.LAPACK.@blasfunc(zunmqr_), LA.LAPACK.liblapack),
            Cvoid,
            (
             Ref{UInt8},
             Ref{UInt8},
             Ref{BlasInt},
             Ref{BlasInt},
             Ref{BlasInt},
             Ptr{ComplexF64},
             Ref{BlasInt},
             Ptr{ComplexF64},
             Ptr{ComplexF64},
             Ref{BlasInt},
             Ptr{ComplexF64},
             Ref{BlasInt},
             Ptr{BlasInt},
            ),
            side,
            trans,
            m,
            n,
            k,
            A,
            max(1, stride(A, 2)),
            tau,
            C,
            max(1, stride(C, 2)),
            work,
            lwork,
            info,
        )
        LA.LAPACK.chklapackerror(info[])
        if i == 1
            lwork = BlasInt(real(work[1]))
            resize!(work, lwork)
        end
    end
    C
end

# This is an adoption of the Julia implementation of ldiv! for QRPivoted
function qr_ldiv!(
    x,
    A::LA.QRPivoted{ComplexF64},
    b::Vector{ComplexF64},
    perm::Vector{Int},
    ormqr_work::Vector{ComplexF64},
      #corank::Int,
)
    mA, nA = size(A.factors)
    nr = min(mA, nA)
    nrhs = length(b)
    nr == 0 && return x

    # The following is equivalent to LA.lmul!(LA.adjoint(A.Q), b) but
    # we can also pass the preallocated work vector
    ormqr!('L', 'C', A.factors, A.τ, b, ormqr_work)
    ldiv_upper!(A.factors, b; singular_exception = false)
    @inbounds for i = 1:nr
        x[i] = b[i]
        perm[i] = A.p[i]
    end
    Base.invpermute!!(x, perm)
    # The following would cover the case of rank deficient qr.
    # This is currently not used but we keep it here since it will be used later.
    # if corank == nr
    #     fill!(x, zero(eltype(x)))
    #     return x
    # end
    # rank = nr - corank
    # if rank ≠ nA || rank ≠ nr
    #     C, τ = LA.LAPACK.tzrzf!(A.factors[1:rank,:])
    #     LA.ldiv!(LA.UpperTriangular(C[1:rank,1:rank]),view(LA.lmul!(LA.adjoint(A.Q), b), 1:rank))
    #     for i in 1:rank
    #         x[i] = b[i]
    #     end
    #     for i in rank+1:nA
    #         x[i] = zero(eltype(x))
    #     end
    #     LA.LAPACK.ormrz!('L', eltype(x)<:Complex ? 'C' : 'T', C, τ, reshape(x, length(x), 1))
    #     invpermute!(x, A.p)
    # end
    return x
end

function LA.ldiv!(x::AbstractVector, WS::MatrixWorkspace, b::AbstractVector)
    WS.factorized[] || factorize!(WS)
    if WS.fact[] == LU_FACT
        lu = WS.lu
        if lu !== nothing
            lu_ldiv!(x, lu, b)
        end
    elseif WS.fact[] == QR_FACT
        copyto!(WS.qr_rhs, b)
    # if issquare(WS.A)
            # # we applied row scaling for square matrices
            # @inbounds for i in eachindex(Jac.b)
            #     Jac.b[i] = b[i] / Jac.D[i]
            # end
    #
        qr_ldiv!(x, WS.qr, WS.qr_rhs, WS.qr_perm, WS.qr_ormqr_work)
    end
    x
end


##########################
## Iterative refinement ##
##########################
"""
    residual!(u, A, x, b, ::Type{T}=eltype(u))

Compute the residual `Ax-b` in precision `T` and store in `u`.
"""
function residual!(
    u::AbstractVector{S},
    A::AbstractMatrix{S},
    x::AbstractVector{S},
    b::AbstractVector{S},
    ::Type{T} = eltype(u),
) where {S,T}
    @boundscheck size(A, 1) == length(b) && size(A, 2) == length(x)
    m, n = size(A)
    @inbounds for i = 1:m
        dot = A[i, 1] * convert(T, x[1])
        for j = 2:n
            dot = muladd(A[i, j], convert(T, x[j]), dot)
        end
        u[i] = dot - b[i]
    end
    u
end

"""
    iterative_refinement_step!([x,] MW::MatrixWorkspace, x̂, b, norm::AbstractNorm=InfNorm(), T=eltype(x̂))

Apply one step of iterative refinement where the residual is computed with precision `T`.
Stores the result in `x` and returns the norm of the update `δx`. If `x` is not provided
`x̂` is updated inplace.
"""
function iterative_refinement_step!(
    WS::MatrixWorkspace,
    x̂,
    b,
    norm::AbstractNorm = InfNorm(),
    ::Type{T} = eltype(x̂),
) where {T}
    iterative_refinement_step!(x̂, WS, x̂, b, norm, T)
end
function iterative_refinement_step!(
    x,
    WS::MatrixWorkspace,
    x̂,
    b,
    norm::AbstractNorm = InfNorm(),
    ::Type{T} = eltype(x̂),
) where {T}
    residual!(WS.ir_r, WS.A, x̂, b, T)
    δx = LA.ldiv!(WS.ir_δx, WS, WS.ir_r)
    for i in eachindex(x)
        x[i] = x̂[i] - δx[i]
    end
    return norm(δx)
end

#########################
## Condition estimates ##
#########################
"""
    rcond(WS::MatrixWorkspace)

Computes the reciprocal condition number of the matrix `A` in `WS` wrt the ∞-norm.
If `A` is a non-square matrix, the reciprocal condition number of `R` in the
QR decomposition is computed.
"""
function rcond(WS::MatrixWorkspace)
    WS.factorized[] || factorize!(WS)
    if WS.fact[] == LU_FACT
        lu = WS.lu
        if lu === nothing
            NaN
        else
            rcond(lu, LA.opnorm(WS.A, Inf), WS.lu_cond_work, WS.lu_cond_rwork)
        end
    else # WS.fact[] == QR_FACT
        rcond(WS.qr, WS.qr_cond_work, WS.qr_cond_rwork)
    end
end
function rcond(LU::LA.LU, anorm, work, rwork)
    A = LU.factors
    n = size(A, 1)
    rcond = Ref{Float64}()
    info = Ref{Int64}()
    lda = max(1, stride(A, 2))
    ccall(
        (LA.LAPACK.@blasfunc(zgecon_), LinearAlgebra.LAPACK.liblapack),
        Cvoid,
        (
         Ref{UInt8},
         Ref{Int64},
         Ptr{ComplexF64},
         Ref{Int64},
         Ref{Float64},
         Ref{Float64},
         Ptr{ComplexF64},
         Ptr{Float64},
         Ptr{Int64},
        ),
        'I',
        n,
        A,
        lda,
        anorm,
        rcond,
        work,
        rwork,
        info,
    )
    rcond[]
end

function rcond(qr::LA.QRPivoted, work, rwork)
    # use trcon to estimate cond of R
    m, n = size(qr.factors)
    rcond = Ref{Float64}(1)
    info = Ref{BlasInt}()
    ccall(
        (LA.LAPACK.@blasfunc(ztrcon_), LA.LAPACK.liblapack),
        Cvoid,
        (
         Ref{UInt8},
         Ref{UInt8},
         Ref{UInt8},
         Ref{BlasInt},
         Ptr{ComplexF64},
         Ref{BlasInt},
         Ref{Float64},
         Ptr{ComplexF64},
         Ptr{Float64},
         Ptr{BlasInt},
        ),
        'I',
        'U',
        'N',
        n,
        qr.factors,
        m,
        rcond,
        work,
        rwork,
        info,
    )
    if m == n
        # Q has condition number 1 in the 2 norm
        # Some experiments indicate that in the Inf norm this is roughly
        # 0.8√n
        rcond[] / (0.8 * √(n))
    else
        rcond[]
    end
end


"""
    inf_norm_est(WS::MatrixWorkspace, g::Union{Nothing,Vector{<:Real}}=nothing, d::Union{Nothing,Vector{<:Real}}=nothing)
    inf_norm_est(lu::LA.LU, g::Union{Nothing,Vector{<:Real}}=nothing, d::Union{Nothing,Vector{<:Real}}=nothing)
    inf_norm_est(lu::LA.LU, g::Union{Nothing,Vector{<:Real}}, d::Union{Nothing,Vector{<:Real}}, work{<:Complex}, rwork::Vector{<:Real})

Estimation of the infinity norm of `diag(d)⁻¹A⁻¹g` where `g` and `d` are optional positive vectors
and `A=LU`. If `g` is `nothing` the infinity norm of `A⁻¹` is estimated.
This uses the 1-norm lapack condition estimator described by Highahm in [^H88].

[^H88]: Higham, Nicholas J. "FORTRAN codes for estimating the one-norm of a real or complex matrix, with applications to condition estimation." ACM Transactions on Mathematical Software (TOMS) 14.4 (1988): 381-396.
"""
function inf_norm_est(
    lu::LA.LU,
    g::Union{Nothing,Vector{<:Real}} = nothing,
    d::Union{Nothing,Vector{<:Real}} = nothing,
)
    n = size(lu.factors, 2)
    work = Vector{eltype(lu.factors)}(undef, n)
    rwork = Vector{real(eltype(lu.factors))}(undef, n)
    inf_norm_est(lu, g, d, work, rwork)
end
function inf_norm_est(
    WS::MatrixWorkspace,
    g::Union{Nothing,Vector{<:Real}} = nothing,
    d::Union{Nothing,Vector{<:Real}} = nothing,
)
    WS.fact[] == LU_FACT || factorization!(WS, LU_FACT)
    WS.factorized[] || factorize!(WS)
    inf_norm_est(WS.lu, g, d, WS.inf_norm_est_work, WS.inf_norm_est_rwork)
end
function inf_norm_est(
    lu::LA.LU,
    g::Union{Nothing,Vector{<:Real}},
    d::Union{Nothing,Vector{<:Real}},
    work::Vector{<:Complex},
    rwork::Vector{<:Real},
)
    z = ξ = y = work
    x = rwork

    n = size(lu.factors, 1)
    x .= inv(n)
    if d !== nothing
        x ./= d
    end
    lu_ldiv_adj!(y, lu, x)
    if g !== nothing
        y .*= g
    end
    γ = sum(abs, y)
    y ./= abs.(y)
    if g !== nothing
        y .*= g
    end
    lu_ldiv!(z, lu, ξ)
    if d !== nothing
        x .= real.(z) ./ d
    else
        x .= real.(z)
    end
    k = 2
    while true
        j, max_xᵢ = 1, abs(x[1])
        for i = 2:n
            abs_xᵢ = abs(x[i])
            if abs_xᵢ > max_xᵢ
                j, max_xᵢ = i, abs_xᵢ
            end
        end
        x .= zero(eltype(x))
        x[j] = one(eltype(x))
        if d !== nothing
            x ./= d
        end
        lu_ldiv_adj!(y, lu, x)
        if g !== nothing
            y .*= g
        end
        γ̄ = γ
        γ = sum(abs, y)
        γ ≤ γ̄ && break
        ξ .= y ./ abs.(y)
        if g !== nothing
            ξ .*= g
        end
        lu_ldiv!(z, lu, ξ)
        if d !== nothing
            x .= real.(z) ./ d
        else
            x .= real.(z)
        end
        k += 1
        if x[j] == LA.norm(x, Inf) || k > 5
            break
        end
    end
    γ
end

"""
    row_scaling!(WS::MatrixWorkspace,
        norm::Union{InfNorm, WeightedNorm{<:InfNorm}},
        residual::Union{Nothing,AbstractVector{<:Real}} = nothing)

Compute a suitable scaling `D` of the rows of `WS.A` and store it in `WS.row_scaling`.
The returned scaling has the property that the from `norm` induced operator for `D⁻¹ * WS.A`
is 1.
The vector `residual` is a an optional vector containing the absolute values of the limiting
residual of Newton's method. The scaling is performaned in such a way that `norm(residual)`
is *not* increase. Unless the from `norm` induced operator norm is less than 1.
Then the residual is increased by the minimal amount necessary.
"""
function row_scaling!(
    WS::MatrixWorkspace{T},
    norm::Union{InfNorm,WeightedNorm{<:InfNorm}},
    r::Union{Nothing,AbstractVector{<:Real}} = nothing,
) where {T}
    d = WS.row_scaling
    d .= one(T)
    m = length(d)
    if isa(norm, WeightedNorm)
        for j = 1:size(WS.A, 2), i = 1:m
            d[i] += abs(WS.A[i, j]) * norm[j]
        end
    else
        for j = 1:size(WS.A, 2), i = 1:m
            d[i] += abs(WS.A[i, j])
        end
    end

    if r !== nothing
        # d̂ = clamp(maximum(d), 1e-4, 1.0)
        d̂ = 1e-4
        for dᵢ in d
            dᵢ ≥ 1.0 && (d̂ = 1.0; break)
            d̂ = max(d̂, dᵢ)
        end

        # compute the norm of r under the scaling of all rows whose norm is larger
        # than d̂, i.e. all rows that either decrease the value of |rᵢ| or
        # or increase it by the minimal abount to still have opnorm 1.
        norm_scaled_r = r[1] / max(d[1], d̂)
        for i = 2:m
            norm_scaled_r = max(norm_scaled_r, r[i] / max(d[i], d̂))
        end
        for i = 1:length(d)
            r̂ᵢ = r[i] / norm_scaled_r
            if isfinite(r̂ᵢ)
                d[i] = max(d[i], r̂ᵢ)
            end
        end
    else
        # Don't scale zero rows -> set dᵢ = 1 if too small
        d_max_eps = sqrt(eps(min(maximum(d), one(T))))
        @inbounds for i = 1:m
            if d[i] < d_max_eps
                d[i] = one(T)
            end
        end
    end
    d
end

"""
    scaled_cond!(WS::MatrixWorkspace, norm, residual)

Perform first a row scaling as described by [`row_scaling!`](@ref) and then
compute the condition number ``κ`` of the row-equilibrated matrix `WS.A`.
This is the best possible condition number under row scaling (e.g. [Higham02, Thm 7.5]).

[Higham02]: Higham, Nicholas J. Accuracy and stability of numerical algorithms. Vol. 80. Siam, 200
"""
function scaled_cond!(
    WS::MatrixWorkspace,
    norm::Union{<:WeightedNorm{<:InfNorm},InfNorm},
    r::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    d = row_scaling!(WS, norm, r)
    # row scaling makes ||A||_∞ = 1 so we don't this here
    if (norm isa WeightedNorm)
        cond_est = inf_norm_est(WS, d, weights(norm))
    else
        cond_est = inf_norm_est(WS, d)
    end
    cond_est
end

#####################
## JacobianMonitor ##
#####################
struct JacobianMonitor{T}
    J::MatrixWorkspace{T}
    cond::Base.RefValue{Float64}
    forward_err::Base.RefValue{Float64}
    # stats
    factorizations::Base.RefValue{Int}
    ldivs::Base.RefValue{Int}
end
function JacobianMonitor(A::AbstractMatrix)
    J = MatrixWorkspace(A)
    cond = Ref(1.0)
    forward_err = Ref(0.0)
    JacobianMonitor(J, cond, forward_err, Ref(0), Ref(0))
end

updated!(JM::JacobianMonitor) = updated!(JM.J)
jacobian(JM::JacobianMonitor) = JM.J

function Base.show(io::IO, JM::JacobianMonitor{T}) where {T}
    println(io, "JacobianMonitor{$T}:")
    println(io, " • cond → ", round(JM.cond[], sigdigits = 5))
    println(io, " • forward_err → ", round(JM.forward_err[], sigdigits = 5))
    println(io, " • # factorizations → ", JM.factorizations[])
    println(io, " • # ldivs → ", JM.ldivs[])
end

"""
    init!(JM::JacobianMonitor)

(Re-)initialize the `JacobianMonitor`.
"""
function init!(JM::JacobianMonitor; keep_stats::Bool = false)
    JM.cond[] = 1.0
    JM.forward_err[] = 0.0
    JM.factorizations[] = 0.0
    JM.ldivs[] = 0.0
    JM
end

"""
    forward_err(JM::JacobianMonitor)

Return the last estimate of the forward error as computed by [`forward_err`](@ref) or [`strong_forward_err`](@ref).
"""
forward_err(JM::JacobianMonitor) = JM.forward_err[]

"""
    forward_err!(x̂::AbstractVector, JM::JacobianMonitor, b::AbstractVector, norm::AbstractNorm)

Compute an estimate of the forward error ||x-x̂||/||x|| and update the value in `JM`.
"""
function forward_err!(
    JM::JacobianMonitor,
    x̂::AbstractVector,
    b::AbstractVector,
    norm::AbstractNorm,
    T = eltype(x̂),
)
    forward_err!(x̂, JM, x̂, b, norm, T)
end
function forward_err!(
    x::AbstractVector,
    JM::JacobianMonitor,
    x̂::AbstractVector,
    b::AbstractVector,
    norm::AbstractNorm,
    T = eltype(x̂),
)
    norm_x̂ = norm(x̂)
    JM.forward_err[] = iterative_refinement_step!(x, JM.J, x̂, b, norm, T) / norm_x̂
end

"""
    strong_forward_err!(x̂::AbstractVector, JM::JacobianMonitor, b::AbstractVector, norm::AbstractNorm)

Compute a more robust estimate of the forward error ||x-x̂||/||x|| using eq. (2.14) in [Demmel, Section 2.4.4]
and update the value in `JM`.
If `jacobian(JM).fact[] == QR_FACT` then this falls back to [`forward_err!`](@ref).

[Demmel]: Demmel, James W. Applied numerical linear algebra. Vol. 56. Siam, 1997.
"""
function strong_forward_err!(
    JM::JacobianMonitor,
    x̂::AbstractVector,
    b::AbstractVector,
    norm::InfNorm,
    T = eltype(x̂),
)
    strong_forward_err!(x̂, JM, x̂, b, norm, T)
end
function strong_forward_err!(
    x::AbstractVector,
    JM::JacobianMonitor,
    x̂::AbstractVector,
    b::AbstractVector,
    norm::Union{InfNorm,WeightedNorm{InfNorm}},
    T = eltype(x̂),
)
    if jacobian(JM).fact[] == QR_FACT
        return JM.forward_err[] = forward_err!(x, JM, x̂, b, norm)
    end
    norm_x̂ = norm(x̂)
    WS = jacobian(JM)
    residual!(WS.ir_r, WS.A, x̂, b, T)
    WS.residual_abs .= abs.(b)
    for j = 1:size(WS.A, 2)
        ax̂ⱼ = abs(x̂[j])
        for i = 1:size(WS.A, 1)
            WS.residual_abs[i] += abs(WS.A[i, j]) * ax̂ⱼ
        end
    end
    WS.residual_abs .*= eps()
    WS.residual_abs .+= abs.(WS.ir_r)

    if norm isa WeightedNorm
        JM.forward_err[] = inf_norm_est(WS, WS.residual_abs, weights(norm)) / norm_x̂
    else
        JM.forward_err[] = inf_norm_est(WS, WS.residual_abs) / norm_x̂
    end
end

"""
    cond!(JM::JacobianMonitor, norm, residual)

Compute the condition number ``κ`` of the Jacobian ``J`` with respect to the
infinity norm and store it. If ``J`` is a square matrix the condition number of the
row-equilibrated Jacobian is computed. See [`scaled_cond!`](@ref) for details.
"""
function cond!(
    JM::JacobianMonitor,
    norm::AbstractNorm = InfNorm(),
    r::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    if jacobian(JM).fact[] == LU_FACT
        JM.cond[] = scaled_cond!(jacobian(JM), norm, r)
    else
        # TODO: The weighted norm is NOT considered currently.
        JM.cond[] = inv(rcond(jacobian(JM)))
    end
end

"""
    apply_row_scaling!(r, JM::JacobianMonitor)

Apply the computed row scaling used in [`cond!`](@ref) resp. [`scaled_cond!`](@ref) to `r`.
"""
apply_row_scaling!(r::AbstractVector, JM::JacobianMonitor) = (r ./= JM.J.row_scaling; r)

"""
    cond(JM::JacobianMonitor)

Returns the with [`cond!`](@ref) computed condition number.
"""
LA.cond(JM::JacobianMonitor) = JM.cond[]

"""
    enum JacobianMonitorUpdates

## Cases

* `JAC_MONITOR_UPDATE_NOTHING`: Do nothing
* `JAC_MONITOR_UPDATE_FERR`: Update the forward error estimate
* `JAC_MONITOR_UPDATE_COND`: Update the condition number estimate
* `JAC_MONITOR_UPDATE_ALL`: Update the forward error and condition number estimate
"""
@enum JacobianMonitorUpdates begin
    JAC_MONITOR_UPDATE_NOTHING
    JAC_MONITOR_UPDATE_FERR
    JAC_MONITOR_UPDATE_COND
    JAC_MONITOR_UPDATE_ALL
end

"""
    ldiv!(x̂::AbstractVector,
                  JM::JacobianMonitor,
                  b::AbstractVector,
                  norm::AbstractNorm,
                  update::JacobianMonitorUpdates=JAC_MONITOR_UPDATE_NOTHING)

solve the linear system `jacobian(JM)x̂=b`. `update` controls the computation of additional
informations, see [`JacobianMonitorUpdates`](@ref).
"""
function LA.ldiv!(
    x̂::AbstractVector,
    JM::JacobianMonitor,
    b::AbstractVector,
    norm::AbstractNorm = InfNorm(),
    update::JacobianMonitorUpdates = JAC_MONITOR_UPDATE_NOTHING,
)
    # stats update
    JM.factorizations[] += !jacobian(JM).factorized[]
    JM.ldivs[] += 1

    LA.ldiv!(x̂, jacobian(JM), b)
    if update == JAC_MONITOR_UPDATE_FERR || update == JAC_MONITOR_UPDATE_ALL
        forward_err!(JM, x̂, b, norm)
    end
    if update == JAC_MONITOR_UPDATE_COND || update == JAC_MONITOR_UPDATE_ALL
        cond!(JM)
    end
    x̂
end

#########################
## Hermite Normal Form ##
#########################

"""
    hnf(A, T=elytpe(A))

Compute the hermite normal form `H` of `A` by overwriting `A` and the corresponding transformation
matrix `U` using precision T. This is `A*U == H` and `H` is a lower triangular matrix.

The implementation follows the algorithm described in [1].

[1] Kannan, Ravindran, and Achim Bachem. "Polynomial algorithms for computing the Smith and Hermite normal forms of an integer matrix." SIAM Journal on Computing 8.4 (1979): 499-507.
"""
function hnf(A, T = eltype(A))
    H = similar(A, T)
    H .= A
    U = similar(A, T)
    hnf!(H, U)
    H, U
end

# use checked arithmethic
⊡(x, y) = Base.checked_mul(x, y)
⊞(x, y) = Base.checked_add(x, y)

"""
    hnf!(H, U, A)

Inplace version of [hnf](@ref).

    hnf!(A, U)

Inplace version of [hnf](@ref) overwriting `A` with `H`.
"""
hnf!(H, U, A) = hnf!(copyto!(H, A), U)
function hnf!(A, U)
    n = size(A, 1)
    U .= 0
    @inbounds for i = 1:n
        U[i, i] = one(eltype(U))
    end
    @inbounds for i = 1:(n-1)
        ii = i ⊞ 1
        for j = 1:i
            if !iszero(A[j, j]) || !iszero(A[j, ii])
                # 4.1
                r, p, q = gcdx(A[j, j], A[j, ii])
                # 4.2
                d_j = -A[j, ii] ÷ r
                d_ii = A[j, j] ÷ r
                for k = 1:n
                    a_kj, a_kii = A[k, j], A[k, ii]
                    A[k, j] = a_kj ⊡ p ⊞ a_kii ⊡ q
                    A[k, ii] = a_kj ⊡ d_j ⊞ a_kii ⊡ d_ii

                    u_kj, u_kii = U[k, j], U[k, ii]
                    U[k, j] = u_kj ⊡ p ⊞ u_kii ⊡ q
                    U[k, ii] = u_kj ⊡ d_j ⊞ u_kii ⊡ d_ii
                end
            end
            # 4.3
            if j > 1
                reduce_off_diagonal!(A, U, j)
            end
        end
        #5
        reduce_off_diagonal!(A, U, ii)
    end
    # Special case for 1 × 1 matrices to guarantee that the diagonal is positive
    # This comes up in polyhedral homotopy for univariate polynomials
    if n == 1
        U[1, 1] = flipsign(U[1, 1], A[1, 1])
        A[1, 1] = flipsign(A[1, 1], A[1, 1])
    end

    nothing
end

@inline function reduce_off_diagonal!(A, U, k)
    n = size(A, 1)
    @inbounds if A[k, k] < 0
        for i = 1:n
            A[i, k] = -A[i, k]
            U[i, k] = -U[i, k]
        end
    end
    @inbounds for z = 1:(k-1)
        if !iszero(A[z, z])
            r = -ceil(eltype(A), A[k, z] / A[z, z])
            for i = 1:n
                U[i, z] = U[i, z] ⊞ r ⊡ U[i, k]
                A[i, z] = A[i, z] ⊞ r ⊡ A[i, k]
            end
        end
    end
    nothing
end
