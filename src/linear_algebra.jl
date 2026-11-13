"""
    MatrixWorkspace

This is a data structure for the efficient repeated solution of a square or
overdetermined linear system `Ax=b`.
"""
struct MatrixWorkspace{M<:AbstractMatrix{ComplexF64}} <: AbstractMatrix{ComplexF64}
    A::M # Matrix
    d::Vector{Float64} # Inverse of scaling factors
    factorized::Base.RefValue{Bool}
    lu::LA.LU{ComplexF64,M,Vector{Int}} # LU Factorization of D * J
    qr::LA.QR{ComplexF64,Matrix{ComplexF64},Vector{ComplexF64}}
    row_scaling::Vector{Float64}
    scaled::Base.RefValue{Bool}
    # mixed precision iterative refinement
    x̄::Vector{ComplexDF64}
    r::Vector{ComplexF64}
    r̄::Vector{ComplexDF64}
    δx::Vector{ComplexF64}
    inf_norm_est_work::Vector{ComplexF64}
    inf_norm_est_rwork::Vector{Float64}
    # real and imaginary parts of the LU factors and of a solve (square matrices)
    lu_re::Matrix{Float64}
    lu_im::Matrix{Float64}
    x_re::Vector{Float64}
    x_im::Vector{Float64}
    # whether `lu.factors` holds the current split factors
    lu_synced::Base.RefValue{Bool}
end

MatrixWorkspace(m::Integer, n::Integer; kwargs...) =
    MatrixWorkspace(zeros(ComplexF64, m, n); kwargs...)
function MatrixWorkspace(Â::AbstractMatrix; optimize_data_structure = true)
    m, n = size(Â)
    m ≥ n || throw(ArgumentError("Expected system with more rows than columns."))

    A = Matrix{ComplexF64}(Â)
    d = ones(m)
    factorized = Ref(false)
    qr = LA.qrfactUnblocked!(copy(A))
    # experiments show that for m > 25 the data layout as a
    # struct array is beneficial
    if m > 25 && optimize_data_structure
        A = StructArrays.StructArray(A)
    end
    row_scaling = ones(m)
    scaled = Ref(false)

    ipiv = zeros(Int, m)
    lu = LinearAlgebra.LU{eltype(A),typeof(A),typeof(ipiv)}(copy(A), ipiv, 0)
    r = zeros(ComplexF64, m)
    r̄ = zeros(ComplexDF64, m)
    x̄ = zeros(ComplexDF64, n)
    δx = zeros(ComplexF64, n)
    inf_norm_est_work = Vector{ComplexF64}(undef, n)
    inf_norm_est_rwork = Vector{Float64}(undef, n)
    k = m == n ? n : 0

    MatrixWorkspace(
        A,
        d,
        factorized,
        lu,
        qr,
        row_scaling,
        scaled,
        x̄,
        r,
        r̄,
        δx,
        inf_norm_est_work,
        inf_norm_est_rwork,
        zeros(k, k),
        zeros(k, k),
        zeros(k),
        zeros(k),
        Ref(true),
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

matrix(M::MatrixWorkspace) = M.A
matrix(M::AbstractMatrix) = M

"""
    updated!(MW::MatrixWorkspace)

Indicate that the matrix `MW` got updated.
"""
function updated!(MW::MatrixWorkspace)
    MW.factorized[] = false
    MW.scaled[] = false
    m, n = size(MW)
    # square matrices are read from `A` by `split_lu!`
    m == n || @inbounds copyto!(MW.qr.factors, MW.A)
    MW
end
updated!(M::AbstractMatrix) = M

"""
    update!(MW::MatrixWorkspace, A::Matrix)

Update the matrix in `MW` with `A`.
"""
@inline function update!(MW::AbstractMatrix, A::Matrix)
    @boundscheck (size(MW.A) == size(A) || throw(ArgumentError("Matrix of invalid size.")))
    @inbounds copyto!(MW.A, A)
    updated!(MW)
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
    ipiv::Union{Vector{I},Nothing} = nothing,
) where {T,I<:Integer}
    m, n = size(A)
    if m == n && T === ComplexF64
        n == 3 && return lu_small!(A, ipiv, Val(3))
        n == 4 && return lu_small!(A, ipiv, Val(4))
    end
    min(m, n) >= 16 && return lu_rank2!(A, ipiv)
    lu_rank1!(A, ipiv)
end

function lu_rank1!(A::AbstractMatrix{T}, ipiv = nothing) where {T}
    m, n = size(A)
    minmn = min(m, n)
    # LU Factorization
    @inbounds begin
        for k = 1:minmn
            # find index max
            kp = k
            if T <: Complex
                amax = abs2(A[k, k])
            else
                amax = abs(A[k, k])
            end
            for i = k+1:m
                if T <: Complex
                    absi = abs2(A[i, k])
                else
                    absi = abs(A[i, k])
                end
                kp, amax = ifelse(absi > amax, (i, absi), (kp, amax))
            end
            if !(ipiv isa Nothing)
                ipiv[k] = kp
            end
            if !iszero(amax)
                if k != kp
                    # Interchange
                    for i = 1:n
                        tmp = A[k, i]
                        A[k, i] = A[kp, i]
                        A[kp, i] = tmp
                    end
                end
                # Scale first column
                Akk = A[k, k]
                for i = k+1:m
                    # we assume the matrix is decently scaled
                    # so that the naive division algorithm doesn't overflow
                    @fastmath A[i, k] = A[i, k] / Akk
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

@generated function lu_small!(A, ipiv, ::Val{N}) where {N}
    body = Expr(:block)
    for k = 1:N
        push!(body.args, :(lu_pivot!(A, ipiv, $k, $N, $N)))
        for j = k+1:N, i = k+1:N
            push!(body.args, :(A[$i, $j] -= A[$i, $k] * A[$k, $j]))
        end
    end
    quote
        @inbounds $body
        A
    end
end

# Factor two columns before streaming over the trailing matrix. This retains
# the existing LU storage and pivot convention (including rectangular inputs).
@inline function lu_pivot!(A, ipiv, k, m, n)
    @inbounds begin
        kp = k
        amax = abs2(A[k, k])
        for i = k+1:m
            a = abs2(A[i, k])
            if a > amax
                kp, amax = i, a
            end
        end
        # Squared magnitudes are cheap in the usual exponent range. Fall back
        # to scaled magnitudes if they cannot distinguish the pivot candidates.
        if iszero(amax) || !isfinite(amax)
            amax = abs(A[k, k])
            kp = k
            for i = k+1:m
                a = abs(A[i, k])
                if a > amax
                    kp, amax = i, a
                end
            end
        end
        ipiv === nothing || (ipiv[k] = kp)
        iszero(amax) && return false
        if kp != k
            for j = 1:n
                A[k, j], A[kp, j] = A[kp, j], A[k, j]
            end
        end
        pivot = A[k, k]
        reciprocal = inv(pivot)
        if isfinite(reciprocal)
            for i = k+1:m
                A[i, k] *= reciprocal
            end
        else
            for i = k+1:m
                A[i, k] /= pivot
            end
        end
    end
    true
end

function lu_rank2!(A, ipiv)
    m, n = size(A)
    @inbounds for k = 1:2:min(m, n)
        lu_pivot!(A, ipiv, k, m, n)
        if k == min(m, n)
            for j = k+1:n, i = k+1:m
                A[i, j] -= A[i, k] * A[k, j]
            end
            break
        end
        # Only the next column is needed to choose the second pivot. Swapping
        # its rows also swaps the first column's multipliers.
        for i = k+1:m
            A[i, k+1] -= A[i, k] * A[k, k+1]
        end
        lu_pivot!(A, ipiv, k + 1, m, n)
        for j = k+2:n
            a = A[k, j]
            b = A[k+1, j] - A[k+1, k] * a
            A[k+1, j] = b
            for i = k+2:m
                A[i, j] = (A[i, j] - A[i, k] * a) - A[i, k+1] * b
            end
        end
    end
    A
end

######################
# Split LU (ComplexF64)
######################
# The factorization and solves of square workspaces work on separate real and
# imaginary parts, so that the inner loops vectorize over rows (ported from wndr's
# SplitLU). Pivots are chosen as in `lu_pivot!`. `lu.factors` receives the factors
# only when a routine reads them directly (condition estimates, inverses), see
# `synced_lu`.

# From this size on, two pivots are eliminated per pass over the trailing matrix.
const SPLIT_LU_RANK2_MIN = 20

@inline function split_lu_pivot!(R, I, ipiv, k, n)
    @inbounds begin
        kp = k
        amax = R[k, k]^2 + I[k, k]^2
        for i = k+1:n
            a = R[i, k]^2 + I[i, k]^2
            kp, amax = ifelse(a > amax, (i, a), (kp, amax))
        end
        # Squared magnitudes over- or underflow outside the usual exponent range.
        squared = !(iszero(amax) || !isfinite(amax))
        if !squared
            kp = k
            amax = hypot(R[k, k], I[k, k])
            for i = k+1:n
                a = hypot(R[i, k], I[i, k])
                if a > amax
                    kp, amax = i, a
                end
            end
        end
        ipiv[k] = kp
        iszero(amax) && return false
        if kp != k
            for j = 1:n
                R[k, j], R[kp, j] = R[kp, j], R[k, j]
                I[k, j], I[kp, j] = I[kp, j], I[k, j]
            end
        end
        pivot = complex(R[k, k], I[k, k])
        # naive reciprocal conj(pivot) / |pivot|² unless |pivot|² is out of range
        ir, ii = squared ? (R[k, k] / amax, -I[k, k] / amax) : reim(inv(pivot))
        if isfinite(ir) && isfinite(ii)
            @simd ivdep for i = k+1:n
                r, m = R[i, k], I[i, k]
                R[i, k] = r * ir - m * ii
                I[i, k] = r * ii + m * ir
            end
        else
            for i = k+1:n
                R[i, k], I[i, k] = reim(complex(R[i, k], I[i, k]) / pivot)
            end
        end
    end
    true
end

# Update of column j (rows k+1:n) by pivot column k
@inline function split_lu_update!(R, I, j, k, n)
    @inbounds begin
        br, bi = R[k, j], I[k, j]
        @simd ivdep for i = k+1:n
            lr, li = R[i, k], I[i, k]
            R[i, j] = muladd(-lr, br, muladd(li, bi, R[i, j]))
            I[i, j] = muladd(-lr, bi, muladd(-li, br, I[i, j]))
        end
    end
end

function split_lu_rank1!(R, I, ipiv, n)
    @inbounds for k = 1:n
        split_lu_pivot!(R, I, ipiv, k, n) || continue
        for j = k+1:n
            split_lu_update!(R, I, j, k, n)
        end
    end
end

# Pivots k and k+1 in one pass over the trailing columns. Every entry receives the
# updates of k and k+1 in that order, as in `split_lu_rank1!`.
function split_lu_rank2!(R, I, ipiv, n)
    k = 1
    @inbounds while k ≤ n
        ok₁ = split_lu_pivot!(R, I, ipiv, k, n)
        k == n && break
        ok₁ && split_lu_update!(R, I, k + 1, k, n)
        ok₂ = split_lu_pivot!(R, I, ipiv, k + 1, n)
        if ok₁ && ok₂
            lr, li = R[k+1, k], I[k+1, k]
            for j = k+2:n
                br, bi = R[k, j], I[k, j]
                R[k+1, j] = muladd(-lr, br, muladd(li, bi, R[k+1, j]))
                I[k+1, j] = muladd(-lr, bi, muladd(-li, br, I[k+1, j]))
            end
            for j = k+2:n
                br, bi, cr, ci = R[k, j], I[k, j], R[k+1, j], I[k+1, j]
                @simd ivdep for i = k+2:n
                    lr, li = R[i, k], I[i, k]
                    mr, mi = R[i, k+1], I[i, k+1]
                    sr = muladd(-lr, br, muladd(li, bi, R[i, j]))
                    si = muladd(-lr, bi, muladd(-li, br, I[i, j]))
                    R[i, j] = muladd(-mr, cr, muladd(mi, ci, sr))
                    I[i, j] = muladd(-mr, ci, muladd(-mi, cr, si))
                end
            end
        else
            for j = k+2:n
                ok₁ && split_lu_update!(R, I, j, k, n)
            end
            for j = k+2:n
                ok₂ && split_lu_update!(R, I, j, k + 1, n)
            end
        end
        k += 2
    end
end

# Factor `D * A`, with D the row scaling if it is applied.
function split_lu!(WS::MatrixWorkspace)
    A, ipiv, R, I = WS.A, WS.lu.ipiv, WS.lu_re, WS.lu_im
    n = size(R, 1)
    if WS.scaled[]
        d = WS.row_scaling
        @inbounds for j = 1:n, i = 1:n
            R[i, j], I[i, j] = reim(A[i, j] * d[i])
        end
    else
        @inbounds for j = 1:n, i = 1:n
            R[i, j], I[i, j] = reim(A[i, j])
        end
    end
    if n ≥ SPLIT_LU_RANK2_MIN
        split_lu_rank2!(R, I, ipiv, n)
    else
        split_lu_rank1!(R, I, ipiv, n)
    end
    WS.lu_synced[] = false
    WS
end

# The LU factorization with its factors in `lu.factors`.
function synced_lu(WS::MatrixWorkspace)
    WS.factorized[] || factorize!(WS)
    if !WS.lu_synced[]
        F, R, I = WS.lu.factors, WS.lu_re, WS.lu_im
        n = size(R, 1)
        @inbounds for j = 1:n, i = 1:n
            F[i, j] = complex(R[i, j], I[i, j])
        end
        WS.lu_synced[] = true
    end
    WS.lu
end

# (br + i bi) / (ar + i ai), naive complex division as in `ldiv_upper!`
@inline function split_cdiv(br, bi, ar, ai)
    d = ar^2 + ai^2
    (br * ar + bi * ai) / d, (bi * ar - br * ai) / d
end

# Solve with the split factors. `b` is multiplied by `scaling` if given.
function split_lu_ldiv!(x, WS::MatrixWorkspace, b, scaling = nothing)
    R, I, ipiv, xr, xi = WS.lu_re, WS.lu_im, WS.lu.ipiv, WS.x_re, WS.x_im
    n = size(R, 1)
    @inbounds begin
        for i = 1:n
            bi = scaling === nothing ? b[i] : scaling[i] * b[i]
            xr[i], xi[i] = reim(bi)
        end
        for i = 1:n
            p = ipiv[i]
            if p != i
                xr[i], xr[p] = xr[p], xr[i]
                xi[i], xi[p] = xi[p], xi[i]
            end
        end
        # unit lower triangular, columns j and j+1 per pass
        j = 1
        while j < n
            br, bi = xr[j], xi[j]
            lr, li = R[j+1, j], I[j+1, j]
            cr = muladd(-lr, br, muladd(li, bi, xr[j+1]))
            ci = muladd(-lr, bi, muladd(-li, br, xi[j+1]))
            xr[j+1], xi[j+1] = cr, ci
            @simd ivdep for i = j+2:n
                lr, li = R[i, j], I[i, j]
                mr, mi = R[i, j+1], I[i, j+1]
                sr = muladd(-lr, br, muladd(li, bi, xr[i]))
                si = muladd(-lr, bi, muladd(-li, br, xi[i]))
                xr[i] = muladd(-mr, cr, muladd(mi, ci, sr))
                xi[i] = muladd(-mr, ci, muladd(-mi, cr, si))
            end
            j += 2
        end
        # upper triangular, columns j and j-1 per pass
        j = n
        while j ≥ 1
            yr, yi = split_cdiv(xr[j], xi[j], R[j, j], I[j, j])
            xr[j], xi[j] = yr, yi
            j == 1 && break
            ur, ui = R[j-1, j], I[j-1, j]
            zr, zi = split_cdiv(
                muladd(-ur, yr, muladd(ui, yi, xr[j-1])),
                muladd(-ur, yi, muladd(-ui, yr, xi[j-1])),
                R[j-1, j-1],
                I[j-1, j-1],
            )
            xr[j-1], xi[j-1] = zr, zi
            @simd ivdep for i = 1:j-2
                ur, ui = R[i, j], I[i, j]
                vr, vi = R[i, j-1], I[i, j-1]
                sr = muladd(-ur, yr, muladd(ui, yi, xr[i]))
                si = muladd(-ur, yi, muladd(-ui, yr, xi[i]))
                xr[i] = muladd(-vr, zr, muladd(vi, zi, sr))
                xi[i] = muladd(-vr, zi, muladd(-vi, zr, si))
            end
            j -= 2
        end
        for i = 1:n
            x[i] = complex(xr[i], xi[i])
        end
    end
    x
end


####################
# QR Factorization #
####################
# This is the generic base implementation with the small change that we introduce
# a @fastmath for the division step since results in an approx 50% speedup

# Elementary reflection similar to LAPACK. The reflector is not Hermitian but
# ensures that tridiagonalization of Hermitian matrices become real. See lawn72
@inline function reflector!(x::AbstractVector)
    n = length(x)
    @inbounds begin
        ξ1 = x[1]
        normu = abs2(ξ1)
        for i = 2:n
            normu += abs2(x[i])
        end
        if iszero(normu)
            return zero(ξ1 / normu)
        end
        normu = sqrt(normu)
        ν = copysign(normu, real(ξ1))
        ξ1 += ν
        x[1] = -ν
        for i = 2:n
            @fastmath x[i] = x[i] / ξ1
        end
    end
    ξ1 / ν
end

# apply reflector from left
@inline function reflectorApply!(x::AbstractVector, τ::Number, A::StridedMatrix)
    m, n = size(A)
    @inbounds begin
        for j = 1:n
            # dot
            vAj = A[1, j]
            for i = 2:m
                vAj += x[i]' * A[i, j]
            end

            vAj = conj(τ) * vAj

            # ger
            A[1, j] -= vAj
            for i = 2:m
                A[i, j] -= x[i] * vAj
            end
        end
    end
    return A
end

function qr!(qr::LA.QR{T}) where {T}
    A = qr.factors
    τ = qr.τ
    m, n = size(A)
    for k = 1:min(m - 1 + !(T <: Real), n)
        x = view(A, k:m, k)
        τk = reflector!(x)
        τ[k] = τk
        reflectorApply!(x, τk, view(A, k:m, k+1:n))
    end
    qr
end


function factorize!(WS::MatrixWorkspace)
    m, n = size(WS)
    if m == n
        split_lu!(WS)
    else
        qr!(WS.qr)
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
@propagate_inbounds function _swap_rows!(B::StridedVector, i::Integer, j::Integer)
    B[i], B[j] = B[j], B[i]
    B
end

@inline function ldiv_upper!(A::AbstractMatrix, b::AbstractVector, x::AbstractVector = b)
    n = size(A, 2)
    # Eliminate two columns per pass over the right-hand side. Keep the update
    # order within each entry, including when b and x are distinct buffers.
    j = n
    @inbounds while j > 1
        xj = x[j] = (@fastmath A[j, j] \ b[j])
        b[j-1] -= A[j-1, j] * xj
        xprev = x[j-1] = (@fastmath A[j-1, j-1] \ b[j-1])
        for i = 1:j-2
            bi = b[i] - A[i, j] * xj
            b[i] = bi - A[i, j-1] * xprev
        end
        j -= 2
    end
    if j == 1
        @inbounds x[1] = (@fastmath A[1, 1] \ b[1])
    end
    b
end
@inline function ldiv_unit_lower!(
    A::AbstractMatrix,
    b::AbstractVector,
    x::AbstractVector = b,
)
    n = size(A, 2)
    j = 1
    @inbounds while j < n
        xj = x[j] = b[j]
        b[j+1] -= A[j+1, j] * xj
        xnext = x[j+1] = b[j+1]
        for i = j+2:n
            bi = b[i] - A[i, j] * xj
            b[i] = bi - A[i, j+1] * xnext
        end
        j += 2
    end
    if j == n
        @inbounds x[n] = b[n]
    end
    x
end

function lu_ldiv!(x, LU::LA.LU, b::AbstractVector)
    x === b || copyto!(x, b)
    _ipiv!(LU, x)
    ldiv_unit_lower!(LU.factors, x)
    ldiv_upper!(LU.factors, x)
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


# QR
# adaption of base implementation
function lmul_Q_adj!(A::LA.QR, b::AbstractVector)
    mA, nA = size(A.factors)
    mB = length(b)
    Afactors = A.factors
    @inbounds begin
        for k = 1:min(mA, nA)
            vBj = b[k]
            for i = k+1:mB
                vBj += conj(Afactors[i, k]) * b[i]
            end
            vBj = conj(A.τ[k]) * vBj
            b[k] -= vBj
            for i = k+1:mB
                b[i] -= Afactors[i, k] * vBj
            end
        end
    end
    b
end
function qr_ldiv!(x, QR::LA.QR, b::AbstractVector)
    # overwrites b
    # assumes QR is a tall matrix
    lmul_Q_adj!(QR, b)
    @inbounds for i = 1:length(x)
        x[i] = b[i]
    end
    ldiv_upper!(QR.factors, x)
    return x
end

function LA.ldiv!(x::AbstractVector, WS::MatrixWorkspace, b::AbstractVector)
    m, n = size(WS)
    if (m, n) == (1, 1)
        x[1] = b[1] / WS[1, 1]
        return x
    end
    WS.factorized[] || factorize!(WS)
    if m == n
        split_lu_ldiv!(x, WS, b, WS.scaled[] ? WS.row_scaling : nothing)
    else
        WS.r .= b
        qr_ldiv!(x, WS.qr, WS.r)
    end
    x
end

function LA.inv!(B::AbstractMatrix, M::MatrixWorkspace)
    lu = synced_lu(M)
    B .= lu.factors
    LA.inv!(lu)
    @inbounds for i in eachindex(B)
        lu.factors[i], B[i] = B[i], lu.factors[i]
    end
    return B
end

"""
    skeel_row_scaling!(W::MatrixWorkspace, c)
    skeel_row_scaling!(d, A, c;
    scaling_threshold = -30.0)

Compute optimal scaling factors `d` for the matrix `A` following Skeel [^S79]
if `c` is approximately of the order of the solution of the linear system
of interest.
The scaling factors are rounded to powers `e` of the base radix 2. Row scaling is only applied to rows with `e - m ≥ scaling_threshold` (this is an inequality of norms at log-scale) to prevent zero rows from being scaled, where `2^m` is approximately the maximum norm of a row of `A`. 

[^S79]: Skeel, Robert D. "Scaling for numerical stability in Gaussian elimination." Journal of the ACM (JACM) 26.3 (1979): 494-526.
"""
function skeel_row_scaling!(
    d::AbstractVector{<:Real},
    A::AbstractMatrix{<:Complex},
    c::AbstractVector{<:Real};
    scaling_threshold::Float64 = -30.0,
)
    n = length(c)
    @inbounds d .= zero(eltype(d))
    @inbounds for j = 1:n
        cj = c[j]
        for i = 1:n
            d[i] += fast_abs(A[i, j]) * cj
        end
    end

    m = maximum(d)
    s = scaling_threshold + m
    @inbounds for i = 1:n
        e = last(frexp(d[i]))
        if e < s
            d[i] = 1.0
        else
            d[i] = exp2(-e)
        end
    end

    d
end

function skeel_row_scaling!(W::MatrixWorkspace, c::AbstractVector{<:Real})
    skeel_row_scaling!(W.row_scaling, W.A, c)
    W
end

function row_scaling!(
    d::AbstractVector{<:Real},
    WS::MatrixWorkspace,
    c::AbstractVector{<:Real},
    scaling_threshold::Float64,
)
    m, n = size(WS)
    if m == n
        skeel_row_scaling!(d, WS.A, c; scaling_threshold = scaling_threshold)
    else
        d .= 1.0
    end
    d
end




"""
    apply_row_scaling!(W::MatrixWorkspace)

Apply the computed row scaling.
"""
function apply_row_scaling!(W::MatrixWorkspace)
    A, d = W.lu.factors, W.row_scaling
    m, n = size(A)
    # square matrices are scaled by `split_lu!`
    if m != n
        @inbounds for j = 1:n, i = 1:m
            A[i, j] = A[i, j] * d[i]
        end
    end
    W.scaled[] = true
    W
end

## Iterative Refinement
"""
    residual!(r, A, x, b)

Compute the residual `Ax-b` and store it in `r`
"""
function residual!(r::AbstractVector, A::AbstractMatrix, x::AbstractVector, b)
    m, n = size(A)
    @boundscheck m == length(b) == length(r) && n == length(x)
    r .= 0
    @inbounds for j = 1:n
        x_j = x[j]
        for i = 1:m
            r[i] += A[i, j] * x_j
        end
    end
    @inbounds for i = 1:m
        r[i] -= b[i]
    end
    r
end

"""
    mixed_precision_iterative_refinement!(x, A, b, norm = nothing)

Perform one step of mixed precision iterative refinement.
If `norm` is an `AbstractNorm` then the normwise relative error before
the refinement step is returned otherwise `x` is returned.
"""
function mixed_precision_iterative_refinement!(
    x::AbstractVector,
    M::MatrixWorkspace,
    b::AbstractVector,
    norm::Union{AbstractNorm,Nothing} = nothing,
)
    M.x̄ .= x
    residual!(M.r̄, M.A, M.x̄, b)
    M.r .= M.r̄
    LA.ldiv!(M.δx, M, M.r)
    x .-= M.δx
    if norm isa Nothing
        x
    else
        norm(M.δx) / norm(x)
    end
end

"""
    fixed_precision_iterative_refinement!(x, A, b, norm = nothing)

Perform one step of mixed precision iterative refinement.
If `norm` is an `AbstractNorm` then the normwise relative error before
the refinement step is returned otherwise `x` is returned.
"""
function fixed_precision_iterative_refinement!(
    x::AbstractVector,
    M::MatrixWorkspace,
    b::AbstractVector,
    norm::Union{AbstractNorm,Nothing} = nothing,
)
    residual!(M.r, M.A, x, b)
    LA.ldiv!(M.δx, M, M.r)
    x .-= M.δx
    if norm isa Nothing
        x
    else
        norm(M.δx) / norm(x)
    end
end



# CONDITION NUMBERS

"""
    inverse_inf_norm_est(WS::MatrixWorkspace,
                         d_l::Union{Nothing,Vector{<:Real}}=nothing
                         d_r::Union{Nothing,Vector{<:Real}}=nothing)

Estimate of the infinity norm of `diag(d_r)⁻¹A⁻¹diag(d_l)⁻¹` where `d_l` and `d_r` are
optional positive vectors.
If `d_l` or `d_r` is `nothing` the all one vector is used.
This uses the 1-norm lapack condition estimator described by Highahm in [^H88].

[^H88]: Higham, Nicholas J. "FORTRAN codes for estimating the one-norm of a real or complex matrix, with applications to condition estimation." ACM Transactions on Mathematical Software (TOMS) 14.4 (1988): 381-396.
"""
function inverse_inf_norm_est(
    WS::MatrixWorkspace,
    d_l::Union{Nothing,Vector{<:Real}} = nothing,
    d_r::Union{Nothing,Vector{<:Real}} = nothing,
)
    lu = synced_lu(WS)
    work, rwork = WS.inf_norm_est_work, WS.inf_norm_est_rwork
    if !WS.scaled[]
        inverse_inf_norm_est(lu, d_l, d_r, nothing, work, rwork)
    else
        inverse_inf_norm_est(lu, d_l, d_r, WS.row_scaling, work, rwork)
    end
end
function inverse_inf_norm_est(
    lu::LA.LU,
    d_l::Union{Nothing,Vector{<:Real}},
    d_r::Union{Nothing,Vector{<:Real}},
    row_scaling::Union{Nothing,Vector{<:Real}},
    work::Vector{<:Complex},
    rwork::Vector{<:Real},
)
    z = ξ = y = work
    x = rwork

    n = size(lu.factors, 1)
    x .= inv(n)
    if d_r !== nothing
        x ./= d_r
    end
    lu_ldiv_adj!(y, lu, x)
    if d_l !== nothing
        y ./= d_l
    end
    if row_scaling !== nothing
        y .*= row_scaling
    end
    γ = sum(fast_abs, y)
    ξ ./= fast_abs.(y)
    if d_l !== nothing
        ξ ./= d_l
    end
    if row_scaling !== nothing
        ξ ./= row_scaling
    end
    lu_ldiv!(z, lu, ξ)

    if d_r !== nothing
        x .= real.(z) ./ d_r
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
        if d_r !== nothing
            x ./= d_r
        end
        lu_ldiv_adj!(y, lu, x)
        if d_l !== nothing
            y ./= d_l
        end
        if row_scaling !== nothing
            y .*= row_scaling
        end
        γ̄ = γ
        γ = sum(fast_abs, y)
        if γ ≤ γ̄
            γ = γ̄
            break
        end
        ξ .= y ./ fast_abs.(y)
        if d_l !== nothing
            ξ ./= d_l
        end
        if row_scaling !== nothing
            ξ ./= row_scaling
        end
        lu_ldiv!(z, lu, ξ)
        if d_r !== nothing
            x .= real.(z) ./ d_r
        else
            x .= real.(z)
        end
        k += 1
        if x[j] == LA.norm(x, Inf) || k > 2
            break
        end
    end
    nanmin(γ, Inf)
end

function inf_norm(
    WS::MatrixWorkspace,
    d_l::Union{Nothing,Vector{<:Real}} = nothing,
    d_r::Union{Nothing,Vector{<:Real}} = nothing,
)
    norm = -Inf
    A = WS.A
    n, m = size(A)
    @inbounds for i = 1:n
        normᵢ = 0.0
        for j = 1:m
            if d_r isa Nothing
                normᵢ += fast_abs(A[i, j])
            else
                normᵢ += fast_abs(A[i, j]) * d_r[j]
            end
        end
        if !(d_l isa Nothing)
            normᵢ *= d_l[i]
        end
        norm = @fastmath max(norm, normᵢ)
    end
    norm
end

function max_min_row(
    WS::MatrixWorkspace,
    d_l::Union{Nothing,Vector{<:Real}} = nothing,
    d_r::Union{Nothing,Vector{<:Real}} = nothing,
)
    max_row = -Inf
    min_row = Inf
    A = WS.A
    n, m = size(A)
    @inbounds for i = 1:n
        normᵢ = 0.0
        for j = 1:m
            if d_r isa Nothing
                normᵢ += fast_abs(A[i, j])
            else
                normᵢ += fast_abs(A[i, j]) * d_r[j]
            end
        end
        if !(d_l isa Nothing)
            normᵢ *= d_l[i]
        end
        max_row = @fastmath max(max_row, normᵢ)
        min_row = @fastmath min(min_row, normᵢ)
    end
    max_row, min_row
end

"""
    cond(A::MatrixWorkspace,
         d_l::Union{Nothing,Vector{<:Real}} = nothing,
         d_r::Union{Nothing,Vector{<:Real}} = nothing)

Compute the condition number w.r.t. the infinity-norm of `diag(d_l) * A * diag(d_r)`.
If `d_l` or `d_r` is `nothing` the all one vector is used.
If `size(A) == (1,1)` then just the norm of the inverse is returned.
"""
function LA.cond(
    WS::MatrixWorkspace,
    d_l::Union{Nothing,Vector{<:Real}} = nothing,
    d_r::Union{Nothing,Vector{<:Real}} = nothing,
)
    m, n = size(WS)
    if m == n == 1
        if isa(d_l, Nothing) && isa(d_r, Nothing)
            inv(abs(WS.A[1, 1]))
        elseif isa(d_l, Nothing)
            inv(abs(WS.A[1, 1]) * d_r[1])
        elseif isa(d_r, Nothing)
            inv(d_l[1] * abs(WS.A[1, 1]))
        else
            inv(d_l[1] * abs(WS.A[1, 1]) * d_r[1])
        end
    elseif m > n
        WS.factorized[] || factorize!(WS)
        rmax, rmin = -Inf, Inf
        for i = 1:n
            rᵢ = fast_abs(WS.qr.factors[i, i]) * d_r[i]
            rmax = max(rmax, rᵢ)
            rmin = min(rmin, rᵢ)
        end
        #@show (rmax, rmin)
        rmax / rmin
    else
        inverse_inf_norm_est(WS, d_l, d_r) * inf_norm(WS, d_l, d_r)
    end
end

function egcond(
    WS::MatrixWorkspace,
    d_l::Union{Nothing,Vector{<:Real}} = nothing,
    d_r::Union{Nothing,Vector{<:Real}} = nothing,
)
    m, n = size(WS)
    if m == n == 1
        a = abs(WS.A[1, 1])
        if d_l !== nothing
            a *= d_l[1]
        end
        if d_r !== nothing
            a *= d_r[1]
        end
        max(a, 1.0) / a
    elseif m > n
        WS.factorized[] || factorize!(WS)
        rmax, rmin = -Inf, Inf
        for i = 1:n
            rᵢ = abs(WS.qr.factors[i, i]) * d_r[i]
            rmax = max(rmax, rᵢ)
            rmin = min(rmin, rᵢ)
        end
        rmax / rmin
    else
        max_row, min_row = max_min_row(WS, d_l, d_r)
        inverse_inf_norm_est(WS, d_l, d_r) * max(max_row, 1 / min_row)
    end
end

#####################
## Jacobian ##
#####################
struct Jacobian{M}
    workspace::MatrixWorkspace{M}
    # stats
    factorizations::Base.RefValue{Int}
    ldivs::Base.RefValue{Int}
end
Jacobian(A::AbstractMatrix) = Jacobian(MatrixWorkspace(A), Ref(0), Ref(0))

updated!(J::Jacobian) = (updated!(J.workspace); J)

workspace(J::Jacobian) = J.workspace
matrix(J::Jacobian) = matrix(workspace(J))

function Base.show(io::IO, J::Jacobian{T}) where {T}
    println(io, "Jacobian{$T}:")
    println(io, " • # factorizations → ", J.factorizations[])
    println(io, " • # ldivs → ", J.ldivs[])
end

"""
    ldiv!(x̂::AbstractVector, J::Jacobian, b::AbstractVector)

solve the linear system `matrix(J)x̂=b`.
"""
function LA.ldiv!(
    x̂::AbstractVector,
    J::Jacobian,
    b::AbstractVector,
    norm = nothing;
    row_scaling::Bool = true,
)
    # stats update
    J.factorizations[] += !workspace(J).factorized[]
    J.ldivs[] += 1
    LA.ldiv!(x̂, workspace(J), b)
    x̂
end
function LA.ldiv!(
    x̂::AbstractVector,
    J::Jacobian,
    b::AbstractVector,
    norm::WeightedNorm;
    row_scaling::Bool = true,
)
    if !workspace(J).factorized[] && row_scaling
        skeel_row_scaling!(workspace(J), weights(norm))
        apply_row_scaling!(workspace(J))
    end
    # stats update
    J.factorizations[] += !workspace(J).factorized[]
    J.ldivs[] += 1
    LA.ldiv!(x̂, workspace(J), b)
    x̂
end

function iterative_refinement!(
    x::AbstractVector,
    J::Jacobian,
    b::AbstractVector,
    norm::AbstractNorm;
    max_iters::Int = 3,
    tol::Float64 = sqrt(eps()),
)
    J.ldivs[] += 1
    δ = mixed_precision_iterative_refinement!(x, workspace(J), b, norm)
    δ < tol && return (accuracy = δ, diverged = false)
    for i = 2:max_iters
        J.ldivs[] += 1
        δ′ = mixed_precision_iterative_refinement!(x, workspace(J), b, norm)
        if δ′ < tol
            return (accuracy = δ′, diverged = false)
        elseif δ′ > 0.5 * δ
            return (accuracy = δ′, diverged = true)
        end
        δ = δ′
    end
    (accuracy = δ, diverged = false)
end


"""
    init!(J::Jacobian)

(Re-)initialize the `Jacobian`.
"""
function init!(J::Jacobian; keep_stats::Bool = false)
    J.factorizations[] = 0
    J.ldivs[] = 0
    J
end

function LA.cond(
    J::Jacobian,
    d_l::Union{Nothing,Vector{<:Real}} = nothing,
    d_r::Union{Nothing,Vector{<:Real}} = nothing,
)
    LA.cond(workspace(J), d_l, d_r)
end

function egcond(
    J::Jacobian,
    d_l::Union{Nothing,Vector{<:Real}} = nothing,
    d_r::Union{Nothing,Vector{<:Real}} = nothing,
)
    egcond(workspace(J), d_l, d_r)
end

@static if VERSION ≥ v"1.7.0-"
    qr_col_norm!(M) = LA.qr!(M, LA.ColumnNorm())
else
    qr_col_norm!(M) = LA.qr!(M, Val(true))
end
