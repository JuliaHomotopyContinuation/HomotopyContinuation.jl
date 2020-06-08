"""
    MatrixWorkspace

This is a data structure for the efficient repeated solution of a square or
overdetermined linear system `Ax=b`.
"""
struct MatrixWorkspace{M<:AbstractMatrix{ComplexF64}} <: AbstractMatrix{ComplexF64}
    A::M # Matrix
    d::Vector{Float64} # Inverse of scaling factors
    factorized::Base.RefValue{Bool}
    lu::LA.LU{ComplexF64,M} # LU Factorization of D * J
    qr::LA.QR{ComplexF64,Matrix{ComplexF64}}
    row_scaling::Vector{Float64}
    scaled::Base.RefValue{Bool}
    # mixed precision iterative refinement
    x̄::Vector{ComplexDF64}
    r::Vector{ComplexF64}
    r̄::Vector{ComplexDF64}
    δx::Vector{ComplexF64}
    inf_norm_est_work::Vector{ComplexF64}
    inf_norm_est_rwork::Vector{Float64}
end

MatrixWorkspace(m::Integer, n::Integer) = MatrixWorkspace(zeros(ComplexF64, m, n))
function MatrixWorkspace(Â::AbstractMatrix)
    m, n = size(Â)
    m ≥ n || throw(ArgumentError("Expected system with more rows than columns."))

    A = Matrix{ComplexF64}(Â)
    d = ones(m)
    factorized = Ref(false)
    qr = LA.qrfactUnblocked!(copy(A))
    # experiments show that for m > 25 the data layout as a
    # struct array is beneficial
    if m > 25
        A = StructArrays.StructArray(A)
    end
    row_scaling = ones(m)
    scaled = Ref(false)

    lu = LinearAlgebra.LU{eltype(A),typeof(A)}(copy(A), zeros(Int, m), 0)
    r = zeros(ComplexF64, m)
    r̄ = zeros(ComplexDF64, m)
    x̄ = zeros(ComplexDF64, n)
    δx = zeros(ComplexF64, n)
    inf_norm_est_work = Vector{ComplexF64}(undef, n)
    inf_norm_est_rwork = Vector{Float64}(undef, n)

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
    if m == n
        @inbounds copyto!(MW.lu.factors, MW.A)
    else
        @inbounds copyto!(MW.qr.factors, MW.A)
    end
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
        lu!(WS.lu.factors, WS.lu.ipiv)
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
    @inbounds for j = n:-1:1
        # singular_exception && iszero(A[j, j]) && throw(LA.SingularException(j))
        xj = x[j] = (@fastmath A[j, j] \ b[j])
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
function qr_ldiv!(x, QR::LA.QR, b::AbstractVector) where {T}
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
        if WS.scaled[]
            x .= WS.row_scaling .* b
            lu_ldiv!(x, WS.lu, x)
        else
            lu_ldiv!(x, WS.lu, b)
        end
    else
        WS.r .= b
        qr_ldiv!(x, WS.qr, WS.r)
    end
    x
end

"""
    skeel_row_scaling!(W::MatrixWorkspace, c)
    skeel_row_scaling!(d, A, c)

Compute optimal scaling factors `d` for the matrix `A` following Skeel [^S79]
if `c` is approximately of the order of the solution of the linear system
of interest.
The scaling factors are rounded to powers of the base radix 2.

[^S79]: Skeel, Robert D. "Scaling for numerical stability in Gaussian elimination." Journal of the ACM (JACM) 26.3 (1979): 494-526.
"""
function skeel_row_scaling!(
    d::AbstractVector{<:Real},
    A::AbstractMatrix{<:Complex},
    c::AbstractVector{<:Real},
)
    n = length(c)
    @inbounds d .= zero(eltype(d))
    @inbounds for j = 1:n
        cj = c[j]
        for i = 1:n
            d[i] += fast_abs(A[i, j]) * cj
        end
    end
    @inbounds for i = 1:n
        d[i] = exp2(-last(frexp(d[i])))
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
)
    m, n = size(WS)
    if m == n
        skeel_row_scaling!(d, WS.A, d)
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
    @inbounds for j = 1:n, i = 1:m
        A[i, j] = A[i, j] * d[i]
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
    WS.factorized[] || factorize!(WS)
    work, rwork = WS.inf_norm_est_work, WS.inf_norm_est_rwork
    if !WS.scaled[]
        inverse_inf_norm_est(WS.lu, d_l, d_r, nothing, work, rwork)
    else
        inverse_inf_norm_est(WS.lu, d_l, d_r, WS.row_scaling, work, rwork)
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
        γ ≤ γ̄ && break
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
    γ
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
jacobian(J::Jacobian) = J.workspace

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
