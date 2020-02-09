struct MatrixWorkspace{M<:AbstractMatrix{ComplexF64}} <: AbstractMatrix{ComplexF64}
    A::M # Matrix
    d::Vector{Float64} # Inverse of scaling factors
    factorized::Base.RefValue{Bool}
    lu::LA.LU{ComplexF64,M} # LU Factorization of D * J
    row_scaling::Vector{Float64}
    scaled::Base.RefValue{Bool}
    # mixed precision iterative refinement
    x̄::Vector{ComplexDF64}
    r::Vector{ComplexF64}
    r̄::Vector{ComplexDF64}
    δx::Vector{ComplexF64}
end

function MatrixWorkspace(Â::AbstractMatrix)
    m, n = size(Â)
    m == n || throw(ArgumentError("Input is not a square matrix"))

    A = Matrix{ComplexF64}(Â)
    d = ones(m)
    factorized = Ref(false)
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

    MatrixWorkspace(A, d, factorized, lu, row_scaling, scaled, x̄, r, r̄, δx)
end

Base.size(MW::MatrixWorkspace) = size(MW.A)

import Base: @propagate_inbounds
@propagate_inbounds Base.getindex(MW::MatrixWorkspace, i::Integer) =
    getindex(MW.A, i)
@propagate_inbounds Base.getindex(MW::MatrixWorkspace, i::Integer, j::Integer) =
    getindex(MW.A, i, j)
@propagate_inbounds Base.setindex!(MW::MatrixWorkspace, x, i::Integer) =
    setindex!(MW.A, x, i)
@propagate_inbounds Base.setindex!(
    MW::MatrixWorkspace,
    x,
    i::Integer,
    j::Integer,
) = setindex!(MW.A, x, i, j)
@propagate_inbounds Base.copyto!(MW::MatrixWorkspace, A::AbstractArray) =
    copyto!(MW.A, A)

"""
    updated!(MW::MatrixWorkspace)

Indicate that the matrix `MW` got updated.
"""
function updated!(MW::MatrixWorkspace)
    MW.factorized[] = false
    MW.scaled[] = false
    @inbounds copyto!(MW.lu.factors, MW.A)
    MW
end

"""
    update!(MW::MatrixWorkspace, A::Matrix)

Update the matrix in `MW` with `A`.
"""
@inline function update!(MW::MatrixWorkspace, A::Matrix)
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
    b::Union{AbstractVector,Nothing} = nothing,
    ipiv::Union{Vector{I},Nothing} = nothing,
) where {T,I<:Integer,Pivot}
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
                    if !(b isa Nothing)
                        b[k], b[kp] = b[kp], b[k]
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

function factorize!(WS::MatrixWorkspace)
    lu!(WS.lu.factors, nothing, WS.lu.ipiv)
    WS.factorized[] = true
    WS
end

##########
## ldiv ##
##########
@inline _ipiv!(A::LA.LU, b::AbstractVector) =
    apply_ipiv!(b, 1:length(A.ipiv), A.ipiv)
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
@propagate_inbounds function _swap_rows!(
    B::StridedVector,
    i::Integer,
    j::Integer,
)
    B[i], B[j] = B[j], B[i]
    B
end

@inline function ldiv_upper!(
    A::AbstractMatrix,
    b::AbstractVector,
    x::AbstractVector = b,
)
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

function ldiv_adj_unit_lower!(
    A::AbstractMatrix,
    b::AbstractVector,
    x::AbstractVector = b,
)
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


function LA.ldiv!(x::AbstractVector, WS::MatrixWorkspace, b::AbstractVector)
    WS.factorized[] || factorize!(WS)
    if WS.scaled[]
        x .= WS.row_scaling .* b
        lu_ldiv!(x, WS.lu, x)
    else
        lu_ldiv!(x, WS.lu, b)
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
    @inbounds for j = 1:m
        x_j = x[j]
        for i = 1:n
            r[i] += A[i, j] * x_j
        end
    end
    r .-= b
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

#####################
## Jacobian ##
#####################
struct Jacobian{M}
    J::MatrixWorkspace{M}
    # stats
    factorizations::Base.RefValue{Int}
    ldivs::Base.RefValue{Int}
end
Jacobian(A::AbstractMatrix) = Jacobian(MatrixWorkspace(A), Ref(0), Ref(0))

updated!(JM::Jacobian) = (updated!(JM.J); JM)
jacobian(JM::Jacobian) = JM.J

function Base.show(io::IO, JM::Jacobian{T}) where {T}
    println(io, "Jacobian{$T}:")
    println(io, " • # factorizations → ", JM.factorizations[])
    println(io, " • # ldivs → ", JM.ldivs[])
end

"""
    ldiv!(x̂::AbstractVector, JM::Jacobian, b::AbstractVector)

solve the linear system `jacobian(JM)x̂=b`.
"""
function LA.ldiv!(
    x̂::AbstractVector,
    JM::Jacobian,
    b::AbstractVector,
    norm = nothing;
    row_scaling::Bool = true,
)
    # stats update
    JM.factorizations[] += !jacobian(JM).factorized[]
    JM.ldivs[] += 1
    LA.ldiv!(x̂, jacobian(JM), b)
    x̂
end
function LA.ldiv!(
    x̂::AbstractVector,
    JM::Jacobian,
    b::AbstractVector,
    norm::WeightedNorm;
    row_scaling::Bool = true,
)
    if !jacobian(JM).factorized[] && row_scaling
        skeel_row_scaling!(jacobian(JM), weights(norm))
        apply_row_scaling!(jacobian(JM))
    end
    # stats update
    JM.factorizations[] += !jacobian(JM).factorized[]
    JM.ldivs[] += 1
    LA.ldiv!(x̂, jacobian(JM), b)
    x̂
end

function iterative_refinement!(
    x::AbstractVector,
    JM::Jacobian,
    b::AbstractVector,
    norm = nothing;
    row_scaling::Bool = true,
    fixed_precision::Bool = false
)
    JM.ldivs[] += 1
    if fixed_precision
        fixed_precision_iterative_refinement!(x, jacobian(JM), b, norm)
    else
        mixed_precision_iterative_refinement!(x, jacobian(JM), b, norm)
    end
end


"""
    init!(JM::Jacobian)

(Re-)initialize the `Jacobian`.
"""
function init!(JM::Jacobian; keep_stats::Bool = false)
    JM.factorizations[] = 0
    JM.ldivs[] = 0
    JM
end
