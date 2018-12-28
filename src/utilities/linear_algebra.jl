export normalized_dot, euclidean_distance, euclidean_norm, rowscaling, rowscaling!, solve!, fast_factorization!, fast_ldiv!, factorization, factorize!

"""
    normalized_dot(u, v)

Compute u⋅v / (||u||*||v||).
"""
function normalized_dot(u::AbstractVector{T}, v::AbstractVector{T}) where T
    @boundscheck length(u) == length(v)
    n = length(u)
    dot = zero(eltype(u))
    norm_u = norm_v = zero(real(eltype(u)))
    for i=1:n
        dot += conj(v[i]) * u[i]
        norm_u += abs2(u[i])
        norm_v += abs2(v[i])
    end

    @fastmath dot / sqrt(norm_u * norm_v)
end

"""
    euclidean_distance(u, v)

Compute ||u-v||₂.
"""
function euclidean_distance(x::AbstractVector{T}, y::AbstractVector{T}) where T
    @boundscheck length(x) == length(y)
    n = length(x)
    @inbounds d = abs2(x[1] - y[1])
    @inbounds for i=2:n
        @fastmath d += abs2(x[i] - y[i])
    end
    sqrt(d)
end

function euclidean_norm(x::AbstractVector)
    out = zero(real(eltype(x)))
    @inbounds for i in eachindex(x)
        out += abs2(x[i])
    end
    sqrt(out)
end

"""
    rowscaling(A)

Apply elementwise row-scaling to a copy of `A`.
"""
rowscaling(A) = rowscaling!(copy(A))

"""
    rowscaling(A, b=nothing)

Apply elementwise row-scaling inplace to `A`.
Also scale the i-th entry of `b` by the same value as the i-th row of `A`.
"""
function rowscaling!(A::AbstractMatrix, b::Union{AbstractVector, Nothing}=nothing)
    @inbounds for i=1:size(A, 1)
        s = zero(eltype(A))
        for j=1:size(A, 2)
            s += abs2(A[i, j])
        end
        s = sqrt(s)
        for j=1:size(A, 2)
            A[i, j] /= s
        end
        if b !== nothing
            b[i] /= s
        end
    end
    A
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
        copyto!(x, b)
    end
    apply_rowscaling!(x, Jac)
    solve!(x, Jac.fac, x)
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
function factorize!(fact::LinearAlgebra.Factorization, A::AbstractMatrix)
    LinearAlgebra.qr!(A)
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
                           b::Union{AbstractVector{T}, Nothing}=nothing,
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
                    absi = abs2(A[i,k])
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




# # Iterative refinement
# function condition(A, x, b)
#     cond = 0.0
#     for k in 1:maxiters
#         residual!(r, A, x, b)
#         norm_r = euclidean_norm(r)
#         if norm_r < tol
#             break
#         end
#     end
# end
#

"""
    residual!(u::AbstractVector{T}, A, x, b, T_res)

Compute the residual `Ax-b` in precision `T_res` and store in precision `T` in `u`.
"""
function residual!(u::AbstractVector{Complex{T}}, A, x, b, ::Type{T_res}) where {T, T_res}
    @boundscheck size(A, 1) == length(b) && size(A,2) == length(x)
    S = Complex{T_res}
    m, n = size(A)
    @inbounds for i in 1:m
        dot = zero(S)
        for j in 1:n
            dot = multiply_add(A[i,j], x[j], dot)
        end
        u[i] = sub(dot, b[i])
    end
    u
end
