export normalized_dot, euclidean_distance, euclidean_norm, rowscaling, rowscaling!, solve!,
    fast_factorization!, fast_ldiv!, factorization, factorize!,
    iterative_refinement!, iterative_refinement_step!, IterativeRefinementResult,
    Jacobian


"""
    Jacobian{T, <:Factorization}

Data structure holding a Jacobian matrix `J` with it's factorization.
Additional contains a vector `D` of scalings which row equilibarate `J`.

    Jacobian(J; rowscaling=false)

Setup a Jacobian with Jacobian matrices of the form `J`. Enable row scaling
if desired.
"""
mutable struct Jacobian{T, F<:LinearAlgebra.Factorization}
    J::Matrix{T} # Jacobian matrix
    D::Vector{Float64} # Scaling factors D
    fac::F # Factorization of D * J
    rowscaling::Bool
    r::Vector{T} # Vector for iterative refinement
end

function Jacobian(A::AbstractMatrix; rowscaling=false)
    m, n = size(A)
    if m == n
        fac = LinearAlgebra.lu(A)
    else
        fac = LinearAlgebra.qr(A)
    end
    J = copy(A)
    D = ones(m)
    r = similar(J, size(J, 1))
    Jacobian(J, D, fac, rowscaling, r)
end

"""
    update!(jacobian, J)

Update the `jacobian` struct with the new Matrix `J`.
This applies rowscaling, if enabled, and computes a factorization of `J`.
"""
function update!(Jac::Jacobian, J::AbstractMatrix)
    copyto!(Jac.J, J)
    Jac.rowscaling && rowscaling!(Jac.D, Jac.J)
    Jac.fac = factorize!(Jac.fac, Jac.J)
    Jac
end

"""
    update_jacobian!(jacobian)

This applies rowscaling, if enabled, and computes a factorization of the stored Jacobian matrix.
Call this instead of `update!` if `jacobian.J` got updated.
"""
function updated_jacobian!(Jac::Jacobian)
    Jac.rowscaling && rowscaling!(Jac.D, Jac.J)
    Jac.fac = factorize!(Jac.fac, Jac.J)
    Jac
end

"""
    apply_rowscaling!(x, jacobian)

Apply the stored row scaling to the vector `x`.
"""
function apply_rowscaling!(x::AbstractVector, Jac::Jacobian)
    if Jac.rowscaling
        @boundscheck length(x) == size(Jac.D)
        @inbounds for i in eachindex(x)
            x[i] *= Jac.D[i]
        end
    end
    x
end

"""
    enable_rowscaling!(jacobian)

Enable row scaling.
"""
enable_rowscaling!(Jac::Jacobian) = begin Jac.rowscaling = true; Jac end

"""
    disable_rowscaling!(jacobian)

Disable row scaling.
"""
disable_rowscaling!(Jac::Jacobian) = begin Jac.rowscaling = true; Jac end


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
    rowscaling!(D::Vector, A)

Apply elementwise row-scaling inplace to `A`.
Stores the scaling factors in D.
"""
function rowscaling!(D::Vector, A::AbstractMatrix)
    m, n = size(A)
    @inbounds for i=1:m
        s = zero(eltype(D))
        for j=1:n
            s += abs2(A[i, j])
        end
        dᵢ = inv(sqrt(s))
        for j=1:n
            A[i, j] *= dᵢ
        end
        D[i] = dᵢ
    end
    D
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


# Iterative refinement
"""
    IterativeRefinementResult

The result of an application of [`iterative_refinement!`](@ref).

## Fields
* `iters::Int`
* `accuracy::Float64` The relative accuarcy of the last update, i.e., `||δx|| / ||x||`.
* `cond::Float64` An estimate of the condition number of the matrix, computed by `||δx|| / eps(||x||)`.
"""
struct IterativeRefinementResult
    iters::Int
    accuracy::Float64
    cond::Float64
end

"""
    solve_with_iterative_refinement!(x, jacobian::Jacobian, b, ::Type{T}, iters=1)

Apply iterative refinement on the solution `x` of the equation `Ax=b`
where `A` is the matrix stored in `jacobian`. `T` is the accuracy with which
the residual is computed.
Returns an [`IterativeRefinementResult`](@ref).
"""
function solve_with_iterative_refinement!(x::AbstractVector, Jac::Jacobian,
                b::AbstractVector, ::Type{T}; iters::Int=1) where {T}
    # we need to apply row scaling sepearetly since we need it the b scaled
    apply_rowscaling!(b, Jac)
    solve!(x, Jac.fac, b)
    cond = 0.0
    accuracy = Inf
    norm_x = maximum(abs, x)
    for iter in 1:iters
        accuracy = iterative_refinement_step!(x, Jac.J, b, Jac.fac, T_res, Jac.r)
        if iters == 1
            cond = norm_δx / (eps() * norm_x)
        end
    end
    IterativeRefinementResult(maxiters, accuracy, cond)
end

"""
    iterative_refinement_step!(x, A, b, fac, T, r)

Apply one step of iterative refinement where the residual is computed with precision `T`.
"""
function iterative_refinement_step!(x, A, b, fac, T, r)
    residual!(r, A, x, b, T)
    δx = solve!(fac, b)
    accuracy = maximum(abs, δx)
    for i in eachindex(r)
        x[i] = Complex{T}(x[i]) - Complex{T}(δx[i])
    end

    return accuracy
end

"""
    residual!(u, A, x, b, [::Type{T}])

Compute the residual `Ax-b` in precision `T` and store in `u`.
"""
residual!(u::AbstractVector, A, x, b) = residual!(u, A, x, b)
function residual!(u::AbstractVector, A, x, b, ::Type{T}) where {T}
    @boundscheck size(A, 1) == length(b) && size(A,2) == length(x)
    m, n = size(A)
    @inbounds for i in 1:m
        dot = zero(Complex{T})
        for j in 1:n
            dot = muladd(Complex{T}(A[i,j]), Complex{T}(x[j]), dot)
        end
        u[i] = dot - Complex{T}(b[i])
    end
    u
end
