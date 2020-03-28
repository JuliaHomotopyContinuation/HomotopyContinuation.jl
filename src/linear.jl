export AffineSubspace,
    AffineExtrinsic,
    AffineIntrinsic,
    intrinsic,
    extrinsic,
    dim,
    codim,
    rand_affine_subspace,
    coord_change,
    Intrinsic,
    Extrinsic,
    IntrinsicStiefel,
    ExtrinsicStiefel

"""
    Coordinates

A type used for encoding the used coordinates and performing coordinate changes.
Currently supported coordinates are:
- [`Intrinsic`](@ref)
- [`Extrinsic`](@ref)
- [`IntrinsicStiefel`](@ref)
- [`ExtrinsicStiefel`](@ref)
"""
struct Coordinates{T} end
Base.broadcastable(C::Coordinates) = Ref(C)

# TODO: DOCUMENTATION
const Intrinsic = Coordinates{:Intrinsic}()
const Extrinsic = Coordinates{:Extrinsic}()
const IntrinsicStiefel = Coordinates{:IntrinsicStiefel}()
const ExtrinsicStiefel = Coordinates{:ExtrinsicStiefel}()

"""
    rand_unitary_matrix(n::Int, T=ComplexF64)

Samples a `n × n` unitary Matrix uniformly from the space of all unitary n × n matrices.

See https://arxiv.org/abs/math-ph/0609050 for a derivation.
"""
function rand_unitary_matrix(n::Int, T::Type = ComplexF64)
    Z = randn(T, n, n) ./ sqrt(2)
    Q, R = LA.qr(Z, Val(true))
    Λ = LA.diagm(0 => [R[i, i] / abs(R[i, i]) for i = 1:n])
    Q * Λ
end

struct AffineExtrinsic{T}
    A::Matrix{T}
    b::Vector{T}
end

(A::AffineExtrinsic)(x::AbstractVector) = A.A * x - A.b

dim(A::AffineExtrinsic) = size(A.A, 2) - size(A.A, 1)
codim(A::AffineExtrinsic) = size(A.A, 1)

function Base.copy!(A::AffineExtrinsic, B::AffineExtrinsic)
    copy!(A.A, B.A)
    copy!(A.b, B.b)
    A
end

struct AffineIntrinsic{T}
    # orthogonal coordinates
    A::Matrix{T}
    b₀::Vector{T}
    # stiefel coordinates
    Y::Matrix{T}
end
AffineIntrinsic(A::Matrix, b₀::Vector) = AffineIntrinsic(A, b₀, stiefel_coordinates(A, b₀))

(A::AffineIntrinsic)(u::AbstractVector, ::Coordinates{:Intrinsic}) = A.A * u + A.b₀
(A::AffineIntrinsic)(u::AbstractVector, ::Coordinates{:IntrinsicStiefel}) = A.Y * u

function stiefel_coordinates(A::AbstractMatrix, b::AbstractVector)
    γ = sqrt(1 + sum(abs2, b))
    n, k = size(A)
    Y = zeros(eltype(A), n + 1, k + 1)
    Y[1:n, 1:k] .= A
    Y[1:n, k+1] .= b ./ γ
    Y[n+1, k+1] = 1 / γ
    Y
end

dim(I::AffineIntrinsic) = size(I.A, 2)
codim(I::AffineIntrinsic) = size(I.A, 1) - size(I.A, 2)

function Base.copy!(A::AffineIntrinsic, B::AffineIntrinsic)
    copy!(A.A, B.A)
    copy!(A.b₀, B.b₀)
    copy!(A.Y, B.Y)
    A
end

function AffineIntrinsic(E::AffineExtrinsic)
    svd = LA.svd(E.A; full = true)
    b₀ = svd \ E.b
    m, n = size(E.A)
    A = Matrix((@view svd.Vt[m+1:end, :])')
    AffineIntrinsic(A, b₀)
end

function AffineExtrinsic(I::AffineIntrinsic)
    svd = LA.svd(I.A; full = true)
    m, n = size(I.A)
    A = copy((@view svd.U[:, n+1:end])')
    b = A * I.b₀
    AffineExtrinsic(A, b)
end


struct AffineSubspace{T}
    extrinsic::AffineExtrinsic{T}
    intrinsic::AffineIntrinsic{T}
end

AffineSubspace(I::AffineIntrinsic) = AffineSubspace(AffineExtrinsic(I), I)
AffineSubspace(E::AffineExtrinsic) = AffineSubspace(E, AffineIntrinsic(E))

function AffineSubspace(A::AbstractMatrix{T}, b::AbstractVector{T}) where {T}
    size(A, 1) == length(b) || throw(ArgumentError("Size of A and b not compatible."))
    0 < size(A, 1) < size(A, 2) ||
    throw(ArgumentError("Affine subspace has to be given in extrinsic coordinates, i.e., by A x = b."))

    AffineSubspace(AffineExtrinsic(Matrix(A), Vector(b)))
end

Base.broadcastable(A::AffineSubspace) = Ref(A)

dim(A::AffineSubspace) = dim(A.intrinsic)
codim(A::AffineSubspace) = codim(A.intrinsic)

intrinsic(A::AffineSubspace) = A.intrinsic
extrinsic(A::AffineSubspace) = A.extrinsic

function Base.copy!(A::AffineSubspace, B::AffineSubspace)
    copy!(A.intrinsic, B.intrinsic)
    copy!(A.extrinsic, B.extrinsic)
    A
end

function (A::AffineSubspace)(x::AbstractVector, ::Coordinates{:Intrinsic})
    intrinsic(A)(x, Intrinsic)
end
function (A::AffineSubspace)(x::AbstractVector, ::Coordinates{:IntrinsicStiefel})
    intrinsic(A)(x, IntrinsicStiefel)
end
function (A::AffineSubspace)(x::AbstractVector, ::Coordinates{:Extrinsic})
    extrinsic(A)(x)
end

function rand_affine_subspace(
    n::Integer;
    dim::Union{Nothing,Integer} = nothing,
    codim::Union{Nothing,Integer} = nothing,
    real::Bool = false,
)
    !isnothing(dim) ||
    !isnothing(codim) ||
    throw(ArgumentError("Neither `dim` or `codim` specified."))

    if !isnothing(dim)
        0 < dim < n || throw(ArgumentError("`dim` has to be between 0 and `n`."))
        k = dim
    else
        0 < codim < n || throw(ArgumentError("`codim` has to be between 0 and `n`."))
        k = n - codim
    end
    Q = rand_unitary_matrix(n, real ? Float64 : ComplexF64)

    # ind = Random.shuffle(1:n)
    A = Q[:, 1:k]
    # This is also an orthonormal vector, although it only has to be orthogonal
    b = Q[:, k+1]
    Y = stiefel_coordinates(A, b)
    AffineSubspace(AffineIntrinsic(A, b, Y))
end


function dist(A::AffineIntrinsic, B::AffineIntrinsic)
    sqrt(sum(σᵢ -> acos(min(σᵢ, 1.0))^2, LA.svdvals(A.Y' * B.Y)))
end
dist(A::AffineSubspace, B::AffineSubspace) = dist(A.intrinsic, B.intrinsic)

# Coordinate changes
coord_change(A::AffineSubspace, ::C, ::C, x) where {C<:Coordinates} = x
coord_change(A::AffineSubspace, ::Coordinates{:Intrinsic}, ::Coordinates{:Extrinsic}, u) =
    A(u, Intrinsic)
function coord_change(
    A::AffineSubspace,
    ::Coordinates{:Intrinsic},
    ::Coordinates{:IntrinsicStiefel},
    u,
)
    [u; 1 / A.intrinsic.Y[end, end]]
end
function coord_change(
    A::AffineSubspace,
    ::Coordinates{:Intrinsic},
    ::Coordinates{:ExtrinsicStiefel},
    u,
)
    v = coord_change(A, Intrinsic(), IntrinsicStiefel(), u)
    coord_change(A, IntrinsicStiefel(), ExtrinsicStiefel(), v)
end

coord_change(A::AffineSubspace, ::Coordinates{:Extrinsic}, ::Coordinates{:Intrinsic}, x) =
    A.intrinsic.A' * (x - A.intrinsic.b₀)
function coord_change(
    A::AffineSubspace,
    ::Coordinates{:Extrinsic},
    ::Coordinates{:IntrinsicStiefel},
    x,
)
    u = coord_change(A, Extrinsic(), Intrinsic(), x)
    coord_change(A, Intrinsic(), IntrinsicStiefel(), u)
end
function coord_change(
    A::AffineSubspace,
    ::Coordinates{:Extrinsic},
    ::Coordinates{:ExtrinsicStiefel},
    x,
)
    [x; 1]
end
function coord_change(
    A::AffineSubspace,
    ::Coordinates{:IntrinsicStiefel},
    ::Coordinates{:Intrinsic},
    u,
)
    γ = A.intrinsic.Y[end, end]
    [u[i] / (γ * u[end]) for i = 1:(length(u)-1)]
end
function coord_change(
    A::AffineSubspace,
    ::Coordinates{:IntrinsicStiefel},
    ::Coordinates{:Extrinsic},
    u,
)
    x = coord_change(A, IntrinsicStiefel(), ExtrinsicStiefel(), u)
    coord_change(A, ExtrinsicStiefel(), Extrinsic(), x)
end
function coord_change(
    A::AffineSubspace,
    ::Coordinates{:IntrinsicStiefel},
    ::Coordinates{:ExtrinsicStiefel},
    u,
)
    A(u, IntrinsicStiefel)
end
function coord_change(
    A::AffineSubspace,
    ::Coordinates{:ExtrinsicStiefel},
    ::Coordinates{:Intrinsic},
    x,
)
    u = coord_change(A, ExtrinsicStiefel(), IntrinsicStiefel(), x)
    coord_change(A, IntrinsicStiefel(), Intrinsic(), u)
end
function coord_change(
    A::AffineSubspace,
    ::Coordinates{:ExtrinsicStiefel},
    ::Coordinates{:Extrinsic},
    x,
)
    x[1:end-1] ./ x[end]
end
function coord_change(
    A::AffineSubspace,
    ::Coordinates{:ExtrinsicStiefel},
    ::Coordinates{:IntrinsicStiefel},
    x,
)
    A.intrinsic.Y' * x
end


# geodesic
geodesic_svd(A::AffineSubspace, B::AffineSubspace) = geodesic_svd(A.intrinsic, B.intrinsic)
function geodesic_svd(A::AffineIntrinsic, B::AffineIntrinsic)
    U, Σ, V = LA.svd(A.Y' * B.Y)
    # Θ = acos.(min.(Σ, 1.0))
    # We have acos(1-eps()) = 2.1073424255447017e-8
    # So this is numerically super unstable if the singular value is off by only eps()
    # We try to correct this fact by treating σ >= 1 - 2eps() as 1.0
    Θ = map(σ -> σ + 2eps() > 1.0 ? 0.0 : acos(σ), Σ)

    n, k = size(A.Y)
    # M = (LA.I - A.Y * A.Y') * B.Y * inv(A.Y' * B.Y)
    # Have to compute an SVD of M s.t. M = Q tanΘ U'
    # Equivalently M * U = Q tan(Θ)
    # We can achieve this by using a *pivoted* QR
    # since then R will be a diagonal matrix s.t. the absolute value of R_ii is θ_{k-i}

    # inv(A.Y' * B.Y) * U = V * LA.diagm(inv.(Σ))
    MU = (LA.I - A.Y * A.Y') * B.Y * V * LA.diagm(inv.(Σ))
    Q, R = LA.qr(MU, Val(true))
    # correct signs and ordering of Q
    Q′ = Q[:, k:-1:1]
    for j = 1:k
        # Look if we need to flip signs
        real(R[k-j+1, k-j+1]) < 0 || continue
        for i = 1:n
            Q′[i, j] = -Q′[i, j]
        end
    end

    Q′, Θ, U
end

geodesic(A::AffineSubspace, B::AffineSubspace) = geodesic(A.intrinsic, B.intrinsic)
function geodesic(A::AffineIntrinsic, B::AffineIntrinsic)
    Q, Θ, U = geodesic_svd(A, B)
    t -> A.Y * U * LA.diagm(cos.(t .* Θ)) * U' + Q * LA.diagm(sin.(t .* Θ)) * U'
end
#
