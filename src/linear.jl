export AffineSubspace,
    AffineExtrinsic,
    AffineIntrinsic,
    Intrinsic,
    Extrinsic,
    intrinsic,
    extrinsic,
    dim,
    codim,
    ambient_dim,
    rand_affine_subspace,
    coord_change,
    geodesic_distance,
    geodesic

# References
const _LKK19 = """Lim, Lek-Heng, Ken Sze-Wai Wong, and Ke Ye. "Numerical algorithms on the affine Grassmannian." SIAM Journal on Matrix Analysis and Applications 40.2 (2019): 371-393"""

"""
    Coordinates

A type used for encoding the used coordinates and for performing coordinate changes.

Currently supported coordinates are:
- [`Intrinsic`](@ref)
- [`Extrinsic`](@ref)
"""
struct Coordinates{T} end
Base.broadcastable(C::Coordinates) = Ref(C)

"""
    Intrinsic <: Coordinates

Indicates the use of the intrinsic description of an affine subspace. See also [`AffineIntrinsic`](@ref).
"""
const Intrinsic = Coordinates{:Intrinsic}()

"""
    Extrinsic <: Coordinates

Indicates the use of the extrinsic description of an affine subspace. See also [`AffineExtrinsic`](@ref).
"""
const Extrinsic = Coordinates{:Extrinsic}()
# These are not part of the public API
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

"""
    AffineExtrinsic

Extrinsic description of an ``m``-dimensional affine subspace ``L`` in ``n``-dimensional space.
That is ``L = { x | A x = b }``.
"""
struct AffineExtrinsic{T}
    A::Matrix{T}
    b::Vector{T}
end
(A::AffineExtrinsic)(x::AbstractVector) = A.A * x - A.b

"""
    dim(A::AffineExtrinsic)

Dimension of the affine subspace `A`.
"""
dim(A::AffineExtrinsic) = size(A.A, 2) - size(A.A, 1)

"""
    codim(A::AffineExtrinsic)

Codimension of the affine subspace `A`.
"""
codim(A::AffineExtrinsic) = size(A.A, 1)

function Base.copy!(A::AffineExtrinsic, B::AffineExtrinsic)
    copy!(A.A, B.A)
    copy!(A.b, B.b)
    A
end

function Base.show(io::IO, mime::MIME"text/plain", A::AffineExtrinsic{T}) where {T}
    println(io, "AffineExtrinsic{$T}:")
    println(io, "A:")
    show(io, mime, A.A)
    println(io, "\nb:")
    show(io, mime, A.b)
end

"""
    AffineIntrinsic

Intrinsic description of an ``m``-dimensional affine subspace ``L`` in ``n``-dimensional space.
That is ``L = { u | A u + b₀ }``. Here, ``A`` and ``b₀`` are in orthogonal coordinates.
That is, the columns of ``A`` are orthonormal and ``A' b₀ = 0``.
"""
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

"""
    dim(A::AffineIntrinsic)

Dimension of the affine subspace `A`.
"""
dim(I::AffineIntrinsic) = size(I.A, 2)

"""
    codim(A::AffineIntrinsic)

Codimension of the affine subspace `A`.
"""
codim(I::AffineIntrinsic) = size(I.A, 1) - size(I.A, 2)

function Base.show(io::IO, mime::MIME"text/plain", A::AffineIntrinsic{T}) where {T}
    println(io, "AffineIntrinsic{$T}:")
    println(io, "A:")
    show(io, mime, A.A)
    println(io, "\nb₀:")
    show(io, mime, A.b₀)
end

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
    A = Matrix((@view svd.U[:, n+1:end])')
    b = A * I.b₀
    AffineExtrinsic(A, b)
end

"""
    AffineSubspace(A, b)

An ``m``-dimensional affine subspace ``L`` in ``n``-dimensional space given by
the extrinsic description ``L = { x | A x = b }``.

```
julia> A = AffineSubspace([1 0 3; 2 1 3], [5, -2])
1-dim. affine subspace {x|Ax=b} with eltype Float64:
A:
2×3 Array{Float64,2}:
 1.0  0.0  3.0
 2.0  1.0  3.0
b:
2-element Array{Float64,1}:
  5.0
 -2.0

julia> dim(A)
1

julia> codim(A)
2

julia> ambient_dim(A)
3
```

An `AffineSubspace` holds always its [`extrinsic`](@ref) description, see also [`AffineIntrinsic`](@ref),
as well as its [`intrinsic`](@ref) description, see [`AffineExtrinsic`](@ref).

```
julia> intrinsic(A)
AffineIntrinsic{Float64}:
A:
3×1 Array{Float64,2}:
 -0.6882472016116853
  0.6882472016116853
  0.22941573387056186
b₀:
3-element Array{Float64,1}:
 -3.0526315789473677
 -3.947368421052632
  2.684210526315789
```

An `AffineSubspace` can be evaluated with either using [`Intrinsic`](@ref)
or [`Extrinsic`](@ref) coordinates.

```
julia> u = [0.5]
1-element Array{Float64,1}:
 0.5

julia> x = A(u, Intrinsic)
3-element Array{Float64,1}:
 -3.3967551797532103
 -3.6032448202467893
  2.79891839325107

julia> A(x, Extrinsic)
  2-element Array{Float64,1}:
   0.0
   0.0
```

To change the used coordinates you can use [`coord_change`](@ref).
```
julia> coord_change(A, Extrinsic, Intrinsic, x)
1-element Array{Float64,1}:
 0.49999999999999994
```
"""
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

    AffineSubspace(AffineExtrinsic(Matrix(float.(A)), Vector(float.(b))))
end

Base.broadcastable(A::AffineSubspace) = Ref(A)

"""
    dim(A::AffineSubspace)

Dimension of the affine subspace `A`.
"""
dim(A::AffineSubspace) = dim(A.intrinsic)

"""
    codim(A::AffineSubspace)

Codimension of the affine subspace `A`.
"""
codim(A::AffineSubspace) = codim(A.intrinsic)

"""
    ambient_dim(A::AffineSubspace)

Dimension of ambient space of the affine subspace `A`.
"""
ambient_dim(A::AffineSubspace) = dim(A) + codim(A)

function Base.show(io::IO, mime::MIME"text/plain", A::AffineSubspace{T}) where {T}
    println(io, "$(dim(A))-dim. affine subspace {x|Ax=b} with eltype $T:")
    println(io, "A:")
    show(io, mime, A.extrinsic.A)
    println(io, "\nb:")
    show(io, mime, A.extrinsic.b)
end

"""
    intrinsic(A::AffineSubspace)

Obtain the intrinsic description of `A`, see also [`AffineIntrinsic`](@ref).
"""
intrinsic(A::AffineSubspace) = A.intrinsic

"""
    extrinsic(A::AffineSubspace)

Obtain the extrinsic description of `A`, see also [`AffineExtrinsic`](@ref).
"""
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

"""
    rand_affine_subspace(n::Integer; dim | codim, real = false)

Generate a random [`AffineSubspace`](@ref) with given dimension `dim` or codimension `codim`
(one of them has to be provided) in ambient space of dimension `n`.
If `real` is `true`, then the extrinsic description is real.
The subspace is generated by drawing each entry of the extrinsic description indepdently
from a normal distribuation using [`randn`](@ref).

## Example
```
julia> rand_affine_subspace(3; dim = 1)
1-dim. affine subspace {x|Ax=b} with eltype Complex{Float64}:
A:
2×3 Array{Complex{Float64},2}:
  -1.73825+1.27987im   -0.0871343+0.840408im  -0.551957+0.106397im
 -0.597132-0.343965im   -0.122543-0.172715im   -1.04949+0.370917im
b:
2-element Array{Complex{Float64},1}:
  0.47083334430689394 + 0.8099804422599071im
 -0.12018696822943896 + 0.11723026326952792im

julia> rand_affine_subspace(4; codim = 1)
3-dim. affine subspace {x|Ax=b} with eltype Complex{Float64}:
A:
1×4 Array{Complex{Float64},2}:
 0.345705+0.0893881im  -0.430867-0.663249im  0.979969-0.569378im  -0.29722-0.192493im
b:
1-element Array{Complex{Float64},1}:
 0.7749708228192062 + 0.9762873764567546im
```
"""
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
    T = real ? Float64 : ComplexF64
    A = randn(T, n - k, n)
    b = randn(T, n - k)

    AffineSubspace(A, b)
end

# Coordinate changes
"""
    coord_change(A::AffineSubspace, C₁::Coordinates, C₂::Coordinates, p)

Given an affine subspace `A` and a point `p` in coordinates `C₁` compute the point
`x` describing p in coordinates `C₂`.

## Example
```
julia> A = AffineSubspace([1 0 3; 2 1 3], [5, -2]);

julia> u = [1.25];

julia> x = coord_change(A, Intrinsic, Extrinsic, u)
3-element Array{Float64,1}:
 -3.9129405809619744
 -3.087059419038025
  2.9709801936539915

julia> A(x, Extrinsic)
2-element Array{Float64,1}:
 0.0
 0.0

julia> x - A(u, Intrinsic)
3-element Array{Float64,1}:
 0.0
 0.0
 0.0
```
"""
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

"""
    geodesic_distance(A::AffineSubspace, B::AffineSubspace)

Compute the geodesic distance between `A` and `B` in the affine Grassmanian `Graff(k, n)`
where `k = dim(A)` and `n` is the amebient dimension.
This follows the derivation in [LKK19].

[LKK19]: $_LKK19.
"""
geodesic_distance(A::AffineSubspace, B::AffineSubspace) =
    geodesic_distance(A.intrinsic, B.intrinsic)
function geodesic_distance(A::AffineIntrinsic, B::AffineIntrinsic)
    sqrt(sum(σᵢ -> acos(min(σᵢ, 1.0))^2, LA.svdvals(A.Y' * B.Y)))
end

# geodesic
"""
    geodesic_svd(A::AffineSubspace, B::AffineSubspace)

Computes the factors ``Q``, ``Θ`` and ``U`` from Corollary 4.3 in [LKK19].
These values are necessary to construct the geodesic between `A` and `B`.

[LKK19]: $_LKK19
"""
geodesic_svd(A::AffineSubspace, B::AffineSubspace) = geodesic_svd(A.intrinsic, B.intrinsic)
function geodesic_svd(A::AffineIntrinsic, B::AffineIntrinsic)
    U, Σ, V = LA.svd(A.Y' * B.Y)
    # Θ = acos.(min.(Σ, 1.0))
    # We have acos(1-eps()) = 2.1073424255447017e-8
    # So this is numerically super unstable if the singular value is off by only eps()
    # We try to correct this fact by treating σ >= 1 - 2eps() as 1.0
    Θ = map(σ -> σ + 2 * eps() > 1.0 ? 0.0 : acos(σ), Σ)

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

"""
    geodesic(A::AffineSubspace, B::AffineSubspace)

Returns the geodesic ``γ(t)`` connecting `A` and `B` in the Grassmanian ``Gr(k+1,n+1)``
where ``k`` is the dimension of ``A`` and ``n`` is the ambient dimension.
See also Corollary 4.3 in [LKK19].

[LKK19]: $_LKK19
"""
geodesic(A::AffineSubspace, B::AffineSubspace) = geodesic(A.intrinsic, B.intrinsic)
function geodesic(A::AffineIntrinsic, B::AffineIntrinsic)
    Q, Θ, U = geodesic_svd(A, B)
    t -> A.Y * U * LA.diagm(cos.(t .* Θ)) * U' + Q * LA.diagm(sin.(t .* Θ)) * U'
end
#
