export LinearSubspace,
    ExtrinsicDescription,
    IntrinsicDescription,
    Coordinates,
    Intrinsic,
    Extrinsic,
    intrinsic,
    is_linear,
    extrinsic,
    dim,
    codim,
    ambient_dim,
    rand_subspace,
    coord_change,
    geodesic_distance,
    geodesic,
    translate

# References
const _LKK19 = """Lim, Lek-Heng, Ken Sze-Wai Wong, and Ke Ye. "Numerical algorithms on the affine Grassmannian." SIAM Journal on Matrix Analysis and Applications 40.2 (2019): 371-393"""

abstract type AbstractSubspace{T} end

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

Indicates the use of the intrinsic description of an (affine) linear subspace. See also [`IntrinsicDescription`](@ref).
"""
const Intrinsic = Coordinates{:Intrinsic}()

"""
    Extrinsic <: Coordinates

Indicates the use of the extrinsic description of an (affine) linear subspace. See also [`ExtrinsicDescription`](@ref).
"""
const Extrinsic = Coordinates{:Extrinsic}()
# These are not part of the public API
const IntrinsicStiefel = Coordinates{:IntrinsicStiefel}()
const ExtrinsicStiefel = Coordinates{:ExtrinsicStiefel}()

"""
    ExtrinsicDescription(A, b)

Extrinsic description of an ``m``-dimensional (affine) linear subspace ``L`` in ``n``-dimensional space.
That is ``L = \\{ x | A x = b \\}``.
Note that internally `A` and `b` will be stored such that the rows of `A` are orthonormal.
"""
struct ExtrinsicDescription{T}
    A::Matrix{T}
    b::Vector{T}

    function ExtrinsicDescription(
        A::Matrix{T},
        b::Vector{T};
        orthonormal::Bool = false,
    ) where {T}
        if orthonormal
            new{T}(A, b)
        else
            svd = LA.svd(A)
            new{T}(svd.Vt, (inv.(svd.S) .* (svd.U' * b)))
        end
    end
end
(A::ExtrinsicDescription)(x::AbstractVector) = A.A * x - A.b

function Base.convert(::Type{ExtrinsicDescription{T}}, A::ExtrinsicDescription) where {T}
    ExtrinsicDescription(
        convert(Matrix{T}, A.A),
        convert(Vector{T}, A.b);
        orthonormal = true,
    )
end

"""
    dim(A::ExtrinsicDescription)

Dimension of the (affine) linear subspace `A`.
"""
dim(A::ExtrinsicDescription) = size(A.A, 2) - size(A.A, 1)

"""
    codim(A::ExtrinsicDescription)

Codimension of the (affine) linear subspace `A`.
"""
codim(A::ExtrinsicDescription) = size(A.A, 1)

Base.:(==)(A::ExtrinsicDescription, B::ExtrinsicDescription) = A.A == B.A && A.b == B.b
function Base.copy!(A::ExtrinsicDescription, B::ExtrinsicDescription)
    copy!(A.A, B.A)
    copy!(A.b, B.b)
    A
end
Base.copy(A::ExtrinsicDescription) =
    ExtrinsicDescription(copy(A.A), copy(A.b); orthonormal = true)

function Base.show(io::IO, A::ExtrinsicDescription{T}) where {T}
    println(io, "ExtrinsicDescription{$T}:")
    println(io, "A:")
    show(io, A.A)
    println(io, "\nb:")
    show(io, A.b)
end

Base.broadcastable(A::ExtrinsicDescription) = Ref(A)

"""
    IntrinsicDescription(A, b₀)

Intrinsic description of an ``m``-dimensional (affine) linear subspace ``L`` in ``n``-dimensional space.
That is ``L = \\{ u | A u + b₀ \\}``. Here, ``A`` and ``b₀`` are in orthogonal coordinates.
That is, the columns of ``A`` are orthonormal and ``A' b₀ = 0``.
"""
struct IntrinsicDescription{T}
    # orthogonal coordinates
    A::Matrix{T}
    b₀::Vector{T}
    # stiefel coordinates
    Y::Matrix{T}
end
IntrinsicDescription(A::Matrix, b₀::Vector) =
    IntrinsicDescription(A, b₀, stiefel_coordinates(A, b₀))

function Base.convert(::Type{IntrinsicDescription{T}}, A::IntrinsicDescription) where {T}
    IntrinsicDescription(
        convert(Matrix{T}, A.A),
        convert(Vector{T}, A.b₀),
        convert(Matrix{T}, A.Y),
    )
end

(A::IntrinsicDescription)(u::AbstractVector, ::Coordinates{:Intrinsic}) = A.A * u + A.b₀
(A::IntrinsicDescription)(u::AbstractVector, ::Coordinates{:IntrinsicStiefel}) = A.Y * u

function stiefel_coordinates(A::AbstractMatrix, b::AbstractVector)
    n, k = size(A)
    Y = zeros(eltype(A), n + 1, k + 1)
    stiefel_coordinates!(Y, A, b)
    Y
end
function stiefel_coordinates!(Y, A::AbstractMatrix, b::AbstractVector)
    γ = sqrt(1 + sum(abs2, b))
    n, k = size(A)
    Y[1:n, 1:k] .= A
    Y[1:n, k+1] .= b ./ γ
    Y[n+1, k+1] = 1 / γ
    Y
end

"""
    dim(A::IntrinsicDescription)

Dimension of the (affine) linear subspace `A`.
"""
dim(I::IntrinsicDescription) = size(I.A, 2)

"""
    codim(A::IntrinsicDescription)

Codimension of the (affine) linear subspace `A`.
"""
codim(I::IntrinsicDescription) = size(I.A, 1) - size(I.A, 2)

function Base.:(==)(A::IntrinsicDescription, B::IntrinsicDescription)
    A.A == B.A && A.b₀ == B.b₀ && A.Y == B.Y
end

function Base.show(io::IO, A::IntrinsicDescription{T}) where {T}
    println(io, "IntrinsicDescription{$T}:")
    println(io, "A:")
    show(io, A.A)
    println(io, "\nb₀:")
    show(io, A.b₀)
end

function Base.copy!(A::IntrinsicDescription, B::IntrinsicDescription)
    copy!(A.A, B.A)
    copy!(A.b₀, B.b₀)
    copy!(A.Y, B.Y)
    A
end
Base.copy(A::IntrinsicDescription) = IntrinsicDescription(copy(A.A), copy(A.b₀), copy(A.Y))

Base.broadcastable(A::IntrinsicDescription) = Ref(A)

function IntrinsicDescription(E::ExtrinsicDescription)
    svd = LA.svd(E.A; full = true)
    m, n = size(E.A)
    A = Matrix((@view svd.Vt[m+1:end, :])')
    if iszero(E.b)
        b₀ = zeros(eltype(E.b), size(A, 1))
    else
        b₀ = svd \ E.b
    end
    IntrinsicDescription(A, b₀)
end

function ExtrinsicDescription(I::IntrinsicDescription)
    svd = LA.svd(I.A; full = true)
    m, n = size(I.A)
    A = Matrix((@view svd.U[:, n+1:end])')
    if iszero(I.b₀)
        b = zeros(eltype(A), size(A, 1))
    else
        b = A * I.b₀
    end
    ExtrinsicDescription(A, b; orthonormal = true)
end

"""
    LinearSubspace(A, b)

An ``m``-dimensional (affine) linear subspace ``L`` in ``n``-dimensional space given by
the extrinsic description ``L = \\{ x | A x = b \\}``.

```julia-repl
julia> A = LinearSubspace([1 0 3; 2 1 3], [5, -2])
1-dim. (affine) linear subspace {x|Ax=b} with eltype Float64:
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

A `LinearSubspace` holds always its [`extrinsic`](@ref) description,
see also [`ExtrinsicDescription`](@ref), as well as its [`intrinsic`](@ref) description,
see [`IntrinsicDescription`](@ref).

```julia-repl
julia> intrinsic(A)
IntrinsicDescription{Float64}:
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

A `LinearSubspace` can be evaluated with either using [`Intrinsic`](@ref)
or [`Extrinsic`](@ref) coordinates.

```julia-repl
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
```julia-repl
julia> coord_change(A, Extrinsic, Intrinsic, x)
1-element Array{Float64,1}:
 0.49999999999999994
```
"""
struct LinearSubspace{T} <: AbstractSubspace{T}
    extrinsic::ExtrinsicDescription{T}
    intrinsic::IntrinsicDescription{T}
end

LinearSubspace(I::IntrinsicDescription) = LinearSubspace(ExtrinsicDescription(I), I)
LinearSubspace(E::ExtrinsicDescription) = LinearSubspace(E, IntrinsicDescription(E))

function LinearSubspace(
    A::AbstractMatrix{T},
    b::AbstractVector{T} = zeros(eltype(A), size(A, 1)),
) where {T}
    size(A, 1) == length(b) || throw(ArgumentError("Size of A and b not compatible."))
    0 < size(A, 1) ≤ size(A, 2) || throw(
        ArgumentError(
            "Affine subspace has to be given in extrinsic coordinates, i.e., by A x = b.",
        ),
    )

    LinearSubspace(ExtrinsicDescription(Matrix(float.(A)), Vector(float.(b))))
end

function Base.convert(::Type{LinearSubspace{T}}, A::LinearSubspace) where {T}
    LinearSubspace(
        convert(ExtrinsicDescription{T}, A.extrinsic),
        convert(IntrinsicDescription{T}, A.intrinsic),
    )
end

Base.broadcastable(A::LinearSubspace) = Ref(A)

"""
    dim(A::LinearSubspace)

Dimension of the (affine) linear subspace `A`.
"""
dim(A::LinearSubspace) = dim(A.intrinsic)

"""
    codim(A::LinearSubspace)

Codimension of the (affine) linear subspace `A`.
"""
codim(A::LinearSubspace) = codim(A.intrinsic)

"""
    ambient_dim(A::LinearSubspace)

Dimension of ambient space of the (affine) linear subspace `A`.
"""
ambient_dim(A::LinearSubspace) = dim(A) + codim(A)

"""
    is_linear(L::LinearSubspace)

Returns `true` if the space is proper linear subspace, i.e., described by
`L = \\{ x | Ax = 0 \\}.
"""
is_linear(A::LinearSubspace) = iszero(extrinsic(A).b)

function Base.show(io::IO, A::LinearSubspace{T}) where {T}
    if is_linear(A)
        println(io, "$(dim(A))-dim. linear subspace {x | Ax=0} with eltype $T:")
        println(io, "A:")
        show(io, A.extrinsic.A)
    else
        println(io, "$(dim(A))-dim. affine linear subspace {x | Ax=b} with eltype $T:")
        println(io, "A:")
        show(io, A.extrinsic.A)
        println(io, "\nb:")
        show(io, A.extrinsic.b)
    end
end

"""
    intrinsic(A::LinearSubspace)

Obtain the intrinsic description of `A`, see also [`IntrinsicDescription`](@ref).
"""
intrinsic(A::LinearSubspace) = A.intrinsic

"""
    extrinsic(A::LinearSubspace)

Obtain the extrinsic description of `A`, see also [`ExtrinsicDescription`](@ref).
"""
extrinsic(A::LinearSubspace) = A.extrinsic

function Base.copy!(A::LinearSubspace, B::LinearSubspace)
    copy!(A.intrinsic, B.intrinsic)
    copy!(A.extrinsic, B.extrinsic)
    A
end
Base.copy(A::LinearSubspace) = LinearSubspace(copy(A.extrinsic), copy(A.intrinsic))

function Base.:(==)(A::LinearSubspace, B::LinearSubspace)
    intrinsic(A) == intrinsic(B) && extrinsic(A) == extrinsic(B)
end
Base.isequal(A::LinearSubspace, B::LinearSubspace) = A === B

function (A::LinearSubspace)(x::AbstractVector, ::Coordinates{:Intrinsic})
    intrinsic(A)(x, Intrinsic)
end
function (A::LinearSubspace)(x::AbstractVector, ::Coordinates{:IntrinsicStiefel})
    intrinsic(A)(x, IntrinsicStiefel)
end
function (A::LinearSubspace)(x::AbstractVector, ::Coordinates{:Extrinsic} = Extrinsic)
    extrinsic(A)(x)
end

"""
    rand_subspace(n::Integer; dim | codim, affine = true, real = false)

Generate a random [`LinearSubspace`](@ref) with given dimension `dim` or codimension `codim`
(one of them has to be provided) in ambient space of dimension `n`.
If `real` is `true`, then the extrinsic description is real.
If `affine` then an affine linear subspace is generated.
The subspace is generated by drawing each entry of the extrinsic description indepdently
from a normal distribuation using `randn`.

    rand_subspace(x::AbstractVector; dim | codim, affine = true)

Generate a random [`LinearSubspace`](@ref) with given dimension `dim` or codimension `codim`
(one of them has to be provided) in ambient space of dimension `length(x)` going through
the given point `x`.

## Example

Construction of a general random subspace:
```julia-repl
julia> rand_subspace(3; dim = 1)
1-dim. (affine) linear subspace {x|Ax=b} with eltype Complex{Float64}:
A:
2×3 Array{Complex{Float64},2}:
  -1.73825+1.27987im   -0.0871343+0.840408im  -0.551957+0.106397im
 -0.597132-0.343965im   -0.122543-0.172715im   -1.04949+0.370917im
b:
2-element Array{Complex{Float64},1}:
  0.47083334430689394 + 0.8099804422599071im
 -0.12018696822943896 + 0.11723026326952792im

julia> rand_subspace(4; codim = 1)
3-dim. (affine) linear subspace {x|Ax=b} with eltype Complex{Float64}:
A:
1×4 Array{Complex{Float64},2}:
 0.345705+0.0893881im  -0.430867-0.663249im  0.979969-0.569378im  -0.29722-0.192493im
b:
1-element Array{Complex{Float64},1}:
 0.7749708228192062 + 0.9762873764567546im
```
"""
function rand_subspace(
    n::Integer;
    dim::Union{Nothing,Integer} = nothing,
    codim::Union{Nothing,Integer} = nothing,
    real::Bool = false,
    affine::Bool = true,
)
    !isnothing(dim) ||
        !isnothing(codim) ||
        throw(ArgumentError("Neither `dim` nor `codim` specified."))

    if !isnothing(dim)
        0 < dim < n || throw(ArgumentError("`dim` has to be between 0 and `n`."))
        k = dim
    else
        0 < codim < n || throw(ArgumentError("`codim` has to be between 0 and `n`."))
        k = n - codim
    end
    T = real ? Float64 : ComplexF64
    A = randn(T, n - k, n)
    if affine
        LinearSubspace(A, randn(T, n - k))
    else
        LinearSubspace(A)
    end

end
rand_subspace(x::AbstractVector{Variable}; kwargs...) = rand_subspace(length(x); kwargs...)
function rand_subspace(
    x::AbstractVector;
    dim::Union{Nothing,Integer} = nothing,
    codim::Union{Nothing,Integer} = nothing,
    affine::Bool = true,
)
    n = length(x)
    !isnothing(dim) ||
        !isnothing(codim) ||
        throw(ArgumentError("Neither `dim` nor `codim` specified."))

    if !isnothing(dim)
        0 < dim < n || throw(ArgumentError("`dim` has to be between 0 and `n`."))
        k = dim
    else
        0 < codim < n || throw(ArgumentError("`codim` has to be between 0 and `n`."))
        k = n - codim
    end

    if affine
        A = randn(eltype(x), n - k, n)
        b = A * x
        LinearSubspace(A, b)
    else
        N = LA.nullspace(Matrix(x'))'
        A = randn(eltype(N), n - k, size(N, 1)) * N
        LinearSubspace(A)
    end
end

# Coordinate changes
"""
    coord_change(A::LinearSubspace, C₁::Coordinates, C₂::Coordinates, p)

Given an (affine) linear subspace `A` and a point `p` in coordinates `C₁` compute the point
`x` describing p in coordinates `C₂`.

## Example
```julia-repl
julia> A = LinearSubspace([1 0 3; 2 1 3], [5, -2]);

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
coord_change(A::LinearSubspace, ::C, ::C, x) where {C<:Coordinates} = x
coord_change(A::LinearSubspace, ::Coordinates{:Intrinsic}, ::Coordinates{:Extrinsic}, u) =
    A(u, Intrinsic)
function coord_change(
    A::LinearSubspace,
    ::Coordinates{:Intrinsic},
    ::Coordinates{:IntrinsicStiefel},
    u,
)
    [
        u
        1 / A.intrinsic.Y[end, end]
    ]
end
function coord_change(
    A::LinearSubspace,
    ::Coordinates{:Intrinsic},
    ::Coordinates{:ExtrinsicStiefel},
    u,
)
    v = coord_change(A, Intrinsic, IntrinsicStiefel, u)
    coord_change(A, IntrinsicStiefel, ExtrinsicStiefel, v)
end

coord_change(A::LinearSubspace, ::Coordinates{:Extrinsic}, ::Coordinates{:Intrinsic}, x) =
    A.intrinsic.A' * (x - A.intrinsic.b₀)
function coord_change(
    A::LinearSubspace,
    ::Coordinates{:Extrinsic},
    ::Coordinates{:IntrinsicStiefel},
    x,
)
    u = coord_change(A, Extrinsic, Intrinsic, x)
    coord_change(A, Intrinsic, IntrinsicStiefel, u)
end
function coord_change(
    A::LinearSubspace,
    ::Coordinates{:Extrinsic},
    ::Coordinates{:ExtrinsicStiefel},
    x,
)
    [x; 1]
end
function coord_change(
    A::LinearSubspace,
    ::Coordinates{:IntrinsicStiefel},
    ::Coordinates{:Intrinsic},
    u,
)
    γ = A.intrinsic.Y[end, end]
    [u[i] / (γ * u[end]) for i = 1:(length(u)-1)]
end
function coord_change(
    A::LinearSubspace,
    ::Coordinates{:IntrinsicStiefel},
    ::Coordinates{:Extrinsic},
    u,
)
    x = coord_change(A, IntrinsicStiefel, ExtrinsicStiefel, u)
    coord_change(A, ExtrinsicStiefel, Extrinsic, x)
end
function coord_change(
    A::LinearSubspace,
    ::Coordinates{:IntrinsicStiefel},
    ::Coordinates{:ExtrinsicStiefel},
    u,
)
    A(u, IntrinsicStiefel)
end
function coord_change(
    A::LinearSubspace,
    ::Coordinates{:ExtrinsicStiefel},
    ::Coordinates{:Intrinsic},
    x,
)
    u = coord_change(A, ExtrinsicStiefel, IntrinsicStiefel, x)
    coord_change(A, IntrinsicStiefel, Intrinsic, u)
end
function coord_change(
    A::LinearSubspace,
    ::Coordinates{:ExtrinsicStiefel},
    ::Coordinates{:Extrinsic},
    x,
)
    x[1:end-1] ./ x[end]
end
function coord_change(
    A::LinearSubspace,
    ::Coordinates{:ExtrinsicStiefel},
    ::Coordinates{:IntrinsicStiefel},
    x,
)
    A.intrinsic.Y' * x
end

"""
    geodesic_distance(A::LinearSubspace, B::LinearSubspace)

Compute the geodesic distance between `A` and `B` in the affine Grassmanian `Graff(k, n)`
where `k = dim(A)` and `n` is the amebient dimension.
This follows the derivation in [^LKK19].

[^LKK19]: $_LKK19.
"""
geodesic_distance(A::LinearSubspace, B::LinearSubspace) =
    geodesic_distance(A.intrinsic, B.intrinsic)
function geodesic_distance(A::IntrinsicDescription, B::IntrinsicDescription)
    sqrt(sum(σᵢ -> acos(min(σᵢ, 1.0))^2, LA.svdvals(A.Y' * B.Y)))
end

# geodesic
"""
    geodesic_svd(A::LinearSubspace, B::LinearSubspace)

Computes the factors ``Q``, ``Θ`` and ``U`` from Corollary 4.3 in [^LKK19].
These values are necessary to construct the geodesic between `A` and `B`.

[^LKK19]: $_LKK19
"""
geodesic_svd(A::LinearSubspace, B::LinearSubspace) = geodesic_svd(A.intrinsic, B.intrinsic)
function geodesic_svd(A::IntrinsicDescription, B::IntrinsicDescription)
    U, Σ, V = LA.svd!(A.Y' * B.Y)
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
    Q, R = LA.qr!(MU, Val(true))
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
    geodesic(A::LinearSubspace, B::LinearSubspace)

Returns the geodesic ``γ(t)`` connecting `A` and `B` in the Grassmanian ``Gr(k+1,n+1)``
where ``k`` is the dimension of ``A`` and ``n`` is the ambient dimension.
See also Corollary 4.3 in [^LKK19].

[^LKK19]: $_LKK19
"""
geodesic(A::LinearSubspace, B::LinearSubspace) = geodesic(A.intrinsic, B.intrinsic)
function geodesic(A::IntrinsicDescription, B::IntrinsicDescription)
    Q, Θ, U = geodesic_svd(A, B)
    t -> A.Y * U * LA.diagm(cos.(t .* Θ,)) * U' + Q * LA.diagm(sin.(t .* Θ,)) * U'
end
#
"""
    translate(L::LinearSubspace, δb, ::Coordinates = Extrinsic)

Translate the (affine) linear subspace `L` by `δb`.
"""
function translate(L::LinearSubspace, δb, ::Coordinates{:Extrinsic} = Extrinsic)
    translate!(copy(L), δb, Extrinsic)
end
function translate!(L::LinearSubspace, δb, ::Coordinates{:Extrinsic} = Extrinsic)
    ext = extrinsic(L)
    ext.b .+= δb
    int = intrinsic(L)
    LA.mul!(int.b₀, ext.A', δb, true, true)
    stiefel_coordinates!(int.Y, int.A, int.b₀)
    L
end

"""
    Base.intersect(L₁::LinearSubspace, L₂::LinearSubspace)

Intersect the two given linear subspaces.
Throws an `ErrorException` if the intersection is the sum of the
codimensions is larger than the ambient dimension.
"""
function Base.intersect(L₁::LinearSubspace, L₂::LinearSubspace)
    ambient_dim(L₁) == ambient_dim(L₂) || error("Ambient dimensions don't match.")
    codim(L₁) + codim(L₂) ≤ ambient_dim(L₁) ||
        error("Sum of codimensions larger than ambient dimension.")
    ext₁ = extrinsic(L₁)
    ext₂ = extrinsic(L₂)
    LinearSubspace([ext₁.A; ext₂.A], [ext₁.b; ext₂.b])
end
