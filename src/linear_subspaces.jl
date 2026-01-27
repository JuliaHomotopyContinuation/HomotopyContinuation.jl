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
            if size(A, 1) == 0
                new{T}(A, b)
            else
                svd = LA.svd(A)
                new{T}(svd.Vt, (inv.(svd.S) .* (svd.U' * b)))
            end
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
    IntrinsicDescription(A, b)

Intrinsic description of an ``m``-dimensional (affine) linear subspace ``L`` in ``n``-dimensional space.
That is ``L = \\{ u | A u + b \\}``. Here, ``A`` and ``b`` are in orthogonal coordinates.
That is, the columns of ``A`` are orthonormal and ``A' b = 0``.
"""
struct IntrinsicDescription{T}
    # orthogonal coordinates
    A::Matrix{T}
    b::Vector{T}
    # stiefel coordinates for A
    X::Matrix{T}
    # stiefel coordinates for [A b]
    Y::Matrix{T}
end
IntrinsicDescription(A::Matrix{T1}, b::Vector{T2}) where {T1,T2} =
    IntrinsicDescription(A, b, stiefel_coordinates_intrinsic(A), stiefel_coordinates_intrinsic(A, b))

function Base.convert(::Type{IntrinsicDescription{T}}, A::IntrinsicDescription) where {T}
    IntrinsicDescription(
        convert(Matrix{T}, A.A),
        convert(Vector{T}, A.b),
    )
end

(A::IntrinsicDescription)(u::AbstractVector, ::Coordinates{:Intrinsic}) = A.A * u + A.b

function stiefel_coordinates_intrinsic(A::AbstractMatrix)
    n, k = size(A)
    SVD = LA.svd!(A)
    SVD.U
end
function stiefel_coordinates_intrinsic!(X, A::AbstractMatrix)
    n, k = size(A)
    SVD = LA.svd!(A)
    X .= SVD.U
    X
end
function stiefel_coordinates_intrinsic(A::AbstractMatrix, b::AbstractVector)
    n, k = size(A)
    Y = zeros(eltype(A), n + 1, k + 1)
    stiefel_coordinates_intrinsic!(Y, A, b)
    Y
end
function stiefel_coordinates_intrinsic!(Y, A::AbstractMatrix, b::AbstractVector)
    γ = sqrt(1 + sum(abs2, b))
    n, k = size(A)
    Y[1:n, 1:k] .= A
    Y[1:n, k+1] .= b ./ γ
    Y[n+1, k+1] = 1 / γ
    SVD = LA.svd!(Y)
    Y .= SVD.U
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
    A.A == B.A && A.b == B.b && A.Y == B.Y
end

function Base.show(io::IO, A::IntrinsicDescription{T}) where {T}
    println(io, "IntrinsicDescription{$T}:")
    println(io, "A:")
    show(io, A.A)
    println(io, "\nb:")
    show(io, A.b)
end

function Base.copy!(A::IntrinsicDescription, B::IntrinsicDescription)
    copy!(A.A, B.A)
    copy!(A.b, B.b)
    copy!(A.Y, B.Y)
    A
end
Base.copy(A::IntrinsicDescription) = IntrinsicDescription(copy(A.A), copy(A.b), copy(A.Y))

Base.broadcastable(A::IntrinsicDescription) = Ref(A)

function IntrinsicDescription(E::ExtrinsicDescription)
    svd = LA.svd(E.A; full = true)
    m, n = size(E.A)
    A = Matrix((@view svd.Vt[m+1:end, :])')
    if iszero(E.b)
        b = zeros(eltype(E.b), size(A, 1))
    else
        b = svd \ E.b
    end
    IntrinsicDescription(A, b)
end

function ExtrinsicDescription(I::IntrinsicDescription)
    svd = LA.svd(I.A; full = true)
    m, n = size(I.A)
    A = Matrix((@view svd.U[:, n+1:end])')
    if iszero(I.b)
        b = zeros(eltype(A), size(A, 1))
    else
        b = A * I.b
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
b:
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
    0 ≤ size(A, 1) ≤ size(A, 2) || throw(
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
function (A::LinearSubspace)(x::AbstractVector, ::Coordinates{:Extrinsic} = Extrinsic)
    extrinsic(A)(x)
end

"""
    rand_subspace(n::Integer; dim | codim, affine = true, real = false)

Generate a random [`LinearSubspace`](@ref) with given dimension `dim` or codimension `codim`
(one of them has to be provided) in ambient space of dimension `n`.
If `real` is `true`, then the extrinsic description is real.
If `affine` then an affine linear subspace is generated.
The subspace is generated by drawing the matrix `A` the extrinsic description independently
from a normal distribuation using `randn`. The offset is drawn from the normal distribution scaled by `sqrt(n)`, where `n` is the number of columns of `A`. 

    rand_subspace(x::AbstractVector; dim | codim, affine = true)

Generate a random [`LinearSubspace`](@ref) with given dimension `dim` or codimension `codim`
(one of them has to be provided) in ambient space of dimension `length(x)` going through
the given point `x`.

## Example

Construction of a general random subspace:
```julia-repl
julia> rand_subspace(3; dim = 1)
1-dim. affine linear subspace {x | Ax=b} with eltype ComplexF64:
A:
ComplexF64[-0.18709056988162928 - 0.021159068437462684im 0.3448082517062497 - 0.5619954950326371im 0.12814427601487394 + 0.716517124796854im; 0.6274330537397367 + 0.4327848781661207im 0.3433265716931795 - 0.4388117939096428im -0.09290449999579178 - 0.3161721696818136im]
b:
ComplexF64[0.1568834271069519 - 0.7893604951164849im, -1.068448434670914 + 3.422846499360692im]

julia> rand_subspace(4; codim = 1)
3-dim. affine linear subspace {x | Ax=b} with eltype ComplexF64:
A:
ComplexF64[-0.16031107752943963 - 0.6372670776364352im 0.4895718108240883 + 0.4093061519132399im -0.32944285959791164 + 0.12384710749057845im -0.17790282756300982 - 0.07388387107969129im]
b:
ComplexF64[-0.09372015908519082 + 0.3051873249600184im]
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
        LinearSubspace(A, sqrt(n) .* randn(T, n - k))
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
coord_change(A::LinearSubspace, ::Coordinates{:Extrinsic}, ::Coordinates{:Intrinsic}, x) =
    A.intrinsic.A' * (x - A.intrinsic.b)


"""
    geodesic_distance(V::LinearSubspace, W::LinearSubspace)

Compute the geodesic distance between `V = {x | Ax = a}` and `W= {x | Bx = b}` as `sqrt(d^2 + ||a-b||^2)`, where `d` is the distance from the columnspan of `A` to the columnspan of `B` in the Grassmannian. This follows the derivation in [^LKK19].

[^LKK19]: $_LKK19.
"""
geodesic_distance(A::LinearSubspace, B::LinearSubspace) =
    geodesic_distance(A.intrinsic, B.intrinsic)
function geodesic_distance(A::IntrinsicDescription, B::IntrinsicDescription)
    sqrt(sum(σᵢ -> acos(min(σᵢ, 1.0))^2, LA.svdvals(A.Y' * B.Y)))
end

# geodesic
"""
    grassmannian_svd(A::LinearSubspace, B::LinearSubspace)

Computes the factors ``Q``, ``Θ`` and ``U`` from Corollary 4.3 in [^LKK19].
These values are necessary to construct the path from `A` and `B`.

[^LKK19]: $_LKK19
"""
grassmannian_svd(A::LinearSubspace, B::LinearSubspace) = grassmannian_svd(extrinsic(A), extrinsic(B))
grassmannian_svd(A::E1, B::E2) where {E1, E2 <: ExtrinsicDescription}  = grassmannian_svd(transpose(A.A), transpose(B.A))
function grassmannian_svd(A::I1, B::I2; 
                        embedded_projective::Bool = false) where {I1, I2 <: IntrinsicDescription}
         if !embedded_projective
            grassmannian_svd(A.X, B.X)
         else
            grassmannian_svd(A.Y, B.Y)
         end
end

function grassmannian_svd(A::AbstractMatrix, B::AbstractMatrix)

    # here A and B are assumed to be Stiefel matrices representing linear spaces.
    n, k = size(A)
    U, Σ, V = LA.svd!(A' * B)
    # M = (LA.I - A * A') * B
    # Have to compute an SVD of M s.t. M = Q sinΘ V'
    # Equivalently M * V = Q sin(Θ)
    # We can achieve this by using a *pivoted* QR
    # since then R will be a diagonal matrix s.t. the absolute value of R_ii is θ_{k-i}
    MV = (LA.I - A * A') * B * V
    
    # Θ = acos.(min.(Σ, 1.0))
    # We have acos(1-eps()) = 2.1073424255447017e-8
    # So this is numerically super unstable if the singular value is off by only eps()
    # We try to correct this fact by treating σ >= 1 - 2eps() as 1.0
    Θ = map(σ -> σ + 2 * eps() > 1.0 ? 0.0 : acos(σ), Σ)

    Q, R = qr_col_norm!(MV)
    # correct signs and ordering of Q
    Q′ = Q[:, k:-1:1]
    for j = 1:k
        # Look if we need to flip signs
        real(R[k-j+1, k-j+1]) < 0 || continue
        for i = 1:n
            Q′[i, j] = -Q′[i, j]
        end
    end

    # scale with c (if c is random complex, then we get generic homotopies)
    Q′, Θ, U
end

"""
    geodesic(V::LinearSubspace, W::LinearSubspace)

    Compute the geodesic distance between `V = {x | Ax = a}` and `W= {x | Bx = b}` as `sqrt(d^2 + ||a-b||^2)`, where `d` is the distance from the columnspan of `A` to the columnspan of `B` in the Grassmannian. This follows the derivation in [^LKK19].

Returns the geodesic ``γ(t)`` by connecting `V = {x | Ax = a}` and `W= {x | Bx = b}`. `A` and `B` are interpolated in the Grassmannian using Stiefel coordinates. See also Corollary 4.3 in [^LKK19]. `a` and `b` are interpolated linearly. Thus, ``γ(t)`` is a geodesic in the Euclidean group.

[^LKK19]: $_LKK19
"""
geodesic(A::LinearSubspace, B::LinearSubspace) = geodesic(A.intrinsic, B.intrinsic)
function geodesic(A::IntrinsicDescription, B::IntrinsicDescription)
    Q, Θ, U = grassmannian_svd(A, B)
    t -> IntrinsicDescription(
            A.X * U * LA.diagm(cos.(t .* Θ,)) * U' + Q * LA.diagm(sin.(t .* Θ,)) * U', 
            t .* A.b + (1-t) .* B.b
            ) |> LinearSubspace
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
    LA.mul!(int.b, ext.A', δb, true, true)
    stiefel_coordinates_intrinsic!(int.X, int.A)
    stiefel_coordinates_intrinsic!(int.Y, int.A, int.b)
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

