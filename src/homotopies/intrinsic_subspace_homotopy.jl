export IntrinsicSubspaceHomotopy, IntrinsicSubspaceStiefelHomotopy, IntrinsicSubspaceOffsetHomotopy, set_subspaces!




struct GrassmannianGeodesic
    Q::Matrix{ComplexF64}
    Q_cos::Matrix{ComplexF64}
    Θ::Vector{Float64}
    U::Matrix{ComplexF64}
    γ1::Matrix{ComplexF64}
end

function GrassmannianGeodesic(start::LinearSubspace, target::LinearSubspace; affine::Bool = true)
    Q, Θ, U = grassmannian_svd(target, start; affine = affine)
    if affine 
        Q_cos = target.intrinsic.X * U
    else
        Q_cos = target.intrinsic.Y * U
    end
    γ1 = similar(Q_cos)
    n, k = size(γ1)
    for j = 1:k
        s, c = sincos(Θ[j])
        for i = 1:n
            γ1[i, j] = Q_cos[i, j] * c + Q[i, j] * s
        end
    end
    G = GrassmannianGeodesic(Q, Q_cos, Θ, U, γ1)
    G
end


const PROJECTIVE_INTRINSIC_LRU = LRU{
    Tuple{LinearSubspace{ComplexF64},LinearSubspace{ComplexF64}},
    GrassmannianGeodesic,
}(
    maxsize = 128,
)

const AFFINE_INTRINSIC_LRU = LRU{
    Tuple{LinearSubspace{ComplexF64},LinearSubspace{ComplexF64}},
    GrassmannianGeodesic,
}(
    maxsize = 128,
)
"""
    IntrinsicSubspaceStiefelHomotopy(F::System, V::LinearSubspace, W::LinearSubspace)
    IntrinsicSubspaceStiefelHomotopy(F::AbstractSystem, V::LinearSubspace, W::LinearSubspace)

Creates a homotopy ``H(x,t) = F(γ(t)x)`` where ``γ(t)`` is a family of affine subspaces
such that ``γ(1) = V`` and ``γ(0) = W``.
Here ``γ(t)`` is the geodesic between `V` and `W` in the affine Grassmanian, i.e.,
it is the curve of minimal length connecting `V` and `W`.
See also [`LinearSubspace`](@ref) and [`geodesic`](@ref) and the references therein.
"""
Base.@kwdef mutable struct IntrinsicSubspaceStiefelHomotopy{S<:AbstractSystem} <: AbstractHomotopy
    system::S

    start::LinearSubspace{ComplexF64}
    target::LinearSubspace{ComplexF64}
    path::GrassmannianGeodesic

    J::Matrix{ComplexF64}
    x::Vector{ComplexF64}
    ẋ::Vector{ComplexF64}
    x_high::Vector{ComplexDF64}
    t_cache::Base.RefValue{ComplexF64}
    # For AD
    taylor_t_cache::Base.RefValue{ComplexF64}
    taylor_γ::NTuple{5,Matrix{ComplexF64}}
    v::Vector{ComplexF64}
    tx⁴::TaylorVector{5,ComplexF64}
    tx³::TaylorVector{4,ComplexF64}
    tx²::TaylorVector{3,ComplexF64}
    tx¹::TaylorVector{2,ComplexF64}

    affine::Bool
end
is_affine(H::IntrinsicSubspaceStiefelHomotopy) = H.affine
get_b(H::IntrinsicSubspaceStiefelHomotopy) = H.start.intrinsic.b

## Implementation details

# The computation is performed using (implicit) Stiefel coordinates. If `size(F) = m, n`
# and `V` and `W` are affine subspaces of dimension `k, then the computation is performed in
# the affine Grassmanian Graff(k,n) embedded in the Grassmanian(k+1,n+1) using the system
# ``[F(γ(t)v); (γ(t)v)[end] - 1]``. Here the `(γ(t)v)[k+1] - 1` condition ensures that
# we are in the affine Grassmanian and that `(γ(t)v)[1:k]` is the correct value in \C^n.

IntrinsicSubspaceHomotopy(
    F::ModelKit.System,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[]
) = IntrinsicSubspaceStiefelHomotopy(fixed(F; compile = compile), start, target; affine = false)
IntrinsicSubspaceStiefelHomotopy(
    F::ModelKit.System,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    kwargs...
) = IntrinsicSubspaceStiefelHomotopy(fixed(F; compile = compile), start, target; kwargs...)


function IntrinsicSubspaceStiefelHomotopy(
    system::AbstractSystem,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    kwargs...
)
    IntrinsicSubspaceStiefelHomotopy(
        system,
        convert(LinearSubspace{ComplexF64}, start),
        convert(LinearSubspace{ComplexF64}, target);
        kwargs...
    )
end

function IntrinsicSubspaceStiefelHomotopy(
    system::AbstractSystem,
    start::LinearSubspace{ComplexF64},
    target::LinearSubspace{ComplexF64};
    affine::Bool = true
)
    if affine
        LRU = AFFINE_INTRINSIC_LRU
    else
        LRU = PROJECTIVE_INTRINSIC_LRU
    end
    path = get!(LRU, (start, target)) do
        GrassmannianGeodesic(start, target; affine = affine)
    end
    Q = path.Q
    tx⁴ = TaylorVector{5}(ComplexF64, size(Q, 1))
    J = affine ? zeros(ComplexF64, size(system)) : zeros(ComplexF64, size(system) .+ (1, 1))
    IntrinsicSubspaceStiefelHomotopy(
        system = system,
        start = start,
        target = target,
        path = GrassmannianGeodesic(start, target; affine = affine),
        J = J,
        x = zeros(ComplexF64, size(Q, 1)),
        ẋ = zeros(ComplexF64, size(Q, 1)),
        x_high = zeros(ComplexDF64, size(Q, 1)),
        t_cache = Ref(complex(NaN, NaN)),
        taylor_t_cache = Ref(complex(NaN, NaN)),
        taylor_γ = tuple((similar(Q) for i = 0:4)...),
        v = zeros(ComplexF64, size(Q, 2)),
        tx⁴ = tx⁴,
        tx³ = TaylorVector{4}(tx⁴),
        tx² = TaylorVector{3}(tx⁴),
        tx¹ = TaylorVector{2}(tx⁴),
        affine = affine
    )
end
@inline function Base.size(H::IntrinsicSubspaceStiefelHomotopy)
    is_affine(H) ? (size(H.system)[1], dim(H.start)) : (size(H.system)[1] + 1, dim(H.start) + 1)
end



"""
    set_subspaces!(H::IntrinsicSubspaceStiefelHomotopy, start::LinearSubspace, target::LinearSubspace)

Update the homotopy `H` to track from the affine subspace `start` to `target`.
"""
function set_subspaces!(
    H::IntrinsicSubspaceStiefelHomotopy,
    start::LinearSubspace,
    target::LinearSubspace,
)
    H.start = start
    H.target = target
    if is_affine(H)
        LRU = AFFINE_INTRINSIC_LRU
    else
        LRU = PROJECTIVE_INTRINSIC_LRU
    end

    H.path = get!(LRU, (start, target)) do
        GrassmannianGeodesic(start, target)
    end
    H
end
start_parameters!(H::IntrinsicSubspaceStiefelHomotopy, p::LinearSubspace) =
    set_subspaces!(H, convert(LinearSubspace{ComplexF64}, p), H.target)
target_parameters!(H::IntrinsicSubspaceStiefelHomotopy, q::LinearSubspace) =
    set_subspaces!(H, H.start, convert(LinearSubspace{ComplexF64}, q))
parameters!(H::IntrinsicSubspaceStiefelHomotopy, p::LinearSubspace, q::LinearSubspace) =
    set_subspaces!(
        H,
        convert(LinearSubspace{ComplexF64}, p),
        convert(LinearSubspace{ComplexF64}, q),
    )

function γ!(H::IntrinsicSubspaceStiefelHomotopy, t::Number)
    H.t_cache[] != t || return first(H.taylor_γ)
    if isreal(t)
        _γ!(H, real(t))
    else
        _γ!(H, t)
    end
    H.t_cache[] = t

    first(H.taylor_γ)
end
@inline function _γ!(H::IntrinsicSubspaceStiefelHomotopy, t::Number)
    @unpack Q, Q_cos, Θ = H.path
    γ = first(H.taylor_γ)
    n, k = size(γ)
    @inbounds for j = 1:k
        Θⱼ = Θ[j]
        s, c = sincos(t * Θⱼ)
        for i = 1:n
            γ[i, j] = Q_cos[i, j] * c + Q[i, j] * s
        end
    end
    γ
end

γ̇!(H::IntrinsicSubspaceStiefelHomotopy, t::Number) = isreal(t) ? _γ̇!(H, real(t)) : _γ̇!(H, t)
@inline function _γ̇!(H::IntrinsicSubspaceStiefelHomotopy, t::Number)
    @unpack Q, Q_cos, Θ = H.path
    _, γ̇ = H.taylor_γ
    n, k = size(γ̇)
    @inbounds for j = 1:k
        Θⱼ = Θ[j]
        s, c = sincos(t * Θⱼ)
        ċ = -s * Θⱼ
        ṡ = c * Θⱼ
        for i = 1:n
            γ̇[i, j] = Q_cos[i, j] * ċ + Q[i, j] * ṡ
        end
    end
    γ̇
end

function set_solution!(u::Vector, H::IntrinsicSubspaceStiefelHomotopy, x::AbstractVector, t)

    if is_affine(H)
        (length(x) == length(H.x)) ||
        throw(ArgumentError("Cannot set solution. Expected extrinsic coordinates."))

        set_solution!(H.x, H.system, x)
        b = get_b(H)
        H.x .= H.x - b
    else
        (length(x) == length(H.x) - 1) ||
            throw(ArgumentError("Cannot set solution. Expected extrinsic coordinates."))

        set_solution!(view(H.x, 1:length(x)), H.system, x)
        H.x[end] = 1
    end

    if isone(t)
        LA.mul!(u, H.path.γ1', H.x)
    elseif iszero(t)
        LA.mul!(u, H.path.Q_cos', H.x)
    else
        LA.mul!(u, γ!(H, t)', H.x)
    end
end

function get_solution(H::IntrinsicSubspaceStiefelHomotopy, u::AbstractVector, t)

    if is_affine(H)
        b = get_b(H)
        if isone(t)
            out = (H.path.γ1 * u) + b
        elseif iszero(t)
            out = (H.path.Q_cos * u) + b
        else
            γ = γ!(H, t)
            out = (γ * u) + b    
        end
    else
        if isone(t)
            out = (@view H.path.γ1[1:end-1, :]) * u
        elseif iszero(t)
            out = (@view H.path.Q_cos[1:end-1, :]) * u
        else
            γ = γ!(H, t)
            out = (@view γ[1:end-1, :]) * u
        end
    end

    out
end

function ModelKit.evaluate!(u, H::IntrinsicSubspaceStiefelHomotopy, v::Vector{ComplexDF64}, t)
    γ = γ!(H, t)
    LA.mul!(H.x_high, γ, v)

    if is_affine(H)
        b = get_b(H)
        H.x_high .= H.x_high .+ b
    else
        n = first(size(H.system))
        u[n+1] = H.x_high[end] - 1.0
    end
    evaluate!(u, H.system, H.x_high)
   
    u
end

function ModelKit.evaluate!(u, H::IntrinsicSubspaceStiefelHomotopy, v::AbstractVector, t)
    γ = γ!(H, t)
    LA.mul!(H.x, γ, v)

    if is_affine(H)
        b = get_b(H)
        H.x .= H.x .+ b
    else
        n = first(size(H.system))
        u[n+1] = H.x[end] - 1.0
    end
    evaluate!(u, H.system, H.x)
   
    u
end

function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::IntrinsicSubspaceStiefelHomotopy,
    v::AbstractVector,
    t,
)
    γ = γ!(H, t)
    LA.mul!(H.x, γ, v)

    if is_affine(H)
        b = get_b(H)
        H.x .= H.x .+ b
    end
    evaluate_and_jacobian!(u, H.J, H.system, H.x)
    LA.mul!(U, H.J, γ)
    
    if !is_affine(H)
        n = first(size(H.system))
        u[n+1] = H.x[end] - 1
        m = length(v)
        for j = 1:m
            U[n+1, j] = γ[end, j]
        end
    end

    nothing
end

function ModelKit.taylor!(
    u,
    ::Val{1},
    H::IntrinsicSubspaceStiefelHomotopy,
    v,
    t,
    incr::Bool = true,
)
    γ = γ!(H, t)
    γ̇ = γ̇!(H, t)

    # apply chain rule
    #    d/dt F(γ(t)v) = J_F(γ(t)v)* γ̇(t)*v

    H.v .= first.(v)
    LA.mul!(H.x, γ, H.v)
    LA.mul!(H.ẋ, γ̇, H.v)

    if is_affine(H)
        b = get_b(H)
        H.x .= H.x .+ b
    end
    evaluate_and_jacobian!(u, H.J, H.system, H.x)
    LA.mul!(u, H.J, H.ẋ)

    if !is_affine(H)
        M = size(H, 1)
        u[M] = H.ẋ[end]
    end

    u
end

function _taylor_γ!(H::IntrinsicSubspaceStiefelHomotopy, t::Number)
    @unpack path, taylor_γ = H
    @unpack Q, Q_cos, Θ, U = path

    γ, γ¹, γ², γ³, γ⁴ = taylor_γ
    n, k = size(γ)
    @inbounds for j = 1:k
        Θⱼ = Θ[j]
        s, c = sincos(t * Θⱼ)
        c¹ = -s * Θⱼ
        s¹ = c * Θⱼ
        Θⱼ_2 = 0.5 * Θⱼ^2
        c² = -c * Θⱼ_2
        s² = -s * Θⱼ_2
        Θⱼ_3 = Θⱼ_2 * Θⱼ / 3
        c³ = s * Θⱼ_3
        s³ = -c * Θⱼ_3
        Θⱼ_4 = 0.25 * Θⱼ_3 * Θⱼ
        c⁴ = c * Θⱼ_4
        s⁴ = s * Θⱼ_4
        for i = 1:n
            γ[i, j] = Q_cos[i, j] * c + Q[i, j] * s
            γ¹[i, j] = Q_cos[i, j] * c¹ + Q[i, j] * s¹
            γ²[i, j] = Q_cos[i, j] * c² + Q[i, j] * s²
            γ³[i, j] = Q_cos[i, j] * c³ + Q[i, j] * s³
            γ⁴[i, j] = Q_cos[i, j] * c⁴ + Q[i, j] * s⁴
        end
    end

end

function taylor_γ!(H::IntrinsicSubspaceStiefelHomotopy, t::Number)
    H.taylor_t_cache[] != t || return H.taylor_γ

    if isreal(t)
        _taylor_γ!(H, real(t))
    else
        _taylor_γ!(H, t)
    end
    H.taylor_t_cache[] = t

    H.taylor_γ
end


function ModelKit.taylor!(
    u,
    ::Val{2},
    H::IntrinsicSubspaceStiefelHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ² = taylor_γ!(H, t)
    x, x¹, x² = vectors(H.tx²)
    v, v¹ = vectors(tv)

    if incr
        x .= H.x
        x¹ .= H.ẋ
    else
        LA.mul!(x, γ, v)
        LA.mul!(x¹, γ¹, v)
    end

    H.v .= v¹
    LA.mul!(H.x, γ, v¹)
    x¹ .+= H.x

    LA.mul!(H.x, γ¹, H.v)
    H.v .= v
    LA.mul!(H.x, γ², H.v, true, true)
    x² .= H.x

    taylor!(u, Val(2), H.system, H.tx²)

    if !is_affine(H)
        n = first(size(H.system))
        u[n+1] = x²[end]
    end

    u
end

function ModelKit.taylor!(
    u,
    ::Val{3},
    H::IntrinsicSubspaceStiefelHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ², γ³ = taylor_γ!(H, t)
    x, x¹, x², x³ = vectors(H.tx³)
    v, v¹, v² = vectors(tv)

    if !incr
        LA.mul!(x, γ, v)
        LA.mul!(x¹, γ¹, v)
        LA.mul!(x¹, γ, v¹, true, true)
        LA.mul!(x², γ², v)
        LA.mul!(x², γ¹, v¹, true, true)
    end

    H.v .= v²
    LA.mul!(H.x, γ, H.v)
    x² .+= H.x

    LA.mul!(H.x, γ¹, H.v)
    H.v .= v
    LA.mul!(H.x, γ³, H.v, true, true)
    H.v .= v¹
    LA.mul!(H.x, γ², H.v, true, true)
    x³ .= H.x

    taylor!(u, Val(3), H.system, H.tx³)

    if !is_affine(H)
        n = first(size(H.system))
        u[n+1] = x³[end]
    end

    u
end

function ModelKit.taylor!(
    u,
    ::Val{4},
    H::IntrinsicSubspaceStiefelHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ², γ³, γ⁴ = taylor_γ!(H, t)
    x, x¹, x², x³, x⁴ = vectors(H.tx⁴)
    v, v¹, v², v³ = vectors(tv)

    if !incr
        LA.mul!(x, γ, v)
        LA.mul!(x¹, γ¹, v)
        LA.mul!(x¹, γ, v¹, true, true)
        LA.mul!(x², γ², v)
        LA.mul!(x², γ¹, v¹, true, true)
        LA.mul!(x², γ, v², true, true)

        LA.mul!(x³, γ³, v)
        LA.mul!(x³, γ², v¹, true, true)
        LA.mul!(x³, γ¹, v², true, true)
    end

    LA.mul!(x³, γ, v³, true, true)

    # Currently H.v .= v³
    LA.mul!(x⁴, γ¹, v³)
    LA.mul!(x⁴, γ², v², true, true)
    LA.mul!(x⁴, γ³, v¹, true, true)
    LA.mul!(x⁴, γ⁴, v, true, true)

    taylor!(u, Val(4), H.system, H.tx⁴)

    if !is_affine(H)
        n = first(size(H.system))
        u[n+1] = x⁴[end]
    end

    u
end


"""
    IntrinsicSubspaceOffsetHomotopy(F::System, V::LinearSubspace, W::LinearSubspace)
    IntrinsicSubspaceOffsetHomotopy(F::AbstractSystem, V::LinearSubspace, W::LinearSubspace)

Creates a homotopy ``H(x,t) = F(Bx + ta + (1-t)b)`` where the matrix part `B` is kept fixed
and the offset interpolates linearly from `a` (start) to `b` (target).

The affine subspaces `V` and `W` represent the start and target configurations:
- `V` represents the subspace `Bx + a`
- `W` represents the subspace `Bx + b`

At ``t=1``, we have ``H(x,1) = F(Bx + a)`` (start configuration).
At ``t=0``, we have ``H(x,0) = F(Bx + b)`` (target configuration).

Unlike `IntrinsicSubspaceStiefelHomotopy` which moves the matrix part `A` to `B`,
this homotopy moves the offset part from `a` to `b` while keeping the matrix fixed.

See also [`LinearSubspace`](@ref).
"""
Base.@kwdef mutable struct IntrinsicSubspaceOffsetHomotopy{S<:AbstractSystem} <: AbstractHomotopy
    system::S

    start::LinearSubspace{ComplexF64}
    target::LinearSubspace{ComplexF64}

    a_minus_b::Vector{ComplexF64}

    J::Matrix{ComplexF64}
    x::Vector{ComplexF64}
    x_high::Vector{ComplexDF64}
    v::Vector{ComplexF64}
    # For AD
    tx⁴::TaylorVector{5,ComplexF64}
    tx³::TaylorVector{4,ComplexF64}
    tx²::TaylorVector{3,ComplexF64}
    tx¹::TaylorVector{2,ComplexF64}
end

IntrinsicSubspaceOffsetHomotopy(
    F::ModelKit.System,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
) = IntrinsicSubspaceOffsetHomotopy(fixed(F; compile = compile), start, target)


function IntrinsicSubspaceOffsetHomotopy(
    system::AbstractSystem,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[]
)
    IntrinsicSubspaceOffsetHomotopy(
        system,
        convert(LinearSubspace{ComplexF64}, start),
        convert(LinearSubspace{ComplexF64}, target);
    )
end

function IntrinsicSubspaceOffsetHomotopy(
    system::AbstractSystem,
    start::LinearSubspace{ComplexF64},
    target::LinearSubspace{ComplexF64}
)
    X = target.intrinsic.X
    n = size(X, 1)
    J = zeros(ComplexF64, size(system)) 
    tx⁴ = TaylorVector{5}(ComplexF64, n)
    a = start.intrinsic.b
    b = target.intrinsic.b
    
    IntrinsicSubspaceOffsetHomotopy(
        system = system,
        start = start,
        target = target,
        a_minus_b = a - b,
        J = J,
        x = zeros(ComplexF64, n),
        x_high = zeros(ComplexDF64, n),
        v = zeros(ComplexF64, size(X, 2)),
        tx⁴ = tx⁴,
        tx³ = TaylorVector{4}(tx⁴),
        tx² = TaylorVector{3}(tx⁴),
        tx¹ = TaylorVector{2}(tx⁴)
    )
end

Base.size(H::IntrinsicSubspaceOffsetHomotopy) = (size(H.system)[1], dim(H.start)) 
get_X(H::IntrinsicSubspaceOffsetHomotopy) = H.target.intrinsic.X

"""
    set_subspaces!(H::IntrinsicSubspaceOffsetHomotopy, start::LinearSubspace, target::LinearSubspace)

Update the homotopy `H` to track from the affine subspace `start` to `target`,
updating both the matrix part `B` and the offsets `a` and `b`.
"""
function set_subspaces!(
    H::IntrinsicSubspaceOffsetHomotopy,
    start::LinearSubspace,
    target::LinearSubspace,
)
    H.start = start
    H.target = target

    H
end

start_parameters!(H::IntrinsicSubspaceOffsetHomotopy, p::LinearSubspace) =
    set_subspaces!(H, convert(LinearSubspace{ComplexF64}, p), H.target)
target_parameters!(H::IntrinsicSubspaceOffsetHomotopy, q::LinearSubspace) =
    set_subspaces!(H, H.start, convert(LinearSubspace{ComplexF64}, q))
parameters!(H::IntrinsicSubspaceOffsetHomotopy, p::LinearSubspace, q::LinearSubspace) =
    set_subspaces!(
        H,
        convert(LinearSubspace{ComplexF64}, p),
        convert(LinearSubspace{ComplexF64}, q),
    )

# Compute the offset at parameter t: offset(t) = t*a + (1-t)*b = t(a-b) + b
offset_at_t(H::IntrinsicSubspaceOffsetHomotopy, t::Number) = t .* H.a_minus_b .+ H.target.intrinsic.b

function set_solution!(u::Vector, H::IntrinsicSubspaceOffsetHomotopy, x::AbstractVector, t)
    (length(x) == length(H.x)) ||
        throw(ArgumentError("Cannot set solution. Expected extrinsic coordinates."))

    set_solution!(H.x, H.system, x)
    b = offset_at_t(H, t)
    H.x .= H.x - b
    X = get_X(H)
    LA.mul!(u, X', H.x)

    nothing
end

function get_solution(H::IntrinsicSubspaceOffsetHomotopy, u::AbstractVector, t)
    # u are coordinates in the subspace
    # Return x in the ambient space via x = X*u + b(t)
    b = offset_at_t(H, t)
    X = get_X(H)
    out = X * u + b
    
    out
end

function ModelKit.evaluate!(u, H::IntrinsicSubspaceOffsetHomotopy, v::Vector{ComplexDF64}, t)
    b = offset_at_t(H, t)
    X = get_X(H)

    LA.mul!(H.x_high, X, v)
    H.x_high .+= b

    evaluate!(u, H.system, H.x_high)

    u
end

function ModelKit.evaluate!(u, H::IntrinsicSubspaceOffsetHomotopy, v::AbstractVector, t)
    b = offset_at_t(H, t)
    X = get_X(H)

    LA.mul!(H.x, X, v)
    H.x .+= b

    evaluate!(u, H.system, H.x)

    u
end

function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::IntrinsicSubspaceOffsetHomotopy,
    v::AbstractVector,
    t,
)
    b = offset_at_t(H, t)
    X = get_X(H)

    LA.mul!(H.x, X, v)
    H.x .+= b
    evaluate_and_jacobian!(u, H.J, H.system, H.x)
    LA.mul!(U, H.J, X)

    nothing
end

function ModelKit.taylor!(
    u,
    ::Val{1},
    H::IntrinsicSubspaceOffsetHomotopy,
    v,
    t,
    incr::Bool = true,
)
    # d/dt H(v,t) = d/dt F(X*v + t*a + (1-t)*b)
    #             = J_F(X*v + t*a + (1-t)*b) * (a - b)

    H.v .= first.(v)
    b = offset_at_t(H, t)
    X = get_X(H)

    LA.mul!(H.x, X, H.v)
    H.x .= H.x .+ b

    evaluate_and_jacobian!(u, H.J, H.system, H.x)

    # u = J_F * (a - b)
    LA.mul!(u, H.J, H.a_minus_b)

    u
end

function ModelKit.taylor!(
    u,
    ::Val{2},
    H::IntrinsicSubspaceOffsetHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    x, x¹, x² = vectors(H.tx²)
    v, v¹ = vectors(tv)
    X = get_X(H)

    if incr
        x .= H.x
    else
        b = offset_at_t(H, t)
       
        LA.mul!(x, X, v)
        x .+= b
    end

    # x¹ = B * v¹ + d/dt[offset]
    H.v .= v¹
    LA.mul!(x¹, X, H.v)
    x¹ .+= H.a_minus_b  # derivative of offset with respect to t

    # x² = B * v
    H.v .= v
    LA.mul!(x², X, H.v)

    taylor!(u, Val(2), H.system, H.tx²)

    u
end

function ModelKit.taylor!(
    u,
    ::Val{3},
    H::IntrinsicSubspaceOffsetHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    x, x¹, x², x³ = vectors(H.tx³)
    v, v¹, v² = vectors(tv)
    X = get_X(H)

    if !incr
        b = offset_at_t(H, t)    
        LA.mul!(x, X, v)
        x .+= b
        LA.mul!(x¹, X, v¹)
        x¹ .+= H.a_minus_b
    end

    # x² = B * v²
    H.v .= v²
    LA.mul!(x², X, H.v)

    # x³ = B * v
    H.v .= v
    LA.mul!(x³, X, H.v)

    taylor!(u, Val(3), H.system, H.tx³)

    u
end

function ModelKit.taylor!(
    u,
    ::Val{4},
    H::IntrinsicSubspaceOffsetHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    x, x¹, x², x³, x⁴ = vectors(H.tx⁴)
    v, v¹, v², v³ = vectors(tv)
    X = get_X(H)

    if !incr
        b = offset_at_t(H, t)
        LA.mul!(x, X, v)
        x .+= b
        LA.mul!(x¹, X, v¹)
        x¹ .+= H.a_minus_b
    end

    # Higher order derivatives
    H.v .= v²
    LA.mul!(x², X, H.v)

    H.v .= v³
    LA.mul!(x³, X, H.v)

    H.v .= v
    LA.mul!(x⁴, X, H.v)

    taylor!(u, Val(4), H.system, H.tx⁴)

    u
end

