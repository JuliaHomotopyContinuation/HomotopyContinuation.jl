export IntrinsicSubspaceCombinedHomotopy

"""
    IntrinsicSubspaceCombinedHomotopy(F::System, V::LinearSubspace, W::LinearSubspace)
    IntrinsicSubspaceCombinedHomotopy(F::AbstractSystem, V::LinearSubspace, W::LinearSubspace)

Creates a homotopy ``H(x,t) = F(γ(t)x + b(t))`` where both the matrix part and the offset
are interpolated simultaneously.

At ``t=1``, we have ``H(x,1) = F(Ax + a)`` where ``V = {Ax + a}``.
At ``t=0``, we have ``H(x,0) = F(Bx + b)`` where ``W = {Bx + b}``.

The matrix part ``γ(t)`` follows a geodesic in the affine Grassmannian from ``A`` to ``B``,
while the offset ``b(t)`` interpolates linearly from ``a`` to ``b``.

See also [`LinearSubspace`](@ref), [`IntrinsicSubspaceStiefelHomotopy`](@ref), 
and [`IntrinsicSubspaceOffsetHomotopy`](@ref).
"""
Base.@kwdef mutable struct IntrinsicSubspaceCombinedHomotopy{S<:AbstractSystem} <: AbstractHomotopy
    system::S

    start::LinearSubspace{ComplexF64}
    target::LinearSubspace{ComplexF64}
    
    # For the matrix part (geodesic interpolation)
    path::GrassmannianGeodesic

    # For the offset part (linear interpolation)
    a_minus_b::Vector{ComplexF64}
    offset::Vector{ComplexF64}
    t_cache::Base.RefValue{ComplexF64}
    offset_t_cache::Base.RefValue{ComplexF64}

    J::Matrix{ComplexF64}
    x::Vector{ComplexF64}
    ẋ::Vector{ComplexF64}
    x_high::Vector{ComplexDF64}
    # For AD
    taylor_t_cache::Base.RefValue{ComplexF64}
    taylor_γ::NTuple{5,Matrix{ComplexF64}}
    v::Vector{ComplexF64}
    tx⁴::TaylorVector{5,ComplexF64}
    tx³::TaylorVector{4,ComplexF64}
    tx²::TaylorVector{3,ComplexF64}
    tx¹::TaylorVector{2,ComplexF64}
end

IntrinsicSubspaceCombinedHomotopy(
    F::ModelKit.System,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[]
) = IntrinsicSubspaceCombinedHomotopy(fixed(F; compile = compile), start, target)

function IntrinsicSubspaceCombinedHomotopy(
    system::AbstractSystem,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[]
)
    IntrinsicSubspaceCombinedHomotopy(
        system,
        convert(LinearSubspace{ComplexF64}, start),
        convert(LinearSubspace{ComplexF64}, target)
    )
end

function IntrinsicSubspaceCombinedHomotopy(
    system::AbstractSystem,
    start::LinearSubspace{ComplexF64},
    target::LinearSubspace{ComplexF64}
)
    # Create geodesic path for the matrix part
    path = GrassmannianGeodesic(start, target; affine = true)
    Q = path.Q
    
    # Prepare offset data for linear interpolation
    a = start.intrinsic.b
    b = target.intrinsic.b
    
    tx⁴ = TaylorVector{5}(ComplexF64, size(Q, 1))
    J = zeros(ComplexF64, size(system))
    
    IntrinsicSubspaceCombinedHomotopy(
        system = system,
        start = start,
        target = target,
        path = path,
        a_minus_b = a - b,
        offset = copy(b),
        t_cache = Ref(complex(NaN, NaN)),
        offset_t_cache = Ref(complex(NaN, NaN)),
        J = J,
        x = zeros(ComplexF64, size(Q, 1)),
        ẋ = zeros(ComplexF64, size(Q, 1)),
        x_high = zeros(ComplexDF64, size(Q, 1)),
        taylor_t_cache = Ref(complex(NaN, NaN)),
        taylor_γ = tuple((similar(Q) for i = 0:4)...),
        v = zeros(ComplexF64, size(Q, 2)),
        tx⁴ = tx⁴,
        tx³ = TaylorVector{4}(tx⁴),
        tx² = TaylorVector{3}(tx⁴),
        tx¹ = TaylorVector{2}(tx⁴)
    )
end

@inline function Base.size(H::IntrinsicSubspaceCombinedHomotopy)
    (size(H.system)[1], dim(H.start))
end

"""
    set_subspaces!(H::IntrinsicSubspaceCombinedHomotopy, start::LinearSubspace, target::LinearSubspace)

Update the homotopy `H` to track from the affine subspace `start` to `target`.
"""
function set_subspaces!(
    H::IntrinsicSubspaceCombinedHomotopy,
    start::LinearSubspace,
    target::LinearSubspace,
)
    H.start = start
    H.target = target
    
    # Update the geodesic path
    H.path = GrassmannianGeodesic(start, target; affine = true)
    
    # Update the offsets
    a = start.intrinsic.b
    b = target.intrinsic.b
    H.a_minus_b .= a .- b
    H.offset .= b
    H.t_cache[] = complex(NaN, NaN)
    H.offset_t_cache[] = complex(NaN, NaN)
    H.taylor_t_cache[] = complex(NaN, NaN)
    
    H
end

start_parameters!(H::IntrinsicSubspaceCombinedHomotopy, p::LinearSubspace) =
    set_subspaces!(H, convert(LinearSubspace{ComplexF64}, p), H.target)
target_parameters!(H::IntrinsicSubspaceCombinedHomotopy, q::LinearSubspace) =
    set_subspaces!(H, H.start, convert(LinearSubspace{ComplexF64}, q))
parameters!(H::IntrinsicSubspaceCombinedHomotopy, p::LinearSubspace, q::LinearSubspace) =
    set_subspaces!(
        H,
        convert(LinearSubspace{ComplexF64}, p),
        convert(LinearSubspace{ComplexF64}, q),
    )

# Helper function for computing γ at parameter t
function γ!(H::IntrinsicSubspaceCombinedHomotopy, t::Number)
    H.t_cache[] != t || return first(H.taylor_γ)
    if isreal(t)
        _γ!(H, real(t))
    else
        _γ!(H, t)
    end
    H.t_cache[] = t

    first(H.taylor_γ)
end

@inline function _γ!(H::IntrinsicSubspaceCombinedHomotopy, t::Number)
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

γ̇!(H::IntrinsicSubspaceCombinedHomotopy, t::Number) = isreal(t) ? _γ̇!(H, real(t)) : _γ̇!(H, t)

@inline function _γ̇!(H::IntrinsicSubspaceCombinedHomotopy, t::Number)
    @unpack Q, Q_cos, Θ = H.path
    _, γ̇ = H.taylor_γ
    n, k = size(γ̇)
    @inbounds for j = 1:k
        Θⱼ = Θ[j]
        s, c = sincos(t * Θⱼ)
        ċ = -s * Θⱼ
        ṡ = c * Θⱼ
        for i = 1:n
            γ̇[i, j] = Q_cos[i, j] * ċ + Q[i, j] * ṡ
        end
    end
    γ̇
end

# Compute the offset at parameter t: offset(t) = t*a + (1-t)*b = t(a-b) + b
@inline function _compute_offset!(H::IntrinsicSubspaceCombinedHomotopy, t::Number)
    H.offset .= H.target.intrinsic.b
    LA.axpy!(t, H.a_minus_b, H.offset)
    H.offset_t_cache[] = t
end

@inline function offset_at_t!(H::IntrinsicSubspaceCombinedHomotopy, t::Number)
    if H.offset_t_cache[] != t
        _compute_offset!(H, t)
    end
    H.offset
end

function set_solution!(u::Vector, H::IntrinsicSubspaceCombinedHomotopy, x::AbstractVector, t)
    (length(x) == length(H.x)) ||
        throw(ArgumentError("Cannot set solution. Expected extrinsic coordinates."))

    set_solution!(H.x, H.system, x)
    offset_at_t!(H, t)
    H.x .= H.x .- H.offset
    
    if isone(t)
        LA.mul!(u, H.path.γ1', H.x)
    elseif iszero(t)
        LA.mul!(u, H.path.Q_cos', H.x)
    else
        LA.mul!(u, γ!(H, t)', H.x)
    end
end

function get_solution(H::IntrinsicSubspaceCombinedHomotopy, u::AbstractVector, t)
    offset_at_t!(H, t)
    
    if isone(t)
        out = (H.path.γ1 * u) + H.offset
    elseif iszero(t)
        out = (H.path.Q_cos * u) + H.offset
    else
        γ = γ!(H, t)
        out = (γ * u) + H.offset
    end
    
    out
end

function ModelKit.evaluate!(u, H::IntrinsicSubspaceCombinedHomotopy, v::Vector{ComplexDF64}, t)
    γ = γ!(H, t)
    offset_at_t!(H, t)
    
    LA.mul!(H.x_high, γ, v)
    LA.axpy!(1, H.offset, H.x_high)

    evaluate!(u, H.system, H.x_high)
   
    u
end

function ModelKit.evaluate!(u, H::IntrinsicSubspaceCombinedHomotopy, v::AbstractVector, t)
    γ = γ!(H, t)
    offset_at_t!(H, t)

    LA.mul!(H.x, γ, v)
    LA.axpy!(1, H.offset, H.x)

    evaluate!(u, H.system, H.x)
   
    u
end

function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::IntrinsicSubspaceCombinedHomotopy,
    v::AbstractVector,
    t,
)
    γ = γ!(H, t)
    offset_at_t!(H, t)
    
    LA.mul!(H.x, γ, v)
    LA.axpy!(1, H.offset, H.x)
    
    evaluate_and_jacobian!(u, H.J, H.system, H.x)
    LA.mul!(U, H.J, γ)

    nothing
end

function ModelKit.taylor!(
    u,
    ::Val{1},
    H::IntrinsicSubspaceCombinedHomotopy,
    v,
    t,
    incr::Bool = true,
)
    # d/dt H(v,t) = d/dt F(γ(t)v + b(t))
    #             = J_F(γ(t)v + b(t)) * (γ̇(t)v + db/dt)

    γ = γ!(H, t)
    γ̇ = γ̇!(H, t)

    H.v .= first.(v)
    LA.mul!(H.x, γ, H.v)
    offset_at_t!(H, t)
    LA.axpy!(1, H.offset, H.x)
    
    LA.mul!(H.ẋ, γ̇, H.v)
    LA.axpy!(1, H.a_minus_b, H.ẋ)  # derivative of offset with respect to t
    
    evaluate_and_jacobian!(u, H.J, H.system, H.x)
    LA.mul!(u, H.J, H.ẋ)

    u
end

function ModelKit.taylor!(
    u,
    ::Val{2},
    H::IntrinsicSubspaceCombinedHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ² = taylor_γ!(H, t)
    x, x¹, x² = vectors(H.tx²)
    v, v¹ = vectors(tv)

    if incr
        x .= H.x
        x¹ .= H.ẋ
    else
        offset_at_t!(H, t)
        LA.mul!(x, γ, v)
        LA.axpy!(1, H.offset, x)
        
        LA.mul!(x¹, γ¹, v)
        LA.axpy!(1, H.a_minus_b, x¹)
    end

    H.v .= v¹
    LA.mul!(H.x, γ, H.v)
    x¹ .+= H.x

    LA.mul!(H.x, γ¹, H.v)
    H.v .= v
    LA.mul!(H.x, γ², H.v, true, true)
    x² .= H.x

    taylor!(u, Val(2), H.system, H.tx²)

    u
end

function ModelKit.taylor!(
    u,
    ::Val{3},
    H::IntrinsicSubspaceCombinedHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ², γ³ = taylor_γ!(H, t)
    x, x¹, x², x³ = vectors(H.tx³)
    v, v¹, v² = vectors(tv)

    if !incr
        offset_at_t!(H, t)
        LA.mul!(x, γ, v)
        LA.axpy!(1, H.offset, x)
        LA.mul!(x¹, γ¹, v)
        LA.axpy!(1, H.a_minus_b, x¹)
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

    u
end

function ModelKit.taylor!(
    u,
    ::Val{4},
    H::IntrinsicSubspaceCombinedHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹, γ², γ³, γ⁴ = taylor_γ!(H, t)
    x, x¹, x², x³, x⁴ = vectors(H.tx⁴)
    v, v¹, v², v³ = vectors(tv)

    if !incr
        offset_at_t!(H, t)
        LA.mul!(x, γ, v)
        LA.axpy!(1, H.offset, x)
        LA.mul!(x¹, γ¹, v)
        LA.axpy!(1, H.a_minus_b, x¹)
        LA.mul!(x¹, γ, v¹, true, true)
        LA.mul!(x², γ², v)
        LA.mul!(x², γ¹, v¹, true, true)
        LA.mul!(x², γ, v², true, true)

        LA.mul!(x³, γ³, v)
        LA.mul!(x³, γ², v¹, true, true)
        LA.mul!(x³, γ¹, v², true, true)
    end

    LA.mul!(x³, γ, v³, true, true)

    LA.mul!(x⁴, γ¹, v³)
    LA.mul!(x⁴, γ², v², true, true)
    LA.mul!(x⁴, γ³, v¹, true, true)
    LA.mul!(x⁴, γ⁴, v, true, true)

    taylor!(u, Val(4), H.system, H.tx⁴)

    u
end

function _taylor_γ!(H::IntrinsicSubspaceCombinedHomotopy, t::Number)
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

function taylor_γ!(H::IntrinsicSubspaceCombinedHomotopy, t::Number)
    H.taylor_t_cache[] != t || return H.taylor_γ

    if isreal(t)
        _taylor_γ!(H, real(t))
    else
        _taylor_γ!(H, t)
    end
    H.taylor_t_cache[] = t

    H.taylor_γ
end
