export IntrinsicSubspaceHomotopy, set_subspaces!


Base.@kwdef struct SubspaceHomotopy <: ModelKit.AbstractHomotopy
    start::LinearSubspace{ComplexF64}
    target::LinearSubspace{ComplexF64}
    extrinsic::Bool
    A::Matrix{ComplexF64}
    b::Vector{ComplexF64}
    ū::Vector{ComplexDF64}
end


"""
    IntrinsicSubspaceHomotopy(F::System, V::LinearSubspace, W::LinearSubspace)
    IntrinsicSubspaceHomotopy(F::AbstractSystem, V::LinearSubspace, W::LinearSubspace)

Creates a homotopy ``H(x,t) = F(γ(t)x)`` where ``γ(t)`` is a family of affine subspaces
such that ``γ(1) = V`` and ``γ(0) = W``.
Here ``γ(t)`` is the geodesic between `V` and `W` in the affine Grassmanian, i.e.,
it is the curve of minimal length connecting `V` and `W`.
See also [`LinearSubspace`](@ref) and [`geodesic`](@ref) and the references therein.
"""
Base.@kwdef mutable struct IntrinsicSubspaceHomotopy{S<:AbstractSystem} <: AbstractHomotopy
    system::S

    start::LinearSubspace{ComplexF64}
    target::LinearSubspace{ComplexF64}

    J::Matrix{ComplexF64}
    x::Vector{ComplexF64}
    ẋ::Vector{ComplexF64}
    x_high::Vector{ComplexDF64}
    t_cache::Base.RefValue{ComplexF64}
    # For AD
    taylor_t_cache::Base.RefValue{ComplexF64}
    A::Matrix{ComplexF64}
    b::Vector{ComplexF64}
    Ȧ::Matrix{ComplexF64}
    ḃ::Vector{ComplexF64}
    v::Vector{ComplexF64}
    tx⁴::TaylorVector{5,ComplexF64}
    tx³::TaylorVector{4,ComplexF64}
    tx²::TaylorVector{3,ComplexF64}
    tx¹::TaylorVector{2,ComplexF64}
end

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
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
) = IntrinsicSubspaceHomotopy(fixed(F; compile = compile), start, target)


function IntrinsicSubspaceHomotopy(
    system::AbstractSystem,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
)
    IntrinsicSubspaceHomotopy(
        system,
        convert(LinearSubspace{ComplexF64}, start),
        convert(LinearSubspace{ComplexF64}, target),
    )
end

function IntrinsicSubspaceHomotopy(
    system::AbstractSystem,
    start::LinearSubspace{ComplexF64},
    target::LinearSubspace{ComplexF64},
)
    A = copy(intrinsic(start).A)
    b = copy(intrinsic(start).b)

    Ȧ = intrinsic(start).A - intrinsic(target).A
    ḃ = intrinsic(start).b - intrinsic(target).b

    tx⁴ = TaylorVector{5}(ComplexF64, size(Q, 1))
    IntrinsicSubspaceHomotopy(
        system = system,
        start = start,
        target = target,
        A = A,
        b = b,
        Ȧ = Ȧ,
        ḃ = ḃ,
        J = zeros(ComplexF64, size(system) .+ (1, 1)),
        x = zeros(ComplexF64, size(Q, 1)),
        ẋ = zeros(ComplexF64, size(Q, 1)),
        x_high = zeros(ComplexDF64, size(Q, 1)),
        t_cache = Ref(complex(NaN, NaN)),
        taylor_t_cache = Ref(complex(NaN, NaN)),
        v = zeros(ComplexF64, size(Q, 2)),
        tx⁴ = tx⁴,
        tx³ = TaylorVector{4}(tx⁴),
        tx² = TaylorVector{3}(tx⁴),
        tx¹ = TaylorVector{2}(tx⁴),
    )
end
Base.size(H::IntrinsicSubspaceHomotopy) = (size(H.system)[1], dim(H.start))

"""
    set_subspaces!(H::IntrinsicSubspaceHomotopy, start::LinearSubspace, target::LinearSubspace)

Update the homotopy `H` to track from the affine subspace `start` to `target`.
"""
function set_subspaces!(
    H::IntrinsicSubspaceHomotopy,
    start::LinearSubspace,
    target::LinearSubspace,
)
    H.start = start
    H.target = target
    H.Ȧ .= intrinsic(start).A .- intrinsic(target).A
    H.ḃ .= intrinsic(start).b .- intrinsic(target).b
    H
end
start_parameters!(H::IntrinsicSubspaceHomotopy, p::LinearSubspace) =
    set_subspaces!(H, convert(LinearSubspace{ComplexF64}, p), H.target)
target_parameters!(H::IntrinsicSubspaceHomotopy, q::LinearSubspace) =
    set_subspaces!(H, H.start, convert(LinearSubspace{ComplexF64}, q))
parameters!(H::IntrinsicSubspaceHomotopy, p::LinearSubspace, q::LinearSubspace) =
    set_subspaces!(
        H,
        convert(LinearSubspace{ComplexF64}, p),
        convert(LinearSubspace{ComplexF64}, q),
    )

function γ!(H::IntrinsicSubspaceHomotopy, t::Number)
    H.t_cache[] != t || return first(H.taylor_γ)
    if isreal(t)
        _γ!(H, real(t))
    else
        _γ!(H, t)
    end
    H.t_cache[] = t

    first(H.taylor_γ)
end
@inline function _γ!(H::IntrinsicSubspaceHomotopy, t::Number)
    A₁ = extrinsic(H.start).A
    A₀ = extrinsic(H.target).A
    b₁ = extrinsic(H.start).b
    b₀ = extrinsic(H.target).b

    H.A .= (1 .- t) .* A₁ .+ t .* A₀
    H.b .= (1 .- t) .* b₁ .+ t .* b₀

    H
end


function set_solution!(u::Vector, H::IntrinsicSubspaceHomotopy, x::AbstractVector, t)
    (length(x) == length(H.x) - 1) ||
        throw(ArgumentError("Cannot set solution. Expected extrinsic coordinates."))

    set_solution!(view(H.x, 1:length(x)), H.system, x)
    H.x[end] = 1

    if isone(t)
        LA.mul!(u, H.geodesic.γ1', H.x)
    elseif iszero(t)
        LA.mul!(u, H.geodesic.Q_cos', H.x)
    else
        LA.mul!(u, γ!(H, t)', H.x)
    end
end

function get_solution(H::IntrinsicSubspaceHomotopy, u::AbstractVector, t)
    if isone(t)
        (@view H.geodesic.γ1[1:end-1, :]) * u
    elseif iszero(t)
        (@view H.geodesic.Q_cos[1:end-1, :]) * u
    else
        γ = γ!(H, t)
        (@view γ[1:end-1, :]) * u
    end
end

function ModelKit.evaluate!(u, H::IntrinsicSubspaceHomotopy, v::Vector{ComplexDF64}, t)
    γ!(H, t)
    LA.mul!(H.x_high, H.A, v)
    H.x_high .+= H.b
    evaluate!(u, H.system, H.x_high)
    u
end

function ModelKit.evaluate!(u, H::IntrinsicSubspaceHomotopy, v::AbstractVector, t)
    γ!(H, t)
    LA.mul!(H.x, H.A, v)
    H.x .+= H.b
    evaluate!(u, H.system, H.x)
    u
end

function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::IntrinsicSubspaceHomotopy,
    v::AbstractVector,
    t,
)
    γ!(H, t)

    LA.mul!(H.x, H.A, v)
    H.x .+= H.b
    evaluate_and_jacobian!(u, H.J, H.system, H.x)
    LA.mul!(U, H.J, γ)

    nothing
end


# function γ̇!(H::IntrinsicSubspaceHomotopy, t_)
#     (t, dt) = t_

#     H.A .+
#     γ̇
# end


function ModelKit.taylor!(
    u,
    ::Val{1},
    H::IntrinsicSubspaceHomotopy,
    v,
    t,
    incr::Bool = true,
)
    γ = γ!(H, t)
    γ̇ = γ̇!(H, t)

    # apply chain rule
    #    d/dt [F(γ(t)v); (γ(t)v)[end] - 1] = [J_F(γ(t)v)* γ̇(t)*v;  (γ̇(t)v)[end]]

    H.v .= first.(v)
    LA.mul!(H.x, γ, H.v)
    LA.mul!(H.ẋ, γ̇, H.v)

    evaluate_and_jacobian!(u, H.J, H.system, H.x)
    LA.mul!(u, H.J, H.ẋ)
    M = size(H, 1)
    u[M] = H.ẋ[end]

    u
end

function _taylor_γ!(H::IntrinsicSubspaceHomotopy, t::Number)
    @unpack geodesic, taylor_γ = H
    @unpack Q, Q_cos, Θ, U = geodesic

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

function taylor_γ!(H::IntrinsicSubspaceHomotopy, t::Number)
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
    H::IntrinsicSubspaceHomotopy,
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

    n = first(size(H.system))
    u[n+1] = x²[end]

    u
end

function ModelKit.taylor!(
    u,
    ::Val{3},
    H::IntrinsicSubspaceHomotopy,
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

    n = first(size(H.system))
    u[n+1] = x³[end]

    u
end

function ModelKit.taylor!(
    u,
    ::Val{4},
    H::IntrinsicSubspaceHomotopy,
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

    n = first(size(H.system))
    u[n+1] = x⁴[end]

    u
end
