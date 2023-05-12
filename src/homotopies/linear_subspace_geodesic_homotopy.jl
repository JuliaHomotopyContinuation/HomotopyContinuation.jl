export LinearSubspaceGeodesicHomotopy, set_subspaces!

struct LinearSubspaceGeodesicInfo
    Q::Matrix{ComplexF64}
    Q_cos::Matrix{ComplexF64}
    Θ::Vector{Float64}
    U::Matrix{ComplexF64}
    γ1::Matrix{ComplexF64}
end

function LinearSubspaceGeodesicInfo(start::LinearSubspace, target::LinearSubspace)
    Q, Θ, U = geodesic_svd(target, start)
    Q_cos = target.intrinsic.Y * U
    γ1 = similar(Q_cos)
    n, k = size(γ1)
    for j = 1:k
        s, c = sincos(Θ[j])
        for i = 1:n
            γ1[i, j] = Q_cos[i, j] * c + Q[i, j] * s
        end
    end
    LinearSubspaceGeodesicInfo(Q, Q_cos, Θ, U, γ1)
end

const GEODESIC_LRU = LRU{
    Tuple{LinearSubspace{ComplexF64},LinearSubspace{ComplexF64}},
    LinearSubspaceGeodesicInfo,
}(
    maxsize = 128,
)

"""
    LinearSubspaceGeodesicHomotopy(F::System, V::LinearSubspace, W::LinearSubspace)
    LinearSubspaceGeodesicHomotopy(F::AbstractSystem, V::LinearSubspace, W::LinearSubspace)

Creates a homotopy ``H(x,t) = F(γ(t)x)`` where ``γ(t)`` is a family of affine subspaces
such that ``γ(1) = V`` and ``γ(0) = W``.
Here ``γ(t)`` is the geodesic between `V` and `W` in the affine Grassmanian, i.e.,
it is the curve of minimal length connecting `V` and `W`.
See also [`LinearSubspace`](@ref) and [`geodesic`](@ref) and the references therein.
"""
Base.@kwdef mutable struct LinearSubspaceGeodesicHomotopy{S<:AbstractSystem} <:
                           AbstractHomotopy
    system::S

    start::LinearSubspace{ComplexF64}
    target::LinearSubspace{ComplexF64}
    geodesic::LinearSubspaceGeodesicInfo

    J::Matrix{ComplexF64}
    x::Vector{ComplexF64}
    ẋ::Vector{ComplexF64}
    x_high::Vector{ComplexDF64}
    t_cache::Base.RefValue{ComplexF64}
    # For AD
    taylor_t_cache::Base.RefValue{ComplexF64}
    taylor_tṫ_cache::Base.RefValue{Tuple{ComplexF64,ComplexF64}}
    taylor_γ::NTuple{5,Matrix{ComplexF64}}
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

LinearSubspaceGeodesicHomotopy(
    F::ModelKit.System,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
) = LinearSubspaceGeodesicHomotopy(fixed(F; compile = compile), start, target)


function LinearSubspaceGeodesicHomotopy(
    system::AbstractSystem,
    start::LinearSubspace,
    target::LinearSubspace;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
)
    LinearSubspaceGeodesicHomotopy(
        system,
        convert(LinearSubspace{ComplexF64}, start),
        convert(LinearSubspace{ComplexF64}, target),
    )
end

function LinearSubspaceGeodesicHomotopy(
    system::AbstractSystem,
    start::LinearSubspace{ComplexF64},
    target::LinearSubspace{ComplexF64},
)
    geodesic = get!(GEODESIC_LRU, (start, target)) do
        LinearSubspaceGeodesicInfo(start, target)
    end
    Q = geodesic.Q
    tx⁴ = TaylorVector{5}(ComplexF64, size(Q, 1))
    LinearSubspaceGeodesicHomotopy(
        system = system,
        start = start,
        target = target,
        geodesic = LinearSubspaceGeodesicInfo(start, target),
        J = zeros(ComplexF64, size(system) .+ (1, 1)),
        x = zeros(ComplexF64, size(Q, 1)),
        ẋ = zeros(ComplexF64, size(Q, 1)),
        x_high = zeros(ComplexDF64, size(Q, 1)),
        t_cache = Ref(complex(NaN, NaN)),
        taylor_t_cache = Ref(complex(NaN, NaN)),
        taylor_tṫ_cache = Ref((complex(NaN, NaN), complex(NaN, NaN))),
        taylor_γ = tuple((similar(Q) for i = 0:4)...),
        v = zeros(ComplexF64, size(Q, 2)),
        tx⁴ = tx⁴,
        tx³ = TaylorVector{4}(tx⁴),
        tx² = TaylorVector{3}(tx⁴),
        tx¹ = TaylorVector{2}(tx⁴),
    )
end
Base.size(H::LinearSubspaceGeodesicHomotopy) = (size(H.system)[1] + 1, dim(H.start) + 1)

"""
    set_subspaces!(H::LinearSubspaceGeodesicHomotopy, start::LinearSubspace, target::LinearSubspace)

Update the homotopy `H` to track from the affine subspace `start` to `target`.
"""
function set_subspaces!(
    H::LinearSubspaceGeodesicHomotopy,
    start::LinearSubspace,
    target::LinearSubspace,
)
    H.start = start
    H.target = target
    H.geodesic = get!(GEODESIC_LRU, (start, target)) do
        LinearSubspaceGeodesicInfo(start, target)
    end
    H
end
start_parameters!(H::LinearSubspaceGeodesicHomotopy, p::LinearSubspace) =
    set_subspaces!(H, convert(LinearSubspace{ComplexF64}, p), H.target)
target_parameters!(H::LinearSubspaceGeodesicHomotopy, q::LinearSubspace) =
    set_subspaces!(H, H.start, convert(LinearSubspace{ComplexF64}, q))
parameters!(H::LinearSubspaceGeodesicHomotopy, p::LinearSubspace, q::LinearSubspace) =
    set_subspaces!(
        H,
        convert(LinearSubspace{ComplexF64}, p),
        convert(LinearSubspace{ComplexF64}, q),
    )

function γ!(H::LinearSubspaceGeodesicHomotopy, t::Number)
    # H.t_cache[] != t || return first(H.taylor_γ)
    if isreal(t)
        _γ!(H, real(t))
    else
        _γ!(H, t)
    end
    H.t_cache[] = t

    first(H.taylor_γ)
end
@inline function _γ!(H::LinearSubspaceGeodesicHomotopy, t::Number)
    @unpack Q, Q_cos, Θ = H.geodesic
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

γ̇!(H::LinearSubspaceGeodesicHomotopy, t::Number) =
    isreal(t) ? _γ̇!(H, real(t)) : _γ̇!(H, t)
@inline function _γ̇!(H::LinearSubspaceGeodesicHomotopy, t::Number)
    @unpack Q, Q_cos, Θ = H.geodesic
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

function set_solution!(u::Vector, H::LinearSubspaceGeodesicHomotopy, x::AbstractVector, t)
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

function get_solution(H::LinearSubspaceGeodesicHomotopy, u::AbstractVector, t)
    if isone(t)
        (@view H.geodesic.γ1[1:end-1, :]) * u
    elseif iszero(t)
        (@view H.geodesic.Q_cos[1:end-1, :]) * u
    else
        γ = γ!(H, t)
        (@view γ[1:end-1, :]) * u
    end
end

function ModelKit.evaluate!(u, H::LinearSubspaceGeodesicHomotopy, v::Vector{ComplexDF64}, t)
    γ = γ!(H, t)
    n = first(size(H.system))
    LA.mul!(H.x_high, γ, v)
    evaluate!(u, H.system, H.x_high)
    u[n+1] = H.x_high[end] - 1.0
end

function ModelKit.evaluate!(u, H::LinearSubspaceGeodesicHomotopy, v::AbstractVector, t)
    γ = γ!(H, t)
    n = first(size(H.system))
    LA.mul!(H.x, γ, v)
    evaluate!(u, H.system, H.x)
    u[n+1] = H.x[end] - 1.0
    u
end

function ModelKit.evaluate_and_jacobian!(
    u,
    U,
    H::LinearSubspaceGeodesicHomotopy,
    v::AbstractVector,
    t,
)
    γ = γ!(H, t)

    LA.mul!(H.x, γ, v)
    evaluate_and_jacobian!(u, H.J, H.system, H.x)
    LA.mul!(U, H.J, γ)

    n = first(size(H.system))
    u[n+1] = H.x[end] - 1

    m = length(v)
    for j = 1:m
        U[n+1, j] = γ[end, j]
    end

    nothing
end

# function ModelKit.taylor!(
#     u,
#     ::Val{1},
#     H::LinearSubspaceGeodesicHomotopy,
#     v,
#     tṫ,
#     incr::Bool = true,
# )
#     t, ṫ = tṫ
#     γ = γ!(H, t)
#     γ̇ = γ̇!(H, t)
#     γ̇ .*= ṫ


#     # apply chain rule
#     #    d/dt [F(γ(t)v); (γ(t)v)[end] - 1] = [J_F(γ(t)v)* γ̇(t)*v;  (γ̇(t)v)[end]]

#     H.v .= first.(v)
#     LA.mul!(H.x, γ, H.v)
#     LA.mul!(H.ẋ, γ̇, H.v)

#     evaluate_and_jacobian!(u, H.J, H.system, H.x)
#     LA.mul!(u, H.J, H.ẋ)
#     M = size(H, 1)
#     if eltype(u) <: TruncatedTaylorSeries
#         u[M, 1] = H.x[end] - 1
#         u[M, 2] = H.ẋ[end]
#     else
#         u[M] = H.ẋ[end]
#     end

#     u
# end

function _taylor_γ!(H::LinearSubspaceGeodesicHomotopy, tṫ)
    @unpack geodesic, taylor_γ = H
    @unpack Q, Q_cos, Θ, U = geodesic

    t, ṫ = tṫ
    γ, γ¹, γ², γ³, γ⁴ = taylor_γ
    n, k = size(γ)
    @inbounds for j = 1:k
        tΘⱼ = t * Θ[j]
        dt_tΘⱼ = ṫ * Θ[j]
        s, c = sincos(tΘⱼ)
        c¹ = -s * dt_tΘⱼ
        s¹ = c * dt_tΘⱼ
        dt_tΘⱼ_2 = 0.5 * dt_tΘⱼ^2
        c² = -c * dt_tΘⱼ_2
        s² = -s * dt_tΘⱼ_2
        dt_tΘⱼ_3 = dt_tΘⱼ_2 * dt_tΘⱼ / 3
        c³ = s * dt_tΘⱼ_3
        s³ = -c * dt_tΘⱼ_3
        dt_tΘⱼ_4 = 0.25 * dt_tΘⱼ_3 * dt_tΘⱼ
        c⁴ = c * dt_tΘⱼ_4
        s⁴ = s * dt_tΘⱼ_4
        for i = 1:n
            γ[i, j] = Q_cos[i, j] * c + Q[i, j] * s
            γ¹[i, j] = Q_cos[i, j] * c¹ + Q[i, j] * s¹
            γ²[i, j] = Q_cos[i, j] * c² + Q[i, j] * s²
            γ³[i, j] = Q_cos[i, j] * c³ + Q[i, j] * s³
            γ⁴[i, j] = Q_cos[i, j] * c⁴ + Q[i, j] * s⁴
        end
    end

end

function taylor_γ!(H::LinearSubspaceGeodesicHomotopy, tṫ::Tuple{<:Number,<:Number})
    # H.taylor_tṫ_cache[] != tṫ || return H.taylor_γ

    _taylor_γ!(H, tṫ)
    H.taylor_tṫ_cache[] = tṫ

    H.taylor_γ
end


function ModelKit.taylor!(
    u,
    ::Val{1},
    H::LinearSubspaceGeodesicHomotopy,
    tv::TaylorVector,
    t,
    incr::Bool = false,
)
    γ, γ¹ = taylor_γ!(H, t)
    x, x¹ = vectors(H.tx¹)
    v, = vectors(tv)


    LA.mul!(x, γ, v)
    LA.mul!(x¹, γ¹, v)
    taylor!(u, Val(1), H.system, H.tx¹)

    n = first(size(H.system))

    if eltype(u) <: TruncatedTaylorSeries
        u[n+1] = H.tx¹[end]
        u[n+1, 1] -= 1
    else
        u[n+1] = x¹[end]
    end

    u
end

function ModelKit.taylor!(
    u,
    ::Val{2},
    H::LinearSubspaceGeodesicHomotopy,
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
    LA.mul!(x¹, γ, v¹, true, true)

    LA.mul!(x², γ¹, v¹)
    LA.mul!(x², γ², v, true, true)

    taylor!(u, Val(2), H.system, H.tx²)

    n = first(size(H.system))

    if eltype(u) <: TruncatedTaylorSeries
        u[n+1] = H.tx²[end]
        u[n+1, 1] -= 1
    else
        u[n+1] = x²[end]
    end

    u
end

function ModelKit.taylor!(
    u,
    ::Val{3},
    H::LinearSubspaceGeodesicHomotopy,
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
    if eltype(u) <: TruncatedTaylorSeries
        u[n+1] = H.tx³[end]
        u[n+1, 1] -= 1
    else
        u[n+1] = x³[end]
    end

    u
end

function ModelKit.taylor!(
    u,
    ::Val{4},
    H::LinearSubspaceGeodesicHomotopy,
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
    if eltype(u) <: TruncatedTaylorSeries
        u[n+1] = H.tx⁴[end]
        u[n+1, 1] -= 1
    else
        u[n+1] = x⁴[end]
    end

    u
end
