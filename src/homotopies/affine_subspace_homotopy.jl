export AffineSubspaceHomotopy

struct AffineSubspaceHomotopy{S<:AbstractSystem} <: AbstractHomotopy
    system::S

    start::AffineSubspace{ComplexF64}
    target::AffineSubspace{ComplexF64}

    J::Matrix{ComplexF64}
    Q::Matrix{ComplexF64}
    Q_cos::Matrix{ComplexF64}
    Θ::Vector{Float64}
    U::Matrix{ComplexF64}
    γ1::Matrix{ComplexF64}
    x_stiefel::Vector{ComplexF64}

    t_cache::Base.RefValue{ComplexF64}
    taylor_γ::NTuple{5,Matrix{ComplexF64}}
    v::Vector{ComplexF64}
    tx⁴::TaylorVector{5,ComplexF64}
    tx³::TaylorVector{4,ComplexF64}
    tx²::TaylorVector{3,ComplexF64}
    tx¹::TaylorVector{2,ComplexF64}
    x::Vector{ComplexF64}
    x_high::Vector{ComplexDF64}
end

AffineSubspaceHomotopy(F::ModelKit.System, start::AffineSubspace, target::AffineSubspace) =
    AffineSubspaceHomotopy(ModelKitSystem(F), start, target)

function AffineSubspaceHomotopy(
    system::AbstractSystem,
    start::AffineSubspace,
    target::AffineSubspace,
)
    J = zeros(ComplexF64, size(system) .+ (1, 1))
    Q, Θ, U = geodesic_svd(target, start)
    Q_cos = target.intrinsic.Y * U
    γ1 = Q_cos * LA.diagm(0 => cos.(Θ)) + Q * LA.diagm(0 => sin.(Θ))
    x_stiefel = zeros(ComplexF64, size(Q, 1))
    t_cache = Ref(Complex(NaN, NaN))
    taylor_γ = tuple((similar(Q) for i = 0:4)...)
    v = zeros(ComplexF64, size(Q, 2))
    tx⁴ = TaylorVector{5}(ComplexF64, size(Q, 1))
    x = zeros(ComplexF64, size(Q, 1))
    x_high = zeros(ComplexDF64, size(Q, 1))
    AffineSubspaceHomotopy(
        system,
        start,
        target,
        J,
        Q,
        Q_cos,
        Θ,
        U,
        γ1,
        x_stiefel,
        t_cache,
        taylor_γ,
        v,
        tx⁴,
        TaylorVector{4}(tx⁴),
        TaylorVector{3}(tx⁴),
        TaylorVector{2}(tx⁴),
        x,
        x_high,
    )
end
Base.size(H::AffineSubspaceHomotopy) = (size(H.system)[1] + 1, dim(H.start) + 1)

function setup!(H::AffineSubspaceHomotopy, start::AffineSubspace, target::AffineSubspace)
    Q, Θ, U = geodesic_svd(target, start)
    copy!(H.start, start)
    copy!(H.target, target)
    H.Q .= Q
    H.Θ .= Θ
    H.U .= U
    LA.mul!(H.Q_cos, target.intrinsic.Y, U)
    LA.mul!(H.γ1, Q, LA.diagm(0 => sin.(Θ)))
    LA.mul!(H.γ1, H.Q_cos, LA.diagm(0 => cos.(Θ)), true, true)
    H
end

function set_solution!(u::Vector, H::AffineSubspaceHomotopy, x::AbstractVector, t)
    (length(x) == length(H.x_stiefel) - 1) ||
    throw(ArgumentError("Cannot set solution. Expected extrinsic coordinates."))
    for i = 1:length(x)
        H.x_stiefel[i] = x[i]
    end
    H.x_stiefel[end] = 1

    if t == 1
        LA.mul!(u, H.γ1', H.x_stiefel)
    elseif t == 0
        LA.mul!(u, H.Q_cos', H.x_stiefel)
    else
        γ, _ = taylor_γ!(H, t)
        LA.mul!(u, γ', H.x_stiefel)
    end
end

function get_solution(H::AffineSubspaceHomotopy, u::AbstractVector, t)
    if t == 1
        (@view H.γ1[1:end-1, :]) * u
    elseif t == 0
        (@view H.Q_cos[1:end-1, :]) * u
    else
        γ, _ = taylor_γ!(H, t)
        (@view γ[1:end-1, :]) * u
    end
end


function _taylor_γ!(H::AffineSubspaceHomotopy, t::Number)
    @unpack Q, Q_cos, Θ, U, taylor_γ = H

    γ, γ¹, γ², γ³, γ⁴ = taylor_γ
    n, k = size(γ)
    @inbounds for j = 1:k
        Θⱼ = Θ[j]
        s, c = sincos(t * Θⱼ)
        c¹ = -s * Θⱼ
        s¹ = c * Θⱼ
        Θⱼ_2 = 0.5 * Θⱼ^2
        s² = -s * Θⱼ_2
        c² = -c * Θⱼ_2
        Θⱼ_3 = Θⱼ_2 * Θⱼ / 3
        s³ = -c * Θⱼ_3
        c³ = s * Θⱼ_3
        Θⱼ_4 = 0.25 * Θⱼ_3 * Θⱼ
        s⁴ = s * Θⱼ_4
        c⁴ = c * Θⱼ_4
        for i = 1:n
            γ[i, j] = Q_cos[i, j] * c + Q[i, j] * s
            γ¹[i, j] = Q_cos[i, j] * c¹ + Q[i, j] * s¹
            γ²[i, j] = Q_cos[i, j] * c² + Q[i, j] * s²
            γ³[i, j] = Q_cos[i, j] * c³ + Q[i, j] * s³
            γ⁴[i, j] = Q_cos[i, j] * c⁴ + Q[i, j] * s⁴
        end
    end

end

function taylor_γ!(H::AffineSubspaceHomotopy, t::Number)
    H.t_cache[] != t || return H.taylor_γ

    if isreal(t)
        _taylor_γ!(H, real(t))
    else
        _taylor_γ!(H, t)
    end

    H.t_cache[] = t

    H.taylor_γ
end

function evaluate!(u, H::AffineSubspaceHomotopy, v::AbstractVector, t)
    γ, _ = taylor_γ!(H, t)
    n = first(size(H.system))
    if eltype(v) isa ComplexDF64
        LA.mul!(H.x_high, γ, v)
        evaluate!(u, H.system, H.x_high)
        u[n+1] = H.x_high[end] - 1
    else
        LA.mul!(H.x, γ, v)
        evaluate!(u, H.system, H.x)
        u[n+1] = H.x[end] - 1
    end
    u
end

function evaluate_and_jacobian!(u, U, H::AffineSubspaceHomotopy, v::AbstractVector, t)
    γ, _ = taylor_γ!(H::AffineSubspaceHomotopy, t)

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

function taylor!(u, ::Val{1}, H::AffineSubspaceHomotopy, tv::TaylorVector, t, incr::Bool)
    γ, γ¹ = taylor_γ!(H, t)
    x, x¹ = vectors(H.tx¹)
    v, = vectors(tv)

    H.v .= v
    LA.mul!(H.x, γ, H.v)
    x .= H.x
    LA.mul!(H.x, γ¹, H.v)
    x¹ .= H.x
    taylor!(u, Val(1), H.system, H.tx¹)

    n = first(size(H.system))
    u[n+1] = x¹[end]

    u
end

function taylor!(u, v::Val{2}, H::AffineSubspaceHomotopy, tv::TaylorVector, t, incr::Bool)
    γ, γ¹, γ², γ³, γ⁴ = taylor_γ!(H, t)
    x, x¹, x² = vectors(H.tx²)
    v, v¹ = vectors(tv)

    if !incr
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

function taylor!(u, v::Val{3}, H::AffineSubspaceHomotopy, tv::TaylorVector, t, incr::Bool)
    γ, γ¹, γ², γ³, γ⁴ = taylor_γ!(H, t)
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

function taylor!(u, v::Val{4}, H::AffineSubspaceHomotopy, tv::TaylorVector, t, incr::Bool)
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
    H.v .= v³
    LA.mul!(H.x, γ, H.v)
    x³ .+= H.x

    LA.mul!(H.x, γ¹, v³)

    H.v .= v
    LA.mul!(H.x, γ⁴, v, true, true)
    H.v .= v¹
    LA.mul!(H.x, γ³, H.v, true, true)
    H.v .= v²
    LA.mul!(H.x, γ², H.v, true, true)
    x⁴ .= H.x

    taylor!(u, Val(4), H.system, H.tx⁴)

    n = first(size(H.system))
    u[n+1] = x⁴[end]

    u
end
