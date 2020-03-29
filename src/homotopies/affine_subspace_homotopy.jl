export AffineSubspaceHomotopy

struct AffineSubspaceHomotopy{S,T} <: AbstractHomotopy
    system::ModelKitSystem{S,T}

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
    taylor_x::NTuple{5,Vector{ComplexF64}}
    x_high::Vector{ComplexDF64}
end

AffineSubspaceHomotopy(F::ModelKit.System, start::AffineSubspace, target::AffineSubspace) =
    AffineSubspaceHomotopy(ModelKitSystem(F), start, target)

function AffineSubspaceHomotopy(
    system::ModelKitSystem,
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
    taylor_x = tuple((zeros(ComplexF64, size(Q, 1)) for i = 0:4)...)
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
        taylor_x,
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

function taylor_γ!(H::AffineSubspaceHomotopy, t::Number)
    H.t_cache[] != t || return H.taylor_γ

    @unpack Q, Q_cos, Θ, U, taylor_γ = H

    γ, γ¹, γ², γ³, γ⁴ = taylor_γ
    n, k = size(γ)
    @inbounds for j = 1:k
        Θⱼ = Θ[j]
        s, c = sincos(t * Θⱼ)
        s¹ = c * Θⱼ
        s² = -s * Θⱼ^2 / 2
        s³ = -c * Θⱼ^3 / 6
        s⁴ = s * (Θⱼ^2)^2 / 24
        c¹ = -s * Θⱼ
        c² = -c * Θⱼ^2 / 2
        c³ = s * Θⱼ^3 / 6
        c⁴ = c * (Θⱼ^2)^2 / 24
        for i = 1:n
            γ[i, j] = Q_cos[i, j] * c + Q[i, j] * s
            γ¹[i, j] = Q_cos[i, j] * c¹ + Q[i, j] * s¹
            γ²[i, j] = Q_cos[i, j] * c² + Q[i, j] * s²
            γ³[i, j] = Q_cos[i, j] * c³ + Q[i, j] * s³
            γ⁴[i, j] = Q_cos[i, j] * c⁴ + Q[i, j] * s⁴
        end
    end

    H.t_cache[] = t

    H.taylor_γ
end

function evaluate!(u, H::AffineSubspaceHomotopy, v::AbstractVector, t)
    γ, _ = taylor_γ!(H, t)
    x, _ = H.taylor_x
    n = first(size(H.system))
    if eltype(v) isa ComplexDF64
        LA.mul!(H.x_high, γ, v)
        evaluate!(u, H.system, H.x_high)
        u[n+1] = H.x_high[end] - 1
    else
        LA.mul!(x, γ, v)
        evaluate!(u, H.system, x)
        u[n+1] = x[end] - 1
    end
    u
end

function evaluate_and_jacobian!(u, U, H::AffineSubspaceHomotopy, v::AbstractVector, t)
    γ, _ = taylor_γ!(H::AffineSubspaceHomotopy, t)
    x, _ = H.taylor_x

    LA.mul!(x, γ, v)
    evaluate_and_jacobian!(u, H.J, H.system, x)
    LA.mul!(U, H.J, γ)

    n = first(size(H.system))
    u[n+1] = x[end] - 1

    m = length(v)
    for j = 1:m
        U[n+1, j] = γ[end, j]
    end

    nothing
end

function diff_t!(u, H::AffineSubspaceHomotopy, v, t, dv::Tuple, incremental::Bool)
    γ, γ¹, γ², γ³, γ⁴ = taylor_γ!(H, t)
    x, x¹, x², x³, x⁴ = H.taylor_x
    n = first(size(H.system))

    if length(dv) == 0
        LA.mul!(x, γ, v)
        LA.mul!(x¹, γ¹, v)
        diff_t!(u, Val(1), H.system, x, (x¹,))
        u[n+1] = x¹[end]
    elseif length(dv) == 1
        v¹, = dv
        if !incremental
            LA.mul!(x, γ, v)
            LA.mul!(x¹, γ¹, v)
        end

        LA.mul!(x¹, γ, v¹, true, true)

        LA.mul!(x², γ², v)
        LA.mul!(x², γ¹, v¹, true, true)

        diff_t!(u, Val(2), H.system, x, (x¹, x²))
        u[n+1] = x²[end]
    elseif length(dv) == 2
        v¹, v² = dv
        if !incremental
            LA.mul!(x, γ, v)
            LA.mul!(x¹, γ¹, v)
            LA.mul!(x¹, γ, v¹, true, true)
            LA.mul!(x², γ², v)
            LA.mul!(x², γ¹, v¹, true, true)
        end

        LA.mul!(x², γ, v², true, true)
        LA.mul!(x³, γ³, v)
        LA.mul!(x³, γ², v¹, true, true)
        LA.mul!(x³, γ¹, v², true, true)

        diff_t!(u, Val(3), H.system, x, (x¹, x², x³))
        u[n+1] = x³[end]
    elseif length(dv) == 3
        v¹, v², v³ = dv
        if !incremental
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
        LA.mul!(x⁴, γ⁴, v)
        LA.mul!(x⁴, γ³, v¹, true, true)
        LA.mul!(x⁴, γ², v², true, true)
        LA.mul!(x⁴, γ¹, v³, true, true)

        diff_t!(u, Val(4), H.system, x, (x¹, x², x³, x⁴))
        u[n+1] = x⁴[end]
    end
    u
end
