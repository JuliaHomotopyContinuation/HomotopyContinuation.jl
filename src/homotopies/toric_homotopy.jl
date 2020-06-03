###
### ToricHomotopy
###
# This is part of the implementation of polyhedral homotopy, namely to reverse the toric
# degeneration from 0 to 1.
# This is only an internal homotopy.

const StructVectorComplexF64 = StructArrays.StructArray{
    Complex{Float64},
    1,
    NamedTuple{(:re, :im),Tuple{Vector{Float64},Vector{Float64}}},
    Int64,
}

Base.@kwdef struct ToricHomotopy{S} <: AbstractHomotopy
    system::ModelKit.CompiledSystem{S}
    system_coeffs::Vector{ComplexF64}
    weights::Vector{Float64}

    #t^w
    t_weights::Vector{Float64}
    #t^w for t complex or negative
    complex_t_weights::StructVectorComplexF64
    coeffs::Vector{ComplexF64}
    dt_coeffs::Vector{ComplexF64}
    x::Vector{ComplexF64}
    # cache value
    t_coeffs::Base.RefValue{ComplexF64}
    t_taylor_coeffs::Base.RefValue{ComplexF64}

    # taylor coefficients
    taylor_coeffs::TaylorVector{5,ComplexF64}
    tc3::TaylorVector{4,ComplexF64}
    tc2::TaylorVector{3,ComplexF64}
end

function ToricHomotopy(
    system::ModelKit.CompiledSystem,
    system_coeffs::Vector{Vector{ComplexF64}},
)
    m = ModelKit.nparameters(system)
    m1 = sum(length, system_coeffs)
    m == m1 ||
        throw(ArgumentError("System parameters and coefficients do not have the same size, got $m and $m1"))
    n = nvariables(system)
    taylor_coeffs = TaylorVector{5}(ComplexF64, m)
    ToricHomotopy(
        system = system,
        system_coeffs = reduce(vcat, system_coeffs),
        weights = zeros(m),
        t_weights = zeros(m),
        complex_t_weights = StructArrays.StructArray(zeros(ComplexF64, m)),
        coeffs = zeros(ComplexF64, m),
        dt_coeffs = zeros(ComplexF64, m),
        x = zeros(ComplexF64, n),
        t_coeffs = Ref(complex(0.0, 0.0)),
        t_taylor_coeffs = Ref(complex(0.0, 0.0)),
        taylor_coeffs = taylor_coeffs,
        tc3 = TaylorVector{4}(taylor_coeffs),
        tc2 = TaylorVector{3}(taylor_coeffs),
    )
end

Base.size(H::ToricHomotopy) = size(H.system)

function update_weights!(
    H::ToricHomotopy,
    support::AbstractVector{<:AbstractMatrix},
    lifting::AbstractVector{<:AbstractVector},
    cell::MixedSubdivisions.MixedCell;
    min_weight::Union{Nothing,Float64} = nothing,
    max_weight::Union{Nothing,Float64} = nothing,
)
    l = 1
    s_max, s_min = 0.0, Inf
    n = length(cell.normal)
    for (i, Aᵢ) in enumerate(support)
        wᵢ = lifting[i]
        βᵢ = cell.β[i]
        mᵢ = size(Aᵢ, 2)
        aᵢ, bᵢ = cell.indices[i]
        for j = 1:mᵢ
            if j == aᵢ || j == bᵢ
                H.weights[l] = 0.0
            else
                sᵢⱼ = wᵢ[j] - βᵢ
                for k = 1:n
                    @inbounds sᵢⱼ += Aᵢ[k, j] * cell.normal[k]
                end
                H.weights[l] = sᵢⱼ
                s_max = max(s_max, sᵢⱼ)
                s_min = min(s_min, sᵢⱼ)
            end
            l += 1
        end
    end

    if min_weight !== nothing
        λ = s_min / min_weight
        H.weights ./= λ
        s_min, s_max = min_weight, s_max / λ
    elseif max_weight !== nothing
        λ = s_max / max_weight
        H.weights ./= λ
        s_min, s_max = s_min / λ, max_weight
    end

    H.t_coeffs[] = NaN
    H.t_taylor_coeffs[] = NaN

    s_min, s_max
end

function evaluate_weights!(u::Vector{Float64}, weights::Vector{Float64}, t::Float64)
    s = log(t)
    n = length(weights)
    LoopVectorization.@avx for i = 1:n
        u[i] = exp(weights[i] * s)
    end
    u
end

function evaluate_weights!(
    u::Vector{Float64},
    v::Vector{Float64},
    weights::Vector{Float64},
    t::ComplexF64,
)
    s = log(t)
    n = length(weights)
    sr, si = reim(s)
    LoopVectorization.@avx for i = 1:n
        r = exp(weights[i] * sr)
        wsi = weights[i] * si
        # TODO: This should actually by s, c = sincos(wsi)
        # However, this it not yet supported by LoopVectorization
        s = sin(wsi)
        c = cos(wsi)
        u[i] = r * c
        v[i] = r * s
    end
end


function coeffs!(H::ToricHomotopy, t)
    # H.t_coeffs[] == t && return H.coeffs
    c = H.coeffs
    if !isreal(t) || real(t) < 0
        evaluate_weights!(
            H.complex_t_weights.re,
            H.complex_t_weights.im,
            H.weights,
            complex(t),
        )
        @inbounds for i = 1:length(c)
            uᵢ = H.system_coeffs[i]
            ctwᵢ = H.complex_t_weights[i]
            c[i] = uᵢ * ctwᵢ
        end
    elseif iszero(t)
        @inbounds for i = 1:length(c)
            wᵢ = H.weights[i]
            uᵢ = H.system_coeffs[i]
            c[i] = ifelse(wᵢ == 0, uᵢ, complex(0.0, 0.0))
        end
    else
        evaluate_weights!(H.t_weights, H.weights, real(t))
        @inbounds for i = 1:length(c)
            uᵢ = H.system_coeffs[i]
            twᵢ = H.t_weights[i]
            c[i] = uᵢ * twᵢ
        end
    end

    H.t_coeffs[] = t
    c
end


function dt_coeffs!(H::ToricHomotopy, t)
    if iszero(t)
        @inbounds for i = 1:length(H.dt_coeffs)
            wᵢ = H.weights[i]
            uᵢ = H.system_coeffs[i]
            H.dt_coeffs[i] = ifelse(wᵢ == 1, uᵢ, complex(0.0, 0.0))
        end
    else
        c = coeffs!(H, t)
        if isreal(t)
            s = real(t)
            @fastmath s⁻¹ = inv(s)
            for i = 1:length(c)
                @inbounds H.dt_coeffs[i] = H.weights[i] * c[i] * s⁻¹
            end
        else
            @fastmath t⁻¹ = inv(t)
            for i = 1:length(c)
                @inbounds H.dt_coeffs[i] = H.weights[i] * c[i] * t⁻¹
            end
        end
    end

    H.dt_coeffs
end

function evaluate!(u, H::ToricHomotopy, x::AbstractVector, t)
    ModelKit.evaluate!(u, H.system, x, coeffs!(H, t))
end

function evaluate_and_jacobian!(u, U, H::ToricHomotopy, x::AbstractVector, t)
    ModelKit.evaluate_and_jacobian!(u, U, H.system, x, coeffs!(H, t))
    nothing
end


function taylor!(u, v::Val{1}, H::ToricHomotopy, x::Vector, t)
    ModelKit.evaluate!(u, H.system, x, dt_coeffs!(H, t))
end

function taylor_coeffs!(H::ToricHomotopy, t::Real)
    H.t_taylor_coeffs[] != t || return H.taylor_coeffs

    tc = H.taylor_coeffs
    # c, c¹, c², c³, c⁴ = H.taylor_coeffs
    if iszero(t)
        for i = 1:length(tc)
            wᵢ = H.weights[i]
            uᵢ = H.system_coeffs[i]
            if wᵢ < 1e-12
                tc[i, 1] = uᵢ
                tc[i, 2] = tc[i, 3] = tc[i, 4] = tc[i, 5] = 0.0
            elseif isapprox(wᵢ, 1)
                tc[i, 2] = uᵢ
                tc[i, 1] = tc[i, 3] = tc[i, 4] = tc[i, 5] = 0.0
            else
                tc[i, 1] = tc[i, 2] = tc[i, 3] = tc[i, 4] = tc[i, 5] = 0.0
            end
        end
    else
        if t !== H.t_coeffs[]
            evaluate_weights!(H.t_weights, H.weights, t)
        end
        t⁻¹ = inv(t)
        @inbounds for i = 1:length(tc)
            wᵢ = H.weights[i]
            uᵢ = H.system_coeffs[i]
            twᵢ = H.t_weights[i]
            tc[i, 1] = uᵢ * twᵢ
            tw1 = wᵢ * twᵢ * t⁻¹
            tc[i, 2] = uᵢ * tw1
            tw2 = 0.5 * (wᵢ - 1) * tw1 * t⁻¹
            tc[i, 3] = uᵢ * tw2
            tw3 = (wᵢ - 2) * tw2 * t⁻¹ / 3
            tc[i, 4] = uᵢ * tw3
            tw4 = 0.25 * (wᵢ - 3) * tw3 * t⁻¹
            tc[i, 5] = uᵢ * tw4
        end

    end

    H.t_taylor_coeffs[] = t

    H.taylor_coeffs
end


function taylor!(u, v::Val{2}, H::ToricHomotopy, tx::TaylorVector, t)
    taylor_coeffs!(H, real(t))
    ModelKit.taylor!(u, v, H.system, tx, H.tc2)
end
function taylor!(u, v::Val{3}, H::ToricHomotopy, tx::TaylorVector, t)
    taylor_coeffs!(H, real(t))
    ModelKit.taylor!(u, v, H.system, tx, H.tc3)
end
function taylor!(u, v::Val{4}, H::ToricHomotopy, tx::TaylorVector, t)
    taylor_coeffs!(H, real(t))
    ModelKit.taylor!(u, v, H.system, tx, H.taylor_coeffs)
end
