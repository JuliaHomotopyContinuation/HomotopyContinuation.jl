export PolyhedralHomotopy

const StructVectorComplexF64 = StructArrays.StructArray{
    Complex{Float64},
    1,
    NamedTuple{(:re, :im),Tuple{Vector{Float64},Vector{Float64}}},
    Int64,
}

"""
    PolyhedralHomotopy(G, F; gamma=exp(i * 2π*rand()))
Construct the homotopy ``H(x, t) = γtG(x) + (1-t)F(x)``.
"""
struct PolyhedralHomotopy{S} <: AbstractHomotopy
    system::ModelKit.CompiledSystem{S}
    system_coeffs::Vector{ComplexF64}
    weights::Vector{Float64}

    t_weights::Vector{Float64}
    complex_t_weights::StructVectorComplexF64
    t_coeffs::Base.RefValue{Float64}
    taylor_coeffs::NTuple{5,Vector{ComplexF64}}
    # these are just here to avoid unnecessary allocations
    dc1::Tuple{Vector{ComplexF64}}
    dc2::NTuple{2,Vector{ComplexF64}}
    dc3::NTuple{3,Vector{ComplexF64}}
    dc4::NTuple{4,Vector{ComplexF64}}
end

function PolyhedralHomotopy(
    system::ModelKit.CompiledSystem,
    system_coeffs::Vector{Vector{ComplexF64}},
)
    m = ModelKit.nparameters(system)
    m1 = sum(length, system_coeffs)
    m == m1 ||
    throw(ArgumentError("System parameters and coefficients do not have the same size, got $m and $m1"))

    n = nvariables(system)

    weights = zeros(m)

    t_weights = zeros(m)
    complex_t_weights = StructArrays.StructArray(zeros(ComplexF64, m))
    t_coeffs = Ref(0.0)
    taylor_coeffs = tuple((zeros(ComplexF64, m) for i = 0:4)...)
    dc1 = (taylor_coeffs[2],)
    dc2 = (taylor_coeffs[2], taylor_coeffs[3])
    dc3 = (taylor_coeffs[2], taylor_coeffs[3], taylor_coeffs[4])
    dc4 = (taylor_coeffs[2], taylor_coeffs[3], taylor_coeffs[4], taylor_coeffs[5])


    PolyhedralHomotopy(
        system,
        reduce(vcat, system_coeffs),
        weights,
        t_weights,
        complex_t_weights,
        t_coeffs,
        taylor_coeffs,
        dc1,
        dc2,
        dc3,
        dc4,
    )
end

Base.size(H::PolyhedralHomotopy) = size(H.system)

function update_weights!(
    H::PolyhedralHomotopy,
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


function taylor_coeffs!(H::PolyhedralHomotopy, t::Real)
    H.t_coeffs[] != t || return H.taylor_coeffs

    c, c¹, c², c³, c⁴ = H.taylor_coeffs

    if t == 0
        for i = 1:length(c)
            wᵢ = H.weights[i]
            uᵢ = H.system_coeffs[i]
            twᵢ = H.t_weights[i]
            if wᵢ == 0
                c[i] = uᵢ
                c¹[i] = c²[i] = c³[i] = c⁴[i] = 0.0
            elseif wᵢ == 1
                c¹[i] = uᵢ
                c[i] = c²[i] = c³[i] = c⁴[i] = 0.0
            elseif wᵢ == 2
                c²[i] = 0.5 * uᵢ
                c[i] = c¹[i] = c³[i] = c⁴[i] = 0.0
            elseif wᵢ == 3
                c³[i] = uᵢ / 6
                c[i] = c¹[i] = c²[i] = c⁴[i] = 0.0
            elseif wᵢ == 4
                c⁴[i] = uᵢ / 24
                c[i] = c¹[i] = c²[i] = c³[i] = 0.0
            else
                c[i] = c¹[i] = c²[i] = c³[i] = c⁴[i] = 0.0
            end
        end
    else
        evaluate_weights!(H.t_weights, H.weights, t)
        t⁻¹ = inv(t)
        for i = 1:length(c)
            wᵢ = H.weights[i]
            uᵢ = H.system_coeffs[i]
            twᵢ = H.t_weights[i]
            c[i] = uᵢ * twᵢ
            tw1 = wᵢ * twᵢ * t⁻¹
            c¹[i] = uᵢ * tw1
            tw2 = 0.5 * (wᵢ - 1) * tw1 * t⁻¹
            c²[i] = uᵢ * tw2
            tw3 = (wᵢ - 2) * tw2 * t⁻¹ / 3
            c³[i] = uᵢ * tw3
            tw4 = 0.25 * (wᵢ - 3) * tw3 * t⁻¹
            c⁴[i] = uᵢ * tw4
        end

    end

    H.t_coeffs[] = t

    H.taylor_coeffs
end

function coeffs!(H::PolyhedralHomotopy, t::Real)
    if t < 0
        c, _ = H.taylor_coeffs
        evaluate_weights!(
            H.complex_t_weights.re,
            H.complex_t_weights.im,
            H.weights,
            complex(t),
        )
        for i = 1:length(c)
            wᵢ = H.weights[i]
            uᵢ = H.system_coeffs[i]
            twᵢ = H.complex_t_weights[i]
            c[i] = uᵢ * twᵢ
        end
        return c
    else
        c, _ = taylor_coeffs!(H, t)
        return c
    end
end

function evaluate!(u, H::PolyhedralHomotopy, x::AbstractVector, t)
    c = coeffs!(H::PolyhedralHomotopy, real(t))
    ModelKit.evaluate!(u, H.system, x, c)
end

function evaluate_and_jacobian!(u, U, H::PolyhedralHomotopy, x::AbstractVector, t)
    c = coeffs!(H::PolyhedralHomotopy, real(t))
    ModelKit.evaluate_and_jacobian!(u, U, H.system, x, c)
    nothing
end

function diff_t!(u, H::PolyhedralHomotopy, x, t, dx::Tuple)
    c, _ = taylor_coeffs!(H::PolyhedralHomotopy, real(t))
    k = length(dx)
    if k == 0
        ModelKit.diff_t!(u, H.system, x, dx, c, H.dc1)
    elseif k == 1
        ModelKit.diff_t!(u, H.system, x, dx, c, H.dc2)
    elseif k == 2
        ModelKit.diff_t!(u, H.system, x, dx, c, H.dc3)
    elseif k == 3
        ModelKit.diff_t!(u, H.system, x, dx, c, H.dc4)
    end
    u
end
