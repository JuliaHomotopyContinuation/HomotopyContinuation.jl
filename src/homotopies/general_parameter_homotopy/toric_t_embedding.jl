const StructVectorComplexF64 = StructArrays.StructVector{
    ComplexF64,
    NamedTuple{(:re, :im),Tuple{Vector{Float64},Vector{Float64}}},
    Int64,
}


struct ToricTEmbedding <: TEmbedding
    # c[i]*t^w[i]
    coeffs::Vector{ComplexF64}
    support::Vector{Matrix{Int}}
    weights::Vector{Float64}
    complex_t_weights::StructVectorComplexF64
    t_cache_weights::Ref{ComplexF64}
end

function ToricTEmbedding(
    coeffs::Vector{ComplexF64},
    support::AbstractVector{<:AbstractMatrix},
)
    ToricTEmbedding(
        coeffs,
        map(s -> convert(Matrix{Int}, s), support),
        zeros(length(coeffs)),
        StructArrays.StructArray(zeros(ComplexF64, length(coeffs))),
        Ref(ComplexF64(NaN)),
    )
end


function update!(
    H::ToricTEmbedding,
    lifting::AbstractVector{<:AbstractVector},
    cell::MixedSubdivisions.MixedCell;
    min_weight::Union{Nothing,Float64} = nothing,
    max_weight::Union{Nothing,Float64} = nothing,
)
    support = H.support
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

    H.t_cache_weights[] = ComplexF64(NaN)

    s_min, s_max
end


function evaluate_weights!(T::ToricTEmbedding, t::Number)
    t == T.t_cache_weights[] && return T.complex_t_weights

    u, v = T.complex_t_weights.re, T.complex_t_weights.im

    weights = T.weights
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

    T.t_cache_weights[] = t

    T.complex_t_weights
end


function ModelKit.evaluate!(c, H::ToricTEmbedding, t)
    if iszero(t)
        @inbounds for i = 1:length(c)
            wᵢ = H.weights[i]
            uᵢ = H.coeffs[i]
            c[i] = ifelse(wᵢ == 0, uᵢ, complex(0.0, 0.0))
        end
    else
        evaluate_weights!(H, t)
        @inbounds for i = 1:length(c)
            uᵢ = H.coeffs[i]
            ctwᵢ = H.complex_t_weights[i]
            c[i] = uᵢ * ctwᵢ
        end
    end

    c
end

function ModelKit.taylor!(tc, H::ToricTEmbedding, tṫ)
    t, ṫ = tṫ

    if iszero(t)
        for i = 1:length(tc)
            wᵢ = H.weights[i]
            uᵢ = H.coeffs[i]
            if wᵢ < 1e-12
                tc[i, 1] = uᵢ
                tc[i, 2] = tc[i, 3] = tc[i, 4] = tc[i, 5] = 0.0
            elseif isapprox(wᵢ, 1)
                tc[i, 2] = ṫ * uᵢ
                tc[i, 1] = tc[i, 3] = tc[i, 4] = tc[i, 5] = 0.0
            else
                tc[i, 1] = tc[i, 2] = tc[i, 3] = tc[i, 4] = tc[i, 5] = 0.0
            end
        end
    else
        evaluate_weights!(H, t)
        t⁻¹ = inv(t)
        @inbounds for i = 1:length(tc)
            wᵢ = H.weights[i]
            uᵢ = H.coeffs[i]
            twᵢ = H.complex_t_weights[i]
            tc[i, 1] = uᵢ * twᵢ
            tw1 = ṫ * wᵢ * twᵢ * t⁻¹
            tc[i, 2] = uᵢ * tw1
            tw2 = 0.5 * (wᵢ - 1) * ṫ * tw1 * t⁻¹
            tc[i, 3] = uᵢ * tw2
            tw3 = (wᵢ - 2) * ṫ * tw2 * t⁻¹ / 3
            tc[i, 4] = uᵢ * tw3
            tw4 = 0.25 * (wᵢ - 3) * ṫ * tw3 * t⁻¹
            tc[i, 5] = uᵢ * tw4
        end

    end

    tc
end
