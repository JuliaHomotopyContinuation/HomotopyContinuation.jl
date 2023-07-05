struct ToricLogRealEmbedding <: TEmbedding
    # (exp(s) * c_start[i] + (1 - exp(s)) * c_target[i])* exp(w[i] * s)
    coeffs::Vector{ComplexF64}
    support::Vector{Matrix{Int}}
    weights::Vector{Float64}
    s_weights::Vector{Float64}
    s_cache_weights::Ref{Float64}
end

function ToricLogRealEmbedding(
    coeffs::Vector{ComplexF64},
    support::AbstractVector{<:AbstractMatrix},
)
    ToricLogRealEmbedding(
        coeffs,
        map(s -> convert(Matrix{Int}, s), support),
        zeros(length(coeffs)),
        zeros(Float64, length(coeffs)),
        Ref(NaN),
    )
end


function update!(
    H::ToricLogRealEmbedding,
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

    H.s_cache_weights[] = ComplexF64(NaN)

    s_min, s_max
end


function evaluate_weights!(T::ToricLogRealEmbedding, t::Number)
    s = real(t)
    # s == T.s_cache_weights[] && return T.s_weights

    u = T.s_weights
    weights = T.weights
    n = length(weights)
    for i = 1:n
        u[i] = exp(weights[i] * s)
    end

    T.s_cache_weights[] = s

    T.s_weights
end


function ModelKit.evaluate!(c, H::ToricLogRealEmbedding, t)
    evaluate_weights!(H, t)
    @inbounds for i = 1:length(c)
        uᵢ = H.coeffs[i]
        ctwᵢ = H.s_weights[i]
        c[i] = uᵢ * ctwᵢ
    end

    c
end

function ModelKit.taylor!(tc, H::ToricLogRealEmbedding, tṫ)
    t, ṫ = tṫ
    evaluate_weights!(H, t)
    ṡ = real(ṫ)
    @inbounds for i = 1:length(tc)
        wᵢ = H.weights[i]
        uᵢ = H.coeffs[i]
        swᵢ = H.s_weights[i]
        ṡ_w = ṡ * wᵢ

        sw1 = uᵢ * swᵢ
        sw2 = ṡ_w * sw1
        sw3 = ṡ_w * sw2 * 0.5
        sw4 = ṡ_w * sw3 / 3
        sw5 = ṡ_w * sw4 * 0.25

        tc[i, 1] = sw1
        tc[i, 2] = sw2
        tc[i, 3] = sw3
        tc[i, 4] = sw4
        tc[i, 5] = sw5
    end

    tc
end
