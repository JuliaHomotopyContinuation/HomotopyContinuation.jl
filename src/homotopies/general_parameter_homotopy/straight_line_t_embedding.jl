
# ConvexCombination(p, q)
# The function γ(t) = t * p + (1 - t) * q 
struct StraightLineTEmbedding <: TEmbedding
    p::Vector{ComplexF64}
    q::Vector{ComplexF64}
end

function ModelKit.evaluate!(u, H::StraightLineTEmbedding, t)
    for i = 1:length(u)
        u[i] = t * H.p[i] + (1.0 - t) * H.q[i]
    end
    u
end

function ModelKit.taylor!(u, H::StraightLineTEmbedding, tṫ::Tuple{<:Number,<:Number})
    (t, ṫ) = tṫ
    @inbounds for i = 1:length(u)
        ptᵢ = t * H.p[i] + (1.0 - t) * H.q[i]
        u[i] = TruncatedTaylorSeries((ptᵢ, ṫ * (H.p[i] - H.q[i])))
    end
    u
end

function update!(H::StraightLineTEmbedding, p, q)
    H.p .= p
    H.q .= q
    H
end