abstract type TEmbedding end

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


struct CacheTEmbedding{E<:TEmbedding}
    embedding::E
    t_cache::Ref{ComplexF64}
    tṫ_cache::Ref{Tuple{ComplexF64,ComplexF64}}
end
CacheTEmbedding(E::TEmbedding) =
    CacheTEmbedding(E, Ref(ComplexF64(NaN)), Ref((ComplexF64(NaN), ComplexF64(NaN))))

function ModelKit.evaluate!(u, H::CacheTEmbedding, t)
    # t == H.t_cache[] && return u
    H.t_cache[] = t
    ModelKit.evaluate!(u, H.embedding, t)
end

function ModelKit.taylor!(u, H::CacheTEmbedding, tṫ::Tuple{<:Number,<:Number})
    # tṫ == H.tṫ_cache[] && return u
    H.tṫ_cache[] = tṫ
    ModelKit.taylor!(u, H.embedding, tṫ)
    u
end

cached(E::TEmbedding) = CacheTEmbedding(E)
function update!(H::CacheTEmbedding, arg1)
    H.t_cache[] = NaN
    H.tṫ_cache[] = (NaN, NaN)
    update!(H.embedding, arg1)
end
function update!(H::CacheTEmbedding, arg1, arg2)
    H.t_cache[] = NaN
    H.tṫ_cache[] = (NaN, NaN)
    update!(H.embedding, arg1, arg2)
end