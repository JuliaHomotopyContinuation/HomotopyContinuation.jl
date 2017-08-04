"""
    projective(x::Vector)

Embeds a vector `x` into projective space via `[1; x]`;
"""
projective(x::Vector{T}) where {T<:Number} = [one(T); x]

"""
    affine(x::Vector)

Maps a projective vector `x` into affine space via `x[2:end]/x[1]`;
"""
affine(x::Vector{T}) where {T<:Number} = x[2:end]/x[1]
