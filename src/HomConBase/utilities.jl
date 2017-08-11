"""
    projective(x, varindex)

Embeds a vector `x` into projective space assuming the "projective variable" has index `varindex`.
"""
function projective(x::Vector{T}, index::Int) where {T<:Number}
    res = copy(x)
    insert!(res, index, one(T))
    res
end

"""
    affine(x, varindex)

Maps a projective vector `x` into affine space assuming the "projective variable" has index `varindex`.
"""
affine(v::Vector{T}, i::Int) where {T<:Number} = [ v[1:i-1] ; v[i+1:end] ] / v[i]
