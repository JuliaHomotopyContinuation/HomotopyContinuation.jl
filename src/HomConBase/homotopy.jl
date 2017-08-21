evaluate(H::AbstractHomotopy, x, t) = error("`evaluate` not implemented for $(typeof(H))<:AbstractHomotopy")
substitute(H::AbstractHomotopy, pair::Pair{Symbol,<:Number}) = error("`substitute` not implemented for $(typeof(H))<:AbstractHomotopy")
startsystem(H::AbstractHomotopy) = error("`startsystem` not implemented for $(typeof(H))<:AbstractHomotopy")
targetsystem(H::AbstractHomotopy) = error("`targetsystem` not implemented for $(typeof(H))<:AbstractHomotopy")
nvariables(H::AbstractHomotopy) = error("`nvars` not implemented for $(typeof(H))<:AbstractHomotopy")
nequations(H::AbstractHomotopy) = error("`nequations` not implemented for $(typeof(H))<:AbstractHomotopy")
degrees(H::AbstractHomotopy) = error("`degrees` not implemented for $(typeof(H))<:AbstractHomotopy")
differentiate(H::AbstractHomotopy) = error("`differentiate` not implemented for $(typeof(H))<:AbstractHomotopy")
dt(H::AbstractHomotopy) = error("`dt` not implemented for $(typeof(H))<:AbstractHomotopy")


ishomogenous(H::AbstractHomotopy) = error("`ishomogenous` not implemented for $(typeof(H))<:AbstractHomotopy")
homogenize(H::AbstractHomotopy) = error("`homogenize` not implemented for $(typeof(H))<:AbstractHomotopy")
homogenized(H::AbstractHomotopy) = error("`homogenized` not implemented for $(typeof(H))<:AbstractHomotopy")

weylnorm(H::AbstractHomotopy, t::Float64) = error("`weyl_norm` not implemented for $(typeof(H)) for AbstractHomotopy")

@inline ∂t(H::AbstractHomotopy) = dt(H)
eltype(H::AbstractHomotopy{T}) where {T<:Number} = T
