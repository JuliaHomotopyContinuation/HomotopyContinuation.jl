evaluate(H::AbstractHomotopy{T}, x::Vector{T}, t::Float64) where T = error("`evaluate` not implemented for $(typeof(H))<:AbstractHomotopy")
startsystem(H::AbstractHomotopy) = error("`startsystem` not implemented for $(typeof(H))<:AbstractHomotopy")
targetsystem(H::AbstractHomotopy) = error("`targetsystem` not implemented for $(typeof(H))<:AbstractHomotopy")
nvars(H::AbstractHomotopy) = error("`nvars` not implemented for $(typeof(H))<:AbstractHomotopy")
vars(H::AbstractHomotopy) = error("`vars` not implemented for $(typeof(H))<:AbstractHomotopy")
nequations(H::AbstractHomotopy) = error("`nequations` not implemented for $(typeof(H))<:AbstractHomotopy")
degrees(H::AbstractHomotopy) = error("`degrees` not implemented for $(typeof(H))<:AbstractHomotopy")
jacobian(H::AbstractHomotopy) = error("`jacobian` not implemented for $(typeof(H))<:AbstractHomotopy")
dt(H::AbstractHomotopy) = error("`dt` not implemented for $(typeof(H))<:AbstractHomotopy")
@inline ∂t(H::AbstractHomotopy) = dt(H)

is_homogenous(H::AbstractHomotopy) = error("`is_homogenous` not implemented for $(typeof(H))<:AbstractHomotopy")
homogenize(H::AbstractHomotopy, ::MP.AbstractVariable) = error("`homogenize` not implemented for $(typeof(H))<:AbstractHomotopy")

# optional
weyl_norm(H::AbstractHomotopy, t::Float64) = error("`weyl_norm` not implemented for $(typeof(H)) for AbstractHomotopy")

eltype(H::AbstractHomotopy{T}) where {T<:Number} = T