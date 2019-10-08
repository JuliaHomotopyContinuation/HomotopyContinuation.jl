function onpatch!(x::AbstractVector, v̄::PVector{T,N}) where {T,N}
    ranges = ProjectiveVectors.dimension_indices(v̄)
    for range in ranges
        λ = zero(eltype(x))
        @inbounds for i in range
            λ += v̄[i] * x[i]
        end
        λ⁻¹ = @fastmath inv(λ)
        for i in range
            x[i] *= λ⁻¹
        end
    end
    x
end

function evaluate_patch!(u, v̄::PVector{S,N}, x::PVector{T,N}) where {S,T,N}
    ranges = ProjectiveVectors.dimension_indices(v̄)
    n = length(u) - N
    for (k, range) in enumerate(ranges)
        out = -one(eltype(x))
        for i in range
            out += v̄[i] * x[i]
        end
        u[n+k] = out
    end
    nothing
end

function jacobian_patch!(U, v̄::PVector{S,N}, x::PVector) where {S,T,N}
    ranges = ProjectiveVectors.dimension_indices(v̄)
    n = size(U, 1) - N
    for j = 1:size(U, 2), i = (n+1):size(U, 1)
        U[i, j] = zero(eltype(U))
    end
    for (k, range) in enumerate(ranges)
        for j in range
            U[n+k, j] = v̄[j]
        end
    end
    nothing
end
