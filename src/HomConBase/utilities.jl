"""
    projective(x::Vector)

Embeds a vector `x` into projective space via `[1; x]`;
"""
projective(x::Vector{T}) where T = [one(T); x] # TODO: shouldn't we have an option for a random patch?

"""
    affine(x::Vector)

Maps a projective vector `x` into affine space via `x[2:end]/x[1]`;
"""
affine(x::AbstractVector)= x[2:end]/x[1]

function total_degree_helper(vars::Vector{Symbol}, degrees::Vector{Int}, b::Vector{T}) where {T<:Complex}
    n = length(vars)


    polys = map(i -> begin
        d = degrees[i]
        exponents = zeros(Int, n, 2)
        coeffs =  zeros(T, 2)
        exponents[i, 1] = d
        coeffs[1] = one(T)
        coeffs[2] = -b[i]
        Poly(exponents, coeffs)
    end, 1:length(degrees))
    F = PolySystem(polys, vars)

    sol_preps = Vector{Vector{T}}()
    for (i, d) in enumerate(degrees)
        sol = [b[i]^(1/d) * exp(2π*im/d)^k for k in 0:d-1]
        push!(sol_preps, sol)
    end

    # TODO: Make solutions an iterator
    solutions = vec(map(collect, Base.product(sol_preps...)))

    F, solutions
end

"""
    totaldegree(T, vars, degrees[, unit_roots=false])

Returns `G`, `solutions` where `G` is a system with the equations ``w^{d_i}-b_i=0`` where ``w`` is a ``d_i``-th unit root.
If `unit_roots=true` then ``b_i=1`` otherwise ``b_i`` is uniformly random drawn out of the interval [0,1].
"""
function totaldegree(::Type{T}, variables::Vector{Symbol}, degrees::Vector{Int}; unit_roots=false ) where {T<:Complex}
    if unit_roots
        b = ones(T, length(degrees))
    else
        b = rand(T, length(degrees))
    end
    total_degree_helper(variables, degrees, b)
end
"""
    totaldegree(F[, unit_roots=false])

Returns `G`, `solutions` where `F` is a system with the equations ``w^{d_i}-b_i=0`` where ``w`` is a ``d_i``-th unit root with ``d_i = \text{deg}F_i``
If `unit_roots=true` then ``b_i=1`` otherwise ``b_i`` is uniformly random drawn out of the interval [0,1].
"""
function totaldegree(G::PolySystem{T}; kwargs...) where {T<:Number}
   totaldegree(promote_type(T, Complex128), variables(G), degrees(G); kwargs...)
end
