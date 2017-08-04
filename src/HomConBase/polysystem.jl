
# Constructors

PolySystem{T} = Vector{Poly{T}}

evaluate(F::PolySystem{T}, x::Vector{T}) where {T<:Number}= map(f -> evaluate(f, x), F)

@inline function jacobian(F::PolySystem{T}) where {T<:Number}
    DF = permutedims(hcat(gradient.(F)...), [2, 1])
    (x::Vector{T}) -> map(f -> evaluate(f, x), DF)
end

@inline is_homogenous(F::PolySystem) = all(is_homogenous, F)
@inline homogenize(F::PolySystem) = homogenize.(F)
@inline degrees(F::PolySystem) = deg.(F)
@inline nequations(F::PolySystem) = size(F, 1)
@inline nvars(F::PolySystem) = nvars(F[1])
@inline vars(F::PolySystem) = MPoly.vars(F[1])
@inline function weyl_dot(F::PolySystem{T},G::PolySystem{T}) where {T<:Complex}
    sum(tup -> weyl_dot(tup[1],tup[2]), zip(F, G))
end

@inline function weyl_norm(F::PolySystem{T}) where {T<:Complex}
    √(weyl_dot(F, F))
end

function total_degree_helper(vars::Vector{Symbol}, degrees::Vector{Int}, b::Vector{T}) where {T<:Complex}
    n = length(vars)

    F = Vector{Poly{T}}()
    for (i, d) in enumerate(degrees)
        exponents = zeros(Int, n, 2)
        coeffs =  zeros(T, 2)
        # if homogenized
        #     exponents[1, 2] = d
        #     exponents[i + 1, 1] = d
        #     coeffs[1] = one(T)
        #     coeffs[2] = -b[i]
        # else
            exponents[i, 1] = d
            coeffs[1] = one(T)
            coeffs[2] = -b[i]
        # end
        push!(F, createpoly(exponents, coeffs, vars))
    end

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
    total_degree(T, vars, degrees[, unit_roots=false])

Returns `G`, `solutions` where `G` is a system with the equations ``w^{d_i}-b_i=0`` where ``w`` is a ``d_i``-th unit root.
If `unit_roots=true` then ``b_i=1`` otherwise ``b_i`` is uniformly random drawn out of the interval [0,1].
"""
function total_degree(::Type{T}, variables::Vector{Symbol}, degrees::Vector{Int}; unit_roots=false ) where {T<:Complex}
    if unit_roots
        b = ones(T, length(degrees))
    else
        b = rand(T, length(degrees))
    end
    total_degree_helper(variables, degrees, b)
end
"""
    total_degree(F[, unit_roots=false])

Returns `G`, `solutions` where `F` is a system with the equations ``w^{d_i}-b_i=0`` where ``w`` is a ``d_i``-th unit root with ``d_i = \text{deg}F_i``
If `unit_roots=true` then ``b_i=1`` otherwise ``b_i`` is uniformly random drawn out of the interval [0,1].
"""
function total_degree(G::PolySystem{T}; kwargs...) where {T<:Complex}
   total_degree(T, vars(G), degrees(G); kwargs...)
end
