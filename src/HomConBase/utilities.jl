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

"""
    totaldegree(F[, unit_roots=false])

Returns `G`, `solutions` where `F` is a system with the equations ``w^{d_i}-b_i=0`` where ``w`` is a ``d_i``-th unit root with ``d_i = \text{deg}F_i``
If `unit_roots=true` then ``b_i=1`` otherwise ``b_i`` is uniformly random drawn out of the interval [0,1].
"""
function totaldegree(G::PolySystem{T}; kwargs...) where {T<:Number}
    if nvariables(G) != length(G)
        return error("In order to create a total degree start system your input system needs " *
            "to have the same number of variables as number of equations. Currently your system has " *
            "$(nvariables(G)) variables and $(length(G)) equations.")
    end
   totaldegree(promote_type(T, Complex128), variables(G), degrees(G); kwargs...)
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
    totaldegree_helper(variables, degrees, b)
end

function totaldegree_helper(vars::Vector{Symbol}, degrees::Vector{Int}, b::Vector{T}) where {T<:Complex}
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
    solutions = vec(map(collect, Base.product(sol_preps...)))

    F, solutions
end


"""
    randomsystem([T<:Complex=Complex128,] nequations, nvars; mindegree=3, maxdegree=8, rng=N(0,1))

Creates a random polynomial system of `nequations` equations with `nvars` vars. Each equation has a degree uniformly drawn
from `[mindegree, maxdegree]`. The coeffiencts are drawn from the given rng. It defaults to N(0,1) where the real and
imaginary part are drawn independently.

    randomsystem([T<:Complex=Complex128,] degrees::Vector{Int}, vars::Vector{Symbol}; rng=N(0,1))

Create a random polynomial system with the given `degrees`.
"""
function randomsystem(::Type{Complex{T}}, nequations::Int, nvars::Int; mindegree=3, maxdegree=8, rng=Base.Random.GLOBAL_RNG) where {T}
    degrees = map(_ -> rand(mindegree:maxdegree), 1:nequations)
    vars = map(i -> Symbol("x_$(i)"), 1:nvars)
    randomsystem(Complex{T}, degrees, vars; rng=rng)
end
function randomsystem(nequations::Int, nvars::Int; kwargs...)
    randomsystem(Complex128, nequations, nvars; kwargs...)
end

function randomsystem{T}(::Type{Complex{T}}, degrees::Vector{Int}, vars::Vector{Symbol}; rng=Base.Random.GLOBAL_RNG)
    polys = map(degrees) do degree
        exponents = create_exponents(degree, length(vars))
        coeffs = map(complex, randn(rng, T, length(exponents)), randn(rng, T, length(exponents)))
        Poly(hcat(exponents...), coeffs)
    end
    PolySystem(polys, vars)
end

function randomsystem(degrees::Vector{Int}, vars::Vector{Symbol}; kwargs...)
    randomsystem(Complex128, degrees, vars; kwargs...)
end

function exponents_helper(curr_sum::Int, target_sum::Int, remaining_elements::Int)::Vector{Vector{Int}}
    if remaining_elements == 0
        return [[]]
    end
    if curr_sum == target_sum
        return [zeros(Int, remaining_elements)]
    end
    if remaining_elements == 1
        return map(x-> [x], 0:(target_sum - curr_sum))
    end

    results = []
    for x=0:(target_sum-curr_sum)
        remaining_results = exponents_helper(curr_sum + x, target_sum, remaining_elements - 1)
        append!(results, map(xs -> [x; xs], remaining_results))
    end
    results
end
create_exponents(total_degree, nvars) = exponents_helper(0, total_degree, nvars)
