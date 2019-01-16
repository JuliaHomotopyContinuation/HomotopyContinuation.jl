export VariableGroups, nvariables, variables, npolynomials, ishomogenous, uniquevar, homogenize, projective_dims,
	remove_zeros!, Composition, maxdegrees, check_zero_dimensional, classify_homogenous_system, check_square_homogenous_system,
	hasparameters

const MPPoly = MP.AbstractPolynomialLike
const MPPolys = Vector{<:MP.AbstractPolynomialLike}
const WeightedVariable = Tuple{<:MP.AbstractVariable, Int}

"""
    VariableGroups{N}

A `VariableGroups` stores an `NTuple` of indices mapping to the indices of the original variables.
"""
struct VariableGroups{N}
    groups::NTuple{N, Vector{Int}}
    dedicated_homvars::Bool
end

"""
    VariableGroups(all_variables, variable_groups, homvars)

# Example
```julia
all_vars = @polyvar a b c x y z
variable_groups = ([c, b, a], [x, y, z])
homvars = (b, y)

vargroups = VariableGroups(all_vars, variable_groups, homvars)
vargroups.groups == ([3, 1, 2], [4, 6, 5])
vargroups.homvars == true
```
"""
function VariableGroups(all_variables,
                       variable_groups::NTuple{N, <:Vector{<:Union{Int, MP.AbstractVariable}}},
                       homvars::Union{Nothing, NTuple{N, <:Union{Int, MP.AbstractVariable}}}=nothing) where N
    groups = ntuple(N) do k
        variable_group = variable_groups[k]
        h = homvars === nothing ? 0 : homvars[k]
        group = Int[]
        homvar = 0
        for v in variable_group
            i = findfirst(w -> w == v, all_variables)
            if v == h
                homvar = i
            else
                push!(group, i)
            end
         end
         if homvar != 0
             push!(group, homvar)
         end
         group
      end
      VariableGroups(groups, homvars !== nothing)
end

function VariableGroups(all_variables, homvar::MP.AbstractVariable)
    VariableGroups(all_variables, (collect(all_variables),), (homvar,))
end
function VariableGroups(all_variables, homvar::Nothing=nothing)
    VariableGroups(all_variables, (collect(all_variables),), homvar)
end
function VariableGroups(nvariables::Int, homvar::Nothing=nothing)
    VariableGroups(collect(1:nvariables), (collect(1:nvariables),), nothing)
end
function VariableGroups(nvariables::Int, homvar::Int)
    VariableGroups(collect(1:nvariables), (collect(1:nvariables),), (homvar,))
end

"""
	projective_dims(variable_groups)

Returns the projective dimension of each variable group.
"""
projective_dims(groups::VariableGroups) = length.(groups.groups) .- 1

nvariables(G::VariableGroups) = sum(length.(G.groups))
Base.length(::VariableGroups{N}) where N = N

##############
# COMPOSITION
##############

abstract type AbstractComposition end
mutable struct Composition{T<:MP.AbstractPolynomialLike} <: AbstractComposition
    polys::Vector{Vector{T}} # polys = [f1, f2, f3] -> f1 ∘ f2 ∘ f3
end

Base.length(C::Composition) = length(C.polys)
Base.:(==)(C::Composition, D::Composition) = C.polys == D.polys

compose(g::Composition, f::Composition) = Composition([g.polys..., f.polys...])
compose(g::Composition, f::Composition, h::Composition...) = compose(Composition([g.polys..., f.polys...]), h...)
compose(g::Composition, fs::MPPolys...) = Composition([g.polys..., fs...])
compose(g::MPPolys, f::Composition) = Composition([g, f.polys...])
compose(g::MPPolys, f::Composition, h::Composition...) = compose(Composition([g, f.polys...]), h...)
compose(g::MPPolys, fs::MPPolys...) = Composition([g, fs...])

import Base: ∘
∘(g::Union{Composition, <:MPPolys}, f::Union{<:MPPolys, <:Composition}) = compose(g, f)
∘(g::Union{Composition, <:MPPolys}, f::MPPolys) = compose(g, f)
∘(g::Union{Composition, <:MPPolys}, f::MPPolys...) = compose(g, f...)
∘(g::Union{Composition, <:MPPolys}, f::Composition...) = compose(g, f...)


"""
	expand(C::Composition; parameters=nothing)

Expand the composition to the polynomial system is originally represents.

## Example
```julia-repl
julia> @polyvar a b c x y z;
julia> f = [a * b * c];
julia> g = [x+y, y + z, x + z];
julia> Utilities.expand(f ∘ g)
1-element Array{DynamicPolynomials.Polynomial{true,Int64},1}:
 x²y + x²z + xy² + 2xyz + xz² + y²z + yz²
```
"""
function expand(C::Composition; parameters=nothing)
	g = nothing
	for f in reverse(C.polys)
		if g === nothing
			g = f
			continue
		end
		vars = variables(f, parameters=parameters)
		g = map(fᵢ -> MP.subs(fᵢ, vars => g), f)
	end
	g
end

"""
	validate(C::Composition; parameters=nothing)

Validates that the composition is well defined.
"""
function validate(C::Composition; parameters=nothing)
	g = nothing
	for f in reverse(C.polys)
		if g === nothing
			g = f
			continue
		end
		nvars = length(variables(f, parameters=parameters))
		if nvars !== length(g)
			return false
		end
		g = f
	end
	true
end


"""
    variables(F; parameters=nothing)

Returns the variables occuring in `F`.
"""
function variables(polys::Union{MPPoly, MPPolys}; parameters=nothing, weights=nothing)
	variables = MP.variables(polys)
    if parameters !== nothing
        setdiff!(variables, parameters)
    end
	if weights === nothing
		variables
	else
		zip(variables, weights)
	end
end
variables(C::Composition; kwargs...) = variables(C.polys[end]; kwargs...)


"""
    hasparameters(F, parameters=nothing)

Returns `true` if the parameters occur in F.
"""
function hasparameters(polys::Union{MPPoly, MPPolys}, parameters=nothing)
	parameters === nothing && return false

	variables = MP.variables(polys)
	for p in parameters
		if p in variables
			return true
		end
	end
	false
end

"""
    nvariables(polys; parameters=nothing)

Returns the number of variables occuring in `polys`.
"""
function nvariables(F::Union{Composition, MPPolys}; parameters=nothing)
	length(variables(F, parameters=parameters))
end

"""
    npolynomials(F)

Returns the number of polynomials occuring in `F`.
"""
npolynomials(F::MPPolys) = length(F)
npolynomials(F::Composition) = length(F.polys[1])


"""
    ishomogenous(f::MP.AbstractPolynomialLike)

Checks whether `f` is homogenous.

    ishomogenous(f::MP.AbstractPolynomialLike, vars)

Checks whether `f` is homogenous in the variables `vars` with possible weights.
"""
ishomogenous(f::MP.AbstractPolynomialLike) = MP.mindegree(f) == MP.maxdegree(f)
function ishomogenous(f::MP.AbstractPolynomialLike, variables)
    d_min, d_max = minmaxdegree(f, variables)
    d_min == d_max
end

"""
    ishomogenous(F::Vector{MP.AbstractPolynomialLike}, variables)

Checks whether each polynomial in `F` is homogenous in the variables `variables`.
"""
function ishomogenous(F::MPPolys, variables)
    all(f -> ishomogenous(f, variables), F)
end
function ishomogenous(F::MPPolys; parameters=nothing)
    if parameters !== nothing
        ishomogenous(F, variables(F, parameters=parameters))
    else
        all(ishomogenous, F)
    end
end
function ishomogenous(C::Composition; kwargs...)
    homogenous_degrees_helper(C; kwargs...) !== nothing
end

function homogenous_degrees_helper(C::Composition; parameters=nothing, weights=nothing)
	for f in reverse(C.polys)
    	weights = homogenous_degrees_helper(f, parameters=parameters, weights=weights)
    	weights === nothing && return nothing
	end
	weights
end
function homogenous_degrees_helper(F::MPPolys; parameters=nothing, weights=nothing)
	vars = variables(F, parameters=parameters, weights=weights)
    degrees = Int[]
    for f in F
        mindeg, maxdeg = minmaxdegree(f, vars)
        if mindeg != maxdeg
            return nothing
        else
            push!(degrees, mindeg)
        end
    end
    degrees
end

"""
    degree(term::MP.AbstractTermLike, vars)

Compute the (weighted) degree of `f` in the variables `vars`.
"""
function degree(term::MP.AbstractTermLike, variables::Vector{<:MP.AbstractVariable})
    sum(MP.degree(term, v) for v in variables)
end
function degree(term::MP.AbstractTermLike, weighted_variables)
    deg = 0
    for (v, w) in weighted_variables
        deg += w * MP.degree(term, v)
    end
    deg
end

"""
    minmaxdegree(f::MP.AbstractPolynomialLike, variables)

Compute the minimum and maximum (weighted total) degree of `f` with respect to the given variables.
"""
function minmaxdegree(f::MP.AbstractPolynomialLike, variables)
    d_min, d_max = typemax(Int), 0
    for term in f
        d = degree(term, variables)
        d_min, d_max = min(d, d_min), max(d, d_max)
    end
    d_min, d_max
end

"""
	maxdegrees(F, parameters=nothing)

Computes the degrees of the polynomials of `F`.
"""
function maxdegrees(F::MPPolys; parameters=nothing, weights=nothing)
	vars = variables(F, parameters=parameters, weights=weights)
	last.(minmaxdegree.(F, Ref(vars)))
end
function maxdegrees(C::Composition; parameters=nothing)
	weights = nothing
	for f in reverse(C.polys)
    	weights = maxdegrees(f, parameters=parameters, weights=weights)
	end
	weights
end

"""
    uniquevar(f::MP.AbstractPolynomialLike, tag=:x0)
    uniquevar(F::MPPolys, tag=:x0)

Creates a unique variable.
"""
uniquevar(f::MP.AbstractPolynomialLike, tag=:x0) = MP.similarvariable(f, gensym(tag))
uniquevar(F::MPPolys, tag=:x0) = uniquevar(F[1], tag)
uniquevar(C::Composition, tag=:x0) = uniquevar(C.polys[1][1], tag)


"""
    homogenize(f::MP.AbstractPolynomial, variable=uniquevar(f))

Homogenize the polynomial `f` by using the given variable `variable`.

    homogenize(F::Vector{<:MP.AbstractPolynomial}, variable=uniquevar(F))

Homogenize each polynomial in `F` by using the given variable `variable`.

	homogenize(f::MP.AbstractPolynomial, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(f))

Homogenize the variables `v` in the polynomial `f` by using the given variable `variable`.

	homogenize(F::Vector{<:MP.AbstractPolynomial}, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(F))

Homogenize the variables `v` in each polynomial in `F` by using the given variable `variable`.
"""
function homogenize(f::MP.AbstractPolynomialLike, var=uniquevar(f); parameters=nothing)
	vars = variables(f; parameters=parameters)
	homogenize(f, vars, var)
end
function homogenize(f::MP.AbstractPolynomialLike, variables::Vector, var::MP.AbstractVariable=uniquevar(f))
    _, d_max = minmaxdegree(f, variables)
    MP.polynomial(map(f) do t
        d = degree(t, variables)
        var^(d_max - d)*t
    end)
end
function homogenize(F::MPPolys, var=uniquevar(F); parameters=nothing)
    homogenize(F, variables(F, parameters=parameters), var)
end
function homogenize(F::MPPolys, variables::Vector, var::MP.AbstractVariable=uniquevar(F))
    map(f -> homogenize(f, variables, var), F)
end
function homogenize(C::Composition, var::MP.AbstractVariable=uniquevar(C.polys[1]); parameters=nothing, weights=nothing)
	polys = map(length(C.polys):-1:1) do k
		f̄, weights = homogenize_degrees(C.polys[k], var; parameters=parameters, weights=weights)
		if k > 1
			push!(f̄, var)
			push!(weights, 1)
		end
    	f̄
	end
	Composition(reverse!(polys))
end


"""
	homogenize_degree(f::MPPoly, variables, var)

Homogenize `f` by using the variable `v`. Returns the homogenized polynomial
and its degree.
"""
function homogenize_degree(f::MPPoly, variables, var::MP.AbstractVariable)
    _, d = minmaxdegree(f, variables)
    MP.polynomial(map(t -> t * var^(d - degree(t, variables)), f)), d
end


"""
	homogenize_degrees(F::MPPolys, var; parameters=nothing, weights=nothing)
	homogenize_degrees(F::MPPolys, variables, var)

Homogenize the system `F` using the variable `var`. Returns the homogenized system
and the degrees of each polynomial.
"""
function homogenize_degrees(F::MPPolys, variables, var::MP.AbstractVariable)
    F̄ = similar(F)
    degrees = Int[]
    for (i,f) in enumerate(F)
        f̄, d = homogenize_degree(f, variables, var)
        F̄[i] = f̄
        push!(degrees, d)
    end

    F̄, degrees
end
function homogenize_degrees(F::MPPolys, var::MP.AbstractVariable; parameters=nothing, weights=nothing)
	allvars = variables(F, parameters=parameters)
	if weights !== nothing
		homogenize_degrees(F, zip(allvars, weights), var)
	else
		homogenize_degrees(F, allvars, var)
	end
end


"""
    remove_zeros!(F::Vector{<:MP.AbstractPolynomialLike})

Remove zero polynomials from the given system `F`.
"""
function remove_zeros!(F::MPPolys)
    filter!(!iszero, F)
    F
end
function remove_zeros!(C::Composition)
    filter!(!iszero, C.polys[1])
    C
end

"""
    check_zero_dimensional(F::Vector{<:MP.AbstractPolynomial})

Check that the given polynomial system can have zero dimensional components.
"""
function check_zero_dimensional(F::Union{MPPolys, Composition})
    N = nvariables(F)
    n = npolynomials(F)

    if n ≥ N || (n == N - 1 && ishomogenous(F))
        return nothing
    end
    error("The input system will not result in a finite number of solutions.")
end

"""
    homogenize_if_necessary(F::MPPolys)

Homogenizes the system `F` if necessary and returns the (new) system `F` its variables
and a subtype of [`AbstractHomogenization`] indicating whether it was homegenized.
If it was homogenized and no then the new variable is the **first one**.
"""
function homogenize_if_necessary(F::Union{MPPolys, Composition}; homvar=nothing, parameters=nothing)
    vars = variables(F, parameters=parameters)

    n, N = npolynomials(F), length(vars)
    if ishomogenous(F; parameters=parameters)
		vargroups = VariableGroups(vars, homvar)
		F, vars[vcat(vargroups.groups...)], vargroups, homvar
    else
        if homvar !== nothing
            error("Input system is not homogenous although `homvar=$(homvar)` was passed.")
        end
        # We create a new variable to homogenize the system
        homvar = uniquevar(F)
        push!(vars, homvar)
        sort!(vars, rev=true)

        F′ = homogenize(F, homvar; parameters=parameters)
		vargroups = VariableGroups(vars, homvar)

		F′, vars[vcat(vargroups.groups...)], vargroups, homvar
    end
end

"""
	classify_homogenous_system(F, vargroups::VariableGroups)

Returns a symbol indicating whether `F` is `:square`, `:overdetermined` or `:underdetermined`.
"""
function classify_homogenous_system(F, vargroups::VariableGroups)
	n = npolynomials(F) - (nvariables(vargroups) - length(vargroups))

	if n == 0
		:square
	elseif n > 0
		:overdetermined
	else
		:underdetermined
	end
end

const overdetermined_error_msg = """
The input system is overdetermined. Therefore it is necessary to provide an explicit start system.
See
    https://www.JuliaHomotopyContinuation.org/guides/latest/overdetermined_tracking/
for details.
"""

"""
	check_square_homogenous_system(F, vargroups::VariableGroups)

Checks whether `F` is a square polynomial system.
"""
function check_square_homogenous_system(F, vargroups::VariableGroups)
	class = classify_homogenous_system(F, vargroups)
	if class == :overdetermined
		error(overdetermined_error_msg)
	elseif class == :underdetermined
		error("Underdetermined polynomial systems are currently not supported.")
	end
	nothing
end
