export VariableGroups, nvariables, variables, ishomogenous, uniquevar, homogenize, projective_dims,
	remove_zeros!, Composition

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


##############
# COMPOSITION
##############

abstract type AbstractComposition end
mutable struct Composition{T<:MP.AbstractPolynomialLike} <: AbstractComposition
    polys::Vector{Vector{T}} # polys = [f1, f2, f3] -> f1 ∘ f2 ∘ f3
end

Base.length(C::Composition) = length(C.polys)

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
function variables(polys::MPPolys; parameters=nothing)
	variables = MP.variables(polys)
    if parameters !== nothing
        setdiff!(variables, parameters)
    end
	variables
end
variables(C::Composition; kwargs...) = variables(C.polys[end]; kwargs...)

"""
    nvariables(polys; parameters=nothing)

Returns the number of variables occuring in `polys`.
"""
function nvariables(F::Union{Composition, MPPolys}; parameters=nothing)
	length(variables(F, parameters=parameters))
end


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
	allvars = variables(F, parameters=parameters)
	vars = weights === nothing ? allvars : zip(allvars, weights)
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
    uniquevar(f::MP.AbstractPolynomialLike, tag=:x0)
    uniquevar(F::MPPolys, tag=:x0)

Creates a unique variable.
"""
uniquevar(f::MP.AbstractPolynomialLike, tag=:x0) = MP.similarvariable(f, gensym(tag))
uniquevar(F::MPPolys, tag=:x0) = uniquevar(F[1], tag)

"""
    homogenize(f::MP.AbstractPolynomial, variable=uniquevar(f))

Homogenize the polynomial `f` by using the given variable `variable`.

    homogenize(F::Vector{<:MP.AbstractPolynomial}, variable=uniquevar(F))

Homogenize each polynomial in `F` by using the given variable `variable`.
"""
function homogenize(f::MP.AbstractPolynomialLike, var=uniquevar(f))
    d = MP.maxdegree(f)
    MP.polynomial(map(t -> var^(d - MP.degree(t)) * t, MP.terms(f)))
end
function homogenize(F::MPPolys, var=uniquevar(F); parameters=nothing)
    if parameters !== nothing
        homogenize(F, setdiff(MP.variables(F), parameters), var)
    else
        homogenize.(F, Ref(var))
    end
end

"""
    homogenize(f::MP.AbstractPolynomial, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(f))

Homogenize the variables `v` in the polynomial `f` by using the given variable `variable`.

    homogenize(F::Vector{<:MP.AbstractPolynomial}, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(F))

Homogenize the variables `v` in each polynomial in `F` by using the given variable `variable`.
"""
function homogenize(f::MP.AbstractPolynomialLike, variables::Vector, var::MP.AbstractVariable=uniquevar(f))
    _, d_max = minmaxdegree(f, variables)
    MP.polynomial(map(f) do t
        d = degree(t, variables)
        var^(d_max - d)*t
    end)
end
function homogenize(F::MPPolys, variables::Vector, var::MP.AbstractVariable=uniquevar(F))
    map(f -> homogenize(f, variables, var), F)
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
    filter!(!iszero, c.g)
    C
end

"""
    check_zero_dimensional(F::Vector{<:MP.AbstractPolynomial})

Check that the given polynomial system can have zero dimensional components.
"""
function check_zero_dimensional(F::Vector{<:MP.AbstractPolynomial})
    N = MP.nvariables(F)
    n = length(F)

    if n ≥ N || (n == N - 1 && ishomogenous(F))
        return nothing
    end
    error("The input system will not result in a finite number of solutions.")
end


const overdetermined_error_msg = """
The input system is overdetermined. Therefore it is necessary to provide an explicit start system.
See
    https://www.JuliaHomotopyContinuation.org/guides/latest/overdetermined_tracking/
for details.
"""

"""
    homogenize_if_necessary(F::MPPolys)

Homogenizes the system `F` if necessary and returns the (new) system `F` its variables
and a subtype of [`AbstractHomogenization`] indicating whether it was homegenized.
If it was homogenized and no then the new variable is the **first one**.
"""
function homogenize_if_necessary(F::MPPolys; homvar=nothing, parameters=nothing)
    variables = MP.variables(F)
    if parameters !== nothing
        variables = setdiff(variables, parameters)
    end

    n, N = length(F), length(variables)
    if ishomogenous(F; parameters=parameters)
        # N = n+1 is the only valid size configuration
        if n + 1 > N
            error(overdetermined_error_msg)
        end
		vargroups = VariableGroups(variables, homvar)
		F, variables[vcat(vargroups.groups...)], vargroups
    else
        if homvar !== nothing
            error("Input system is not homogenous although `homvar` was passed.")
        end
        # We create a new variable to homogenize the system
        homvar = uniquevar(F)
        push!(variables, homvar)
        sort!(variables, rev=true)

        F′ = homogenize(F, homvar; parameters=parameters)
		vargroups = VariableGroups(variables, homvar)

		F′, variables[vcat(vargroups.groups...)], vargroups
    end
end
