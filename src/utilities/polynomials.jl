export is_homogeneous, homogenize, uniquevar, Composition, expand, compose, validate,
		precondition, normalize_coefficients, linear_system

const MPPoly{T} = MP.AbstractPolynomialLike{T}
const MPPolys{T} = Vector{<:MP.AbstractPolynomialLike{T}}
const WeightedVariable = Tuple{<:MP.AbstractVariable, Int}

##############
# COMPOSITION
##############

abstract type AbstractComposition end

"""
    Composition

A `Composition` is a composition of polynomial systems. This is the result of [`compose`](@ref).
"""
mutable struct Composition{T<:MP.AbstractPolynomialLike} <: AbstractComposition
    polys::Vector{Vector{T}} # polys = [f1, f2, f3] -> f1 ∘ f2 ∘ f3
end

Base.length(C::Composition) = length(C.polys)
Base.:(==)(C::Composition, D::Composition) = C.polys == D.polys

"""
    compose(g, f)::Composition

Compose the polynomial systems `g` and `f`.
You can also use the infix operator `∘` (written by \\circ).

```julia-repl
julia> @polyvar a b c x y z;
julia> g = [a * b * c];
julia> f = [x+y, y + z, x + z];
julia> expand(compose(g, f))
1-element Array{DynamicPolynomials.Polynomial{true,Int64},1}:
 x²y + x²z + xy² + 2xyz + xz² + y²z + yz²
```
"""
compose(g::Composition, f::Composition) = Composition([g.polys..., f.polys...])
compose(g::Composition, f::Composition, h::Composition...) = compose(Composition([g.polys..., f.polys...]), h...)
compose(g::Composition, fs::MPPolys...) = Composition([g.polys..., fs...])
compose(g::MPPolys, f::Composition) = Composition([g, f.polys...])
compose(g::MPPolys, f::Composition, h::Composition...) = compose(Composition([g, f.polys...]), h...)
compose(g::MPPolys, fs::MPPolys...) = Composition([g, fs...])

import Base: ∘
∘(g::Union{Composition, <:MPPolys}, f::MPPolys) = compose(g, f)
∘(g::Union{Composition, <:MPPolys}, f::MPPolys...) = compose(g, f...)
∘(g::Union{Composition, <:MPPolys}, f::Composition...) = compose(g, f...)

"""
    scale(C::Composition, λ)

Scale a composition by λ.
"""
scale(C::Composition, λ) = Composition([[C.polys[1] .* λ]; C.polys[2:end]])

function subs_into_composition(C::Composition, args...)
	p = map(pᵢ -> MP.subs(pᵢ, args...), C.polys[end])
	Composition([C.polys[1:end-1]; [p]])
end

"""
    expand(C::Composition; parameters=nothing)

Expand the composition to the polynomial system is originally represents.

## Example
```julia-repl
julia> @polyvar a b c x y z;
julia> f = [a * b * c];
julia> g = [x+y, y + z, x + z];
julia> expand(f ∘ g)
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


#############################
# HomogenizationInformation
#############################

"""
	HomogenizationInformation(;homvar=nothing, homvars=nothing, variable_groups=nothing)

This captures the information about the desired homogenization of the problem.
`homvar(s) = nothing` indicates that the problem should be homogenized.
`variable_groups` indicate the desired variable grouping for multi-homogeneous homotopies.

## Fields
* `homvars`
* `vargroups`
"""
struct HomogenizationInformation{N, T}
	# one of the two fields is always not Nothing
	homvars::Union{Nothing, NTuple{N, T}}
	vargroups::Union{Nothing, NTuple{N, Vector{T}}}
end

function HomogenizationInformation(;homvar=nothing, homvars=nothing, variable_groups=nothing)
	if homvars === nothing && homvar !== nothing
		homvars = (homvar,)
		if variable_groups !== nothing && length(variable_groups) > 1
			error("You provided more than 1 variables group but only declared `homvar=$homvar`. Use instead `homvars=...`.")
		end
	end

	if variable_groups === nothing && homvars === nothing
		return nothing
	elseif variable_groups === nothing && length(homvars) > 1
		error("Currently we cannot find variable groups just from given `homvars`. Please provide explicit variable groups using `variable_groups=...`.")
	elseif variable_groups === nothing
		HomogenizationInformation(homvars, nothing)
	else
		HomogenizationInformation(homvars, tuple(collect.(variable_groups)...))
	end
end

function add_variable_groups(HI::HomogenizationInformation, F::MPPolys; parameters=nothing)
	if HI.vargroups === nothing
		if length(HI.homvars) == 1
			vargroups = (variables(F; parameters=parameters),)
			return HomogenizationInformation(HI.homvars, vargroups)
		else
			error("Cannot add variable groups")
		end
	end
	HI
end
add_variable_groups(::Nothing, F; parameters=nothing) = nothing

homvars(H::HomogenizationInformation) = H.homvars
homvars(::Nothing) = nothing

"""
    VariableGroups{N, Homvars, V}

A `VariableGroups` stores an `NTuple` of indices mapping to the indices of the original variables
of type `V`. `Homvars` is true if there have been homogenization variables declared.
"""
struct VariableGroups{N, Homvars, V<:Union{MP.AbstractVariable, Int}}
	variables::Vector{V}
    groups::NTuple{N, Vector{Int}}
	pull_back_mapping::Vector{Tuple{Int,Int}}
end

function VariableGroups(variables, groups::NTuple{N, Vector{Int}}, dedicated_homvars::Bool) where {N}
	hom_vars = last.(groups)
	elements = vcat(groups...)
	pull_back_mapping = Tuple{Int,Int}[]
	proj_homvar_indices = cumsum([length(g) for g in groups])
	for i in 1:sum(length, groups)
		i ∉ hom_vars || continue
		l = 0
		for (gᵢ, group) in enumerate(groups), k in group
			l += 1
			if i == k
				push!(pull_back_mapping, (l, proj_homvar_indices[gᵢ]))
				break
			end
		end
	end
	VariableGroups{N, dedicated_homvars, eltype(variables)}(variables, groups, pull_back_mapping)
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
      VariableGroups(all_variables, groups, homvars !== nothing)
end
function VariableGroups(all_variables, homvar_info::HomogenizationInformation)
	if homvar_info.vargroups === nothing
		VariableGroups(all_variables, (collect(all_variables),), homvar_info.homvars)
	else
		VariableGroups(all_variables, homvar_info.vargroups, homvar_info.homvars)
	end
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
function VariableGroups(nvariables::Int, hominfo::HomogenizationInformation)
	if hominfo.vargroups === nothing
		VariableGroups(collect(1:nvariables), (collect(1:nvariables),), hominfo.homvars)
	else
		VariableGroups(collect(1:nvariables), hominfo.vargroups, hominfo.homvars)
	end
end

"""
	has_dedicated_homvars(variable_groups)

Returns `true` if there have been homogenization variables declared.
"""
has_dedicated_homvars(groups::VariableGroups{N, Homvars}) where {N, Homvars} = Homvars

"""
	projective_dims(variable_groups)

Returns the projective dimension of each variable group.
"""
function projective_dims(groups::VariableGroups)
	length.(groups.groups) .- 1
end

nvariables(G::VariableGroups) = sum(length.(G.groups))
Base.length(::VariableGroups{N}) where N = N

"""
	variable_groups(VG::VariableGroups)

Returns the variable groups.
"""
variable_groups(VG::VariableGroups) = map(g -> variables(VG)[g], VG.groups)

"""
	flattened_groups(VG::VariableGroups)

Returns the variable groups.
"""
flattened_variable_groups(VG::VariableGroups) = vcat(variable_groups(VG)...)

"""
	variables(VG::VariableGroups)

Returns all variables in their original order.
"""
variables(VG::VariableGroups) = VG.variables

"""
	ngroups(VG::VariableGroups)

Group the given variables in their corresponding groups.
"""
ngroups(VG::VariableGroups{M}) where M = M


function embed_projective(VG::VariableGroups{M}, v::AbstractVector) where {M}
	dims = projective_dims(VG)
	x = PVector(zeros(eltype(v), sum(dims) + M), dims)
	embed_projective!(x, VG, v)
end

function embed_projective!(x::PVector{<:Number, M}, VG::VariableGroups{M}, v::PVector{<:Number, M}) where {M}
	x.data .= v.data
	x
end
function embed_projective!(x::PVector{<:Number, M}, VG::VariableGroups{M}, v::AbstractVector) where {M}
	if sum(projective_dims(VG)) == length(v)
		k = 1
		for group in VG.groups
			for i in 1:(length(group)-1)
				x[k] = v[group[i]]
				k += 1
			end
			x[k] = one(eltype(x))
			k += 1
		end
	else
		k = 1
		for group in VG.groups
			for i in group
				x[k] = v[i]
				k += 1
			end
		end
	end
	x
end


##############
# POLYNOMIALS
#############
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
    is_homogeneous(f::MP.AbstractPolynomialLike)

Checks whether `f` is homogeneous.

    is_homogeneous(f::MP.AbstractPolynomialLike, vars)

Checks whether `f` is homogeneous in the variables `vars` with possible weights.
"""
is_homogeneous(f::MP.AbstractPolynomialLike) = MP.mindegree(f) == MP.maxdegree(f)
function is_homogeneous(f::MP.AbstractPolynomialLike, variables)
    d_min, d_max = minmaxdegree(f, variables)
    d_min == d_max
end

"""
    is_homogeneous(F::Vector{MP.AbstractPolynomialLike}, variables)

Checks whether each polynomial in `F` is homogeneous in the variables `variables`.
"""
function is_homogeneous(F::MPPolys, variables)
    all(f -> is_homogeneous(f, variables), F)
end
function is_homogeneous(F::MPPolys, ::Nothing=nothing; parameters=nothing)
    if parameters !== nothing
        is_homogeneous(F, variables(F, parameters=parameters))
    else
        all(is_homogeneous, F)
    end
end
function is_homogeneous(F::MPPolys, hominfo::HomogenizationInformation; parameters=nothing)
	if hominfo.vargroups !== nothing
		all(vars -> is_homogeneous(F, vars), hominfo.vargroups)
	elseif parameters !== nothing
        is_homogeneous(F, variables(F, parameters=parameters))
    else
        all(is_homogeneous, F)
    end
end
function is_homogeneous(C::Composition, hominfo=nothing; kwargs...)
    homogeneous_degrees_helper(C; kwargs...) !== nothing
end
function is_homogeneous(C::Composition, hominfo::HomogenizationInformation; kwargs...)
	if hominfo.vargroups && length(hominfo.vargroups) > 1
		error("`variable_groups` and compositions are currently not supported.")
	end

    homogeneous_degrees_helper(C; kwargs...) !== nothing
end

function homogeneous_degrees_helper(C::Composition; parameters=nothing, weights=nothing)
    for f in reverse(C.polys)
        weights = homogeneous_degrees_helper(f, parameters=parameters, weights=weights)
        weights === nothing && return nothing
    end
    weights
end
function homogeneous_degrees_helper(F::MPPolys; parameters=nothing, weights=nothing)
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
function degree(term::MP.AbstractTermLike, variables::NTuple{N,<:MP.AbstractVariable}) where {N}
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
    for term in MP.terms(f)
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
	multidegrees(F, variable_groups)

Computes the multi-degrees of the polynomials of `F` with respect to the given
variable groups. Returns an integer matrix, each column is the degree with respect
to the corresponding variable group and each row is one of the polynomials in the system.

## Example
```
julia> @polyvar x y;
julia> multidegrees([x*y-2, x^2-4], ([x], [y]))
2×2 Array{Int64,2}:
 1  2
 1  0
 ```
"""
multidegrees(F, vargroups) = [minmaxdegree(f, vars)[2] for vars in vargroups, f in F]
multidegrees(F, VG::VariableGroups) = multidegrees(F, variable_groups(VG))


"""
	support_coefficients(F::MPPolys, vars=variables(F))

Returns a tuple `(A, c)` where `A` is the vector of supports of `F` and `c` is the vector
of corresponding coefficients.
"""
function support_coefficients(F::MPPolys, vars=variables(F))
	terms₁ = MP.terms(F[1])
    c₁ = float.(MP.coefficient.(terms₁))
    A₁ = [Int32(MP.degree(term, var)) for var in vars, term in terms₁]
	A = [A₁]
	c = [c₁]
	for i in 2:length(F)
		termsᵢ = MP.terms(F[i])
		cᵢ = float.(MP.coefficient.(termsᵢ))
	   	Aᵢ = [Int32(MP.degree(term, var)) for var in vars, term in termsᵢ]
		push!(A, Aᵢ)
		push!(c, cᵢ)
	end
	A, c
end

"""
    compute_numerically_degrees(F::AbstractSystem)

This tries to compute numerically the degrees of `F`. If successfull
this returns a vector containing the degrees. If not `nothign` is returned.
"""
function compute_numerically_degrees(F)
    n, N = size(F)
    # The number of variables match, but it still cannot be homogeneous.
    # We evaluate the system with y:=rand(N) and 2y. If homogeneous then the output
    # scales accordingly to the degrees which we can obtain by taking logarithms.
    x = PVector(rand(ComplexF64, N))
    system_cache = cache(F, x)
    y = evaluate(F, x, system_cache)
    LinearAlgebra.rmul!(x, 2)
    y2 = evaluate(F, x, system_cache)

    degrees = Int[]
	for (y2ᵢ, yᵢ) in zip(y2, y)
        # y2ᵢ = 2^dᵢ yᵢ
        float_dᵢ = log2(abs(y2ᵢ / yᵢ))
        dᵢ = round(Int, float_dᵢ)
        if abs(dᵢ - float_dᵢ) > 1e-10
			return nothing
        end
        push!(degrees, dᵢ)
    end
    degrees
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
function homogenize(f::MP.AbstractPolynomialLike, hominfo::HomogenizationInformation)
	g = f
	for (vargroup, homvar) in zip(hominfo.vargroups, hominfo.homvars)
		g = homogenize(g, vargroup, homvar)
	end
	g
end

function homogenize(F::MPPolys, var=uniquevar(F); parameters=nothing)
    homogenize(F, variables(F, parameters=parameters), var)
end
function homogenize(F::MPPolys, variables::Vector, var::MP.AbstractVariable=uniquevar(F))
    map(f -> homogenize(f, variables, var), F)
end
function homogenize(F::MPPolys, hominfo::HomogenizationInformation)
	homogenize.(F, Ref(hominfo))
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
    has_constant_polynomial(F::Vector{<:MP.AbstractPolynomialLike})

Returns true if the system `F` contains a constant polynomial.
"""
function has_constant_polynomial(F::MPPolys)
    for f in F
		MP.nterms(f) == 1 && MP.isconstant(first(MP.terms(f))) && return true
	end
	false
end
has_constant_polynomial(C::Composition) = has_constant_polynomial(C.polys[1])

permute(F::MPPolys, perm) = F[perm]
permute(C::Composition, perm) = Composition([[C.polys[1][perm]]; C.polys[2:end]])

"""
    check_zero_dimensional(F::Vector{<:MP.AbstractPolynomial})

Check that the given polynomial system can have zero dimensional components.
"""
function check_zero_dimensional(F::Union{MPPolys, Composition})
    N = nvariables(F)
    n = npolynomials(F)

    if n ≥ N || (n == N - 1 && is_homogeneous(F))
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
function homogenize_if_necessary(F::Union{MPPolys, Composition}, hominfo::Union{Nothing, HomogenizationInformation}; parameters=nothing)
    vars = variables(F, parameters=parameters)

	# This fills in the simple variable group (allvars,)
	hominfo = add_variable_groups(hominfo, F; parameters=parameters)

    if is_homogeneous(F, hominfo; parameters=parameters)
		vargroups = VariableGroups(vars, hominfo)
		F, vargroups, homvars(hominfo)
    elseif hominfo === nothing
		homvar = uniquevar(F)
        push!(vars, homvar)

        F′ = homogenize(F, homvar; parameters=parameters)
        vargroups = VariableGroups(vars, homvar)

		F′, vargroups, (homvar,)
	else
        if hominfo.homvars !== nothing
            error("Input system is not homogeneous although `homvars=$(hominfo.homvars)` was passed.")
        end

		# We create a new variable for each variable group to homogenize the system
		new_homvars = map(_ -> uniquevar(F), hominfo.vargroups)
		push!(vars, new_homvars...)
		new_vargroups = map(hominfo.vargroups, new_homvars) do group, v
			vcat(group, v)
		end
		hominfo = HomogenizationInformation(new_homvars, new_vargroups)

        F′ = homogenize(F, hominfo)
		vargroups = VariableGroups(vars, hominfo)

		F′, vargroups, new_homvars
    end
end

"""
    classify_system(F, vargroups::VariableGroups; affine_tracking=false)

Returns a symbol indicating whether `F` is `:square`, `:overdetermined` or `:underdetermined`.
"""
function classify_system(F, vargroups::VariableGroups; affine_tracking=false)
	n = npolynomials(F) - nvariables(vargroups)
	if !affine_tracking
		n += ngroups(vargroups)
	end

    if n == 0
        :square
    elseif n > 0
        :overdetermined
    else
        :underdetermined
    end
end

exponent(term::MP.AbstractTermLike, vars) = [MP.degree(term, v) for v in vars]

"""
    weyldot(f::Polynomial, g::Polynomial)

Compute the [Bombieri-Weyl dot product](https://en.wikipedia.org/wiki/Bombieri_norm).
Note that this is only properly defined if `f` and `g` are homogeneous.

    weyldot(f::Vector{Polynomial}, g::Vector{Polynomial})

Compute the dot product for vectors of polynomials.
"""
function weyldot(f::MP.AbstractPolynomialLike{T}, g::MP.AbstractPolynomialLike{S}, vars=variables([f, g])) where {T,S}
    if f === g
        return sum(f) do term
            abs2(MP.coefficient(term)) / multinomial(exponent(term, vars))
        end
    end
    result = zero(promote_type(T,S))
    for term_f in f
        c_f = MP.coefficient(term_f)
        exp_f = exponent(term_f, vars)
        for term_g in g
            c_g = MP.coefficient(term_g)
            exp_g = exponent(term_g, vars)
            if exp_f == exp_g
                result += (c_f * conj(c_g)) / multinomial(exp_f)
                break
            end
        end
    end
    result
end
function weyldot(F::MPPolys, G::MPPolys, vars=variables(F))
    sum(fg -> weyldot(fg..., vars), zip(F, G))
end

"""
    weylnorm(F, vars=variables(F))

Compute the [Bombieri-Weyl norm](https://en.wikipedia.org/wiki/Bombieri_norm).
Note that this is only properly defined if `f` is homogeneous.
"""
weylnorm(F::MPPolys, vars=variables(F)) = √(sum(f -> weyldot(f, f, vars), F))
weylnorm(f::MPPoly, vars=variables(f)) = √(weyldot(f, f, vars))

"Computes the multinomial coefficient (|k| \\over k)"
function multinomial(k::Vector{Int})
    s = 0
    result = 1
    for i in k
        s += i
        result *= binomial(s, i)
    end
    result
end

"""
	precondition(f::Union{MPPolys,Composition}, variables=MP.variables(f))

Rescale the equations and variables of `f` following the heuristic described in [1].
Returns the scaled system and a vector of scaling factors of the variables.

[1]: HOM4PS-2.0: a software package for solving polynomial systems by the polyhedral homotopy continuation method,
	Lee, T.L., Li, T.Y. & Tsai, C.H. Computing (2008) 83: 109.
	https://doi.org/10.1007/s00607-008-0015-6
"""
function precondition(f::Composition, variables=MP.variables(f))
	vars, equations = preconditioning_factors(expand(f), variables)
	y = vars .* variables
	scale(subs_into_composition(f, variables => y), equations), vars
end

function precondition(f::MPPolys, variables=MP.variables(f))
	vars, equations = preconditioning_factors(f, variables)
	y = vars .* variables
	map(equations, f) do c, fᵢ
		c * MP.subs(fᵢ, variables => y)
	end, vars
end
function preconditioning_factors(f::MPPolys, vars)
	nvars = length(vars)
	m = sum(nᵢ -> nᵢ + div(nᵢ * (nᵢ - 1), 2), MP.nterms.(f))
	n = nvars + length(f)
	A, b = zeros(Int, m, n), zeros(Float64, m)

	row = 1
	for (i, fᵢ) in enumerate(f)
		termsᵢ = MP.terms(fᵢ)
		ntermsᵢ = length(termsᵢ)
		blockᵢ_start = row - 1
		for t in termsᵢ
			for (j, vⱼ) in enumerate(vars)
				A[row, j] = MP.degree(t, vⱼ)
			end
			A[row, nvars + i] = 1
			b[row] = -log10(convert(Float64, abs(MP.coefficient(t))))
			row += 1
		end

		for k in 1:ntermsᵢ, l in (k+1):ntermsᵢ
			for j in 1:nvars
				A[row, j] = A[blockᵢ_start+k, j] - A[blockᵢ_start+l, j]
			end
			b[row] = b[blockᵢ_start+k] - b[blockᵢ_start+l]
			row += 1
		end
	end
	c = ones(n)
	try
		c = LinearAlgebra.qr(A, Val(true)) \ b
	catch e
		if !isa(e, LinearAlgebra.LAPACKException)
			rethrow(e)
		end
	end
	(variables = exp10.(c[1:nvars]), equations = exp10.(c[nvars+1:end]))
end


"""
	normalize_coefficients(f::Union{MPPolys,Composition})

Rescale the equations of `f` such that each coefficient has absolute value between 0 and 1.
"""
function normalize_coefficients(f::Composition)
  	Λ = maximum.(abs, MP.coefficients.(expand(f)))
	scale(f, Λ)
end
function normalize_coefficients(f::MPPolys)
	map(f) do fᵢ
	    cmax = maximum(abs, MP.coefficients(fᵢ))
	    map(MP.terms(fᵢ)) do t
	        c = MP.coefficient(t)
	        m = MP.monomial(t)
	        c / cmax * m
	    end |> MP.polynomial
	end
end


"""
    linear_system(f::Vector{<:MP.AbstractPolynomialLike})

Given a polynomial system which represents a linear system ``Ax=b`` return
`A` and `b`. If `f ` is not a linear system `nothing` is returned.

### Example
```
julia> @polyvar x y z;
julia> f = [2x+3y+z+5, -1x+2y+4z-2];
julia> A, b = linear_system(f);
julia> A
2×3 Array{Int64,2}:
  2  3  1
 -1  2  4

julia> b
2-element Array{Int64,1}:
 -5
  2
```
"""
function linear_system(f::Vector{<:MP.AbstractPolynomialLike}, vars = MP.variables(f))
    n = length(vars)
    A = zeros(MP.coefficienttype(f[1]), length(f), n)
    b = zeros(eltype(A), length(f))
    for (i, fᵢ) in enumerate(f)
        for t in MP.terms(fᵢ)
            constant = true
            for (j, v) in enumerate(vars)
                d = MP.degree(t, v)
                d ≤ 1 || return nothing

                if d == 1
                    A[i,j] = MP.coefficient(t)
                    constant = false
                    break
                end
            end
            if constant
                b[i] = -MP.coefficient(t)
            end
        end
    end
    A, b
end
