export VariableGroups, allvariables, nvariables, ishomogenous, uniquevar, homogenize, projective_dims

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


"""
    homogenize_if_necessary(F::Vector{<:MP.AbstractPolynomialLike})

Homogenizes the system `F` if necessary and returns the (new) system `F` its variables
and a subtype of [`AbstractHomogenization`] indicating whether it was homegenized.
If it was homogenized and no then the new variable is the **first one**.
"""
function homogenize_if_necessary(F::Vector{<:MP.AbstractPolynomialLike}; homvar=nothing, parameters=nothing)
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

"""
    allvariables(polys)

Returns a sorted list of all variables occuring in `polys`.
"""
function allvariables(polys::Vector{<:MP.AbstractPolynomialLike})
    sort!(union(Iterators.flatten(MP.variables.(polys))), rev=true)
end

"""
    nvariables(polys)

Returns the number of variables occuring in `polys`.
"""
nvariables(polys::Vector{<:MP.AbstractPolynomialLike}) = length(allvariables(polys))


"""
    ishomogenous(f::MP.AbstractPolynomialLike)

Checks whether `f` is homogenous.

    ishomogenous(polys::Vector{MP.AbstractPolynomialLike})

Checks whether each polynomial in `polys` is homogenous.
"""
ishomogenous(f::MP.AbstractPolynomialLike) = MP.mindegree(f) == MP.maxdegree(f)
function ishomogenous(F::Vector{<:MP.AbstractPolynomialLike}; parameters=nothing)
    if parameters !== nothing
        ishomogenous(F, setdiff(MP.variables(F), parameters))
    else
        all(ishomogenous, F)
    end
end

"""
    ishomogenous(f::MP.AbstractPolynomialLike, v::Vector{<:MP.AbstractVariable})

Checks whether `f` is homogenous in the variables `v`.

    ishomogenous(polys::Vector{<:MP.AbstractPolynomialLike}, v::Vector{<:MP.AbstractVariable})

Checks whether each polynomial in `polys` is homogenous in the variables `v`.
"""
function ishomogenous(f::MP.AbstractPolynomialLike, variables::Vector{T}) where {T<:MP.AbstractVariable}
    d_min, d_max = minmaxdegree(f, variables)
    d_min == d_max
end

function ishomogenous(F::Vector{<:MP.AbstractPolynomialLike}, variables::Vector{T}) where {T<:MP.AbstractVariable}
    all(f -> ishomogenous(f, variables), F)
end


"""
    minmaxdegree(f::MP.AbstractPolynomialLike, variables)

Compute the minimum and maximum (total) degree of `f` with respect to the given variables.
"""
function minmaxdegree(f::MP.AbstractPolynomialLike, variables)
    d_min, d_max = typemax(Int), 0
    for t in f
        d = sum(MP.degree(t, v) for v in variables)
        d_min = min(d, d_min)
        d_max = max(d, d_max)
    end
    d_min, d_max
end

"""
    uniquevar(f::MP.AbstractPolynomialLike, tag=:x0)
    uniquevar(F::Vector{<:MP.AbstractPolynomialLike}, tag=:x0)

Creates a unique variable.
"""
uniquevar(f::MP.AbstractPolynomialLike, tag=:x0) = MP.similarvariable(f, gensym(tag))
uniquevar(F::Vector{<:MP.AbstractPolynomialLike}, tag=:x0) = uniquevar(F[1], tag)

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
function homogenize(F::Vector{<:MP.AbstractPolynomialLike}, var=uniquevar(F); parameters=nothing)
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
function homogenize(f::MP.AbstractPolynomialLike, variables::Vector{T}, var=uniquevar(f)) where {T<:MP.AbstractVariable}
    _, d_max = minmaxdegree(f, variables)
    MP.polynomial(map(f) do t
        d = sum(MP.degree(t, v) for v in variables)
        var^(d_max - d)*t
    end)
end
function homogenize(F::Vector{<:MP.AbstractPolynomialLike}, variables::Vector{T}, var=uniquevar(F)) where {T<:MP.AbstractVariable}
    map(f -> homogenize(f, variables, var), F)
end
