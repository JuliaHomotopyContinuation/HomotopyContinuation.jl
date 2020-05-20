export Result,
    ResultStatistics,
    statistics,
    seed,
    path_results,
    results,
    nresults,
    solutions,
    real_solutions,
    nonsingular,
    singular,
    at_infinity,
    failed,
    nfinite,
    nsolutions,
    nsingular,
    nat_infinity,
    nexcess_solutions,
    nfailed,
    nnonsingular,
    nreal,
    ntracked

"""
    Result

The result of [`solve`](@ref). This is a wrapper around the results of each single path
([`PathResult`](@ref)) and it contains some additional information like a random seed to
replicate the result.
"""
struct Result
    path_results::Vector{PathResult}
    tracked_paths::Int
    seed::Union{Nothing,UInt32}
end
Result(
    path_results;
    seed::Union{Nothing,UInt32} = nothing,
    tracked_paths = length(path_results),
) = Result(path_results, tracked_paths, seed)

Base.size(r::Result) = (length(r),)
Base.length(r::Result) = length(r.path_results)
Base.getindex(r::Result, I) = getindex(r.path_results, I)

Base.iterate(r::Result) = iterate(r.path_results)
Base.iterate(r::Result, state) = iterate(r.path_results, state)
Base.lastindex(r::Result) = lastindex(r.path_results)
Base.eltype(r::Type{Result}) = PathResult
Base.keys(r::Result) = 1:length(r.path_results)
Base.values(r::Result) = r.path_results

"""
    seed(::Result)

Returns the seed to replicate the result.
"""
seed(R::Result) = R.seed

"""
    path_results(::Result)

Returns the stored [`PathResult`](@ref)s.
"""
path_results(R::Result) = R.path_results


#TODO: is_multiple_result(r::PathResult, r::Result)
is_multiple_result(path_result::PathResult, result) = false


Base.@kwdef struct ResultStatistics
    total::Int
    nonsingular::Int
    singular::Int
    singular_with_multiplicity::Int = singular
    real::Int
    real_nonsingular::Int
    real_singular::Int
    real_singular_with_multiplicity::Int = real_singular
    at_infinity::Int
    excess_solution::Int = 0
    failed::Int = 0
end
Base.show(io::IO, stats::ResultStatistics) = print_fieldnames(io, stats)

function ResultStatistics(
    result::Result;
    real_tol::Float64 = 1e-6,
    singular_tol::Float64 = 1e10,
)
    failed = at_infinity = excess_solution = 0
    nonsingular = singular = real_nonsingular = real_singular = 0
    singular_with_multiplicity = real_singular_with_multiplicity = 0
    for r in result
        is_multiple_result(r, result) && continue
        if is_failed(r)
            failed += 1
        elseif is_at_infinity(r)
            at_infinity += 1
        elseif is_excess_solution(r)
            excess_solution += 1
        elseif is_singular(r, singular_tol)
            if is_real(r, real_tol)
                real_singular += 1
                real_singular_with_multiplicity += unpack(multiplicity(r), 1)
            end
            singular += 1
            singular_with_multiplicity += unpack(multiplicity(r), 1)
        else # finite, nonsingular
            if is_real(r, real_tol)
                real_nonsingular += 1
            end
            nonsingular += 1
        end
    end
    ResultStatistics(
        nonsingular = nonsingular,
        singular = singular,
        singular_with_multiplicity = singular_with_multiplicity,
        real_nonsingular = real_nonsingular,
        real_singular = real_singular,
        real_singular_with_multiplicity = real_singular_with_multiplicity,
        real = real_nonsingular + real_singular,
        at_infinity = at_infinity,
        excess_solution = excess_solution,
        failed = failed,
        total = result.tracked_paths,
    )
end

"""
    statistics(R::Result; real_tol = 1e-6, singular_tol = 1e10)

Statistic about the number of (real) singular and non-singular solutions etc.
"""
statistics(r; kwargs...) = ResultStatistics(r; kwargs...)

const Results = Union{Result,AbstractVector{<:PathResult}}

"""
    results(
        result;
        only_real = false,
        real_tol = 1e-6,
        only_nonsingular = false,
        only_singular = false,
        singular_tol = 1e10,
        only_finite = true,
        multiple_results = false,
    )
    results(f, result; options...)

Return all [`PathResult`](@ref)s for which satisfy the given conditions and apply if provided
the function `f`.
"""
results(R::Results; kwargs...) = results(identity, R; kwargs...)
function results(
    f::Function,
    R::Results;
    only_real::Bool = false,
    real_tol::Float64 = 1e-6,
    only_nonsingular::Bool = false,
    only_singular::Bool = false,
    singular_tol::Float64 = 1e10,
    onlyfinite::Bool = true, # deprecated
    only_finite::Bool = onlyfinite,
    multiple_results::Bool = false,
)
    [
        f(r)
        for
        r in R if
        (!only_real || is_real(r, real_tol)) &&
            (!only_nonsingular || is_nonsingular(r, singular_tol)) &&
            (!only_singular || is_singular(r, singular_tol)) &&
            (!only_finite || is_finite(r)) &&
            (multiple_results || !is_multiple_result(r, R))
    ]
end

"""
    nresults(
        result;
        only_real = false,
        real_tol = 1e-6,
        only_nonsingular = false,
        only_singular = false,
        singular_tol = 1e10,
        only_finite = true,
        multiple_results = false,
    )

Count the number of solutions which satisfy the corresponding conditions. See also
[`results`](@ref).
"""
function nresults(
    R::Results;
    only_real::Bool = false,
    real_tol::Float64 = 1e-6,
    only_nonsingular::Bool = false,
    only_singular::Bool = false,
    singular_tol::Float64 = 1e10,
    onlyfinit::Bool = true, # deprecated
    only_finite::Bool = onlyfinit,
    multiple_results::Bool = false,
)
    count(R) do r
        (!only_real || is_real(r, real_tol)) &&
            (!only_nonsingular || is_nonsingular(r, singular_tol)) &&
            (!only_singular || is_singular(r, singular_tol)) &&
            (!only_finite || isfinite(r)) &&
            (multiple_results || !is_multiple_result(r, R))
    end
end

"""
    solutions(result; conditions...)

Returns all solutions for which the given conditions apply, see [`results`](@ref) for the
possible conditions.

## Example
```julia-repl
julia> @var x y
julia> F = System([(x-2)y, y+x+3]);
julia> solutions(solve(F))
2-element Array{Array{Complex{Float64},1},1}:
 [2.0 + 0.0im, -5.0 + 0.0im]
 [-3.0 + 0.0im, 0.0 + 0.0im]
```
"""
solutions(result::Results; kwargs...) = results(solution, result; kwargs...)

"""
    real_solutions(result; tol=1e-6, conditions...)

Return all real solution for which the given conditions apply.
For the possible `conditions` see [`results`](@ref).
Note that `only_real` is always `true` and `real_tol` is now `tol`.

## Example
```julia-repl
julia> @var x y;
julia> F = System([(x-2)y, y+x+3]);
julia> real_solutions(solve(F))
2-element Array{Array{Float64,1},1}:
 [2.0, -5.0]
 [-3.0, 0.0]
```
"""
function real_solutions(result::Results; tol::Float64 = 1e-6, kwargs...)
    results(real ∘ solution, result; only_real = true, real_tol = tol, kwargs...)
end


"""
    nonsingular(result; conditions...)

Return all [`PathResult`](@ref)s for which the solution is non-singular.
This is just a shorthand for `results(R; only_nonsingular=true, conditions...)`.
For the possible `conditions` see [`results`](@ref).
"""
nonsingular(R::Results; kwargs...) = results(R; only_nonsingular = true, kwargs...)

"""
    singular(result; tol=1e10, multiple_results=false, kwargs...)

Return all [`PathResult`]s for which the solution is singular.
A solution is labeled singular
if the condition number is greater than `singular_tol`, or if the winding number is > 1.
If `multiple_results=false` only one point from each cluster of multiple solutions is returned.
If If `multiple_results=true` all singular solutions in `R` are returned.
For the possible `kwargs` see [`results`](@ref).
"""
function singular(R::Results; tol::Float64 = 1e10, kwargs...)
    results(R; only_singular = true, singular_tol = tol, kwargs...)
end

"""
    real(result, tol=1e-6)

Get all results where the solutions are real with the given tolerance `tol`.
See [`is_real`](@ref) for details regarding the determination of 'realness'.
"""
Base.real(R::Results; tol::Float64 = 1e-6) = filter(r -> is_real(r, tol), path_results(R))

"""
    failed(result)
Get all results where the path tracking failed.
"""
failed(R::Results) = filter(is_failed, path_results(R))

"""
    at_infinity(result)

Get all results where the solutions is at infinity.
"""
at_infinity(R::Results) = filter(is_at_infinity, path_results(R))

# """
#     multiplicities(result; tol =1 e-6)
#
# Returns a `Vector` of `Vector{PathResult}`s grouping the `PathResult`s whose solutions
# appear with multiplicities *greater* than 1 in 'result'.
# Two solutions are regarded as equal, when their pairwise distance is less than 'tol'.
# """
# function multiplicities(results::Results; tol::Float64 = 1e-6)
#     map(i -> results[i], multiplicities(solution, results; tol = tol))
# end

"""
    nsolutions(result)

The number of solutions.
"""
nsolutions(R::Results) = nresults(R)

"""
    nsingular(
        result;
        singular_tol = 1e10,
        counting_multiplicities = false,
        kwargs...,
    )

The number of singular solutions. A solution is considered singular if its winding number is
larger than 1 or the condition number is larger than `tol`.
If `counting_multiplicities=true` the number of singular solutions times their
multiplicities is returned.
"""
function nsingular(
    R::Results;
    singular_tol::Float64 = 1e10,
    counting_multiplicities::Bool = false,
    kwargs...,
)
    S = results(
        R;
        only_singular = true,
        multiple_results = false,
        singular_tol = singular_tol,
        kwargs...,
    )
    isempty(S) && return 0
    counting_multiplicities && return sum(multiplicity, S)
    length(S)
end


"""
    nat_infinity(result)

The number of solutions at infinity.
"""
nat_infinity(R::Results) = count(is_at_infinity, R)

"""
    nexcess_solutions(result)

The number of exess solutions. See also [`excess_solution_check`](@ref).
"""
nexcess_solutions(R::Results) = count(is_excess_solution, R)

"""
    nfailed(result)

The number of failed paths.
"""
nfailed(R::Results) = count(is_failed, R)

"""
    nnonsingular(result; tol=1e-10)

The number of non-singular solutions. See also [`is_singular`](@ref).
"""
nnonsingular(R::Result; tol = 1e10) = count(r -> is_nonsingular(r, tol), R)

"""
    nreal(result; tol=1e-6)

The number of real solutions. See also [`is_real`](@ref).
"""
nreal(R::Results; tol = 1e-6) = count(r -> is_real(r, tol), R)


"""
    ntracked(result)

Returns the total number of paths tracked.
"""
ntracked(R::Result) = R.tracked_paths

###
### Show
####
plural(singularstr, n) = n == 1 ? singularstr : singularstr * "s"

function Base.show(io::IO, x::Result)
    s = statistics(x)
    header = "Result with $(s.nonsingular + s.singular) solutions"
    println(io, header)
    println(io, "="^(length(header)))
    println(
        io,
        "• $(s.nonsingular) non-singular $(plural("solution", s.nonsingular)) ",
        "($(s.real_nonsingular) real)",
    )
    println(
        io,
        "• $(s.singular) singular $(plural("solution", s.singular)) ",
        "($(s.real_singular) real)",
    )
    s.at_infinity > 0 &&
        println(io, "• $(s.at_infinity) $(plural("solution", s.at_infinity)) at infinity")
    s.excess_solution > 0 && println(
        io,
        "• $(s.excess_solution) excess $(plural("solution", s.excess_solution))",
    )
    println(io, "• $(ntracked(x)) paths tracked")
    print(io, "• random seed: ")
    show(io, seed(x))
    println()
    if s.singular > 0
        println(io, "• multiplicity table of singular solutions:")
        #     singular_multiplicities_table(io, x, s)
    end
end

function TreeViews.treelabel(io::IO, x::Result, ::MIME"application/prs.juno.inline")
    s = statistics(x)
    print(
        io,
        "<span><span class=\"syntax--support syntax--type syntax--julia\">Result</span> with $(s.nonsingular + s.singular) solutions<span>",
    )
end
TreeViews.hastreeview(::Result) = true
TreeViews.numberofnodes(::Result) = 8
function TreeViews.nodelabel(io::IO, x::Result, i::Int, ::MIME"application/prs.juno.inline")
    s = statistics(x)
    if i == 1
        print(io, "Paths tracked")
    elseif i == 2 && s.nonsingular > 0
        print(io, "$(s.nonsingular) non-singular ($(s.real_nonsingular) real)")
    elseif i == 3 && s.singular > 0
        print(io, "$(s.singular) singular ($(s.real_singular) real)")
    elseif i == 4 && (s.real_nonsingular + s.real_singular) > 0
        print(io, "$(s.real_nonsingular+s.real_singular) real")
    elseif i == 5 && s.at_infinity > 0
        print(io, "$(s.at_infinity) at_infinity")
    elseif i == 6 && s.failed > 0
        print(io, "$(s.failed) failed")
    elseif i == 7
        print(io, "Random seed used")
    elseif i == 8 && s.singular > 0
        print(io, "  multiplicity table of singular solutions: \n")
        # singular_multiplicities_table(io, x, s)
    end
end

function TreeViews.treenode(r::Result, i::Integer)
    s = statistics(r)
    if i == 1
        return ntracked(r)
    elseif i == 2 && s.nonsingular > 0
        return results(r, only_nonsingular = true)
    elseif i == 3 && s.singular > 0
        return results(r, only_singular = true)
    elseif i == 4 && (s.real_nonsingular + s.real_singular) > 0
        return results(r, only_real = true)
    elseif i == 5 && s.at_infinity > 0
        return at_infinity(r)
    elseif i == 6 && s.failed > 0
        return failed(r)
    elseif i == 7
        return seed(r)
    elseif i == 8 && s.singular > 0
        return missing
    end
    missing
end




# function singular_multiplicities_table(io, result::Result, stats = statistics(result))
#     M = result.multiplicity_info[].multiplicities
#     if isempty(M)
#         n_higher_mults_total = 0
#     else
#         n_higher_mults_total = sum(length, values(M))
#     end
#
#     headers = ["mult.", "total", "# real", "# non-real"]
#     mult_1_exists = n_higher_mults_total < stats.singular
#     data = Matrix{Int}(undef, mult_1_exists + length(M), 4)
#
#     n_higher_real = 0
#     i = mult_1_exists + 1
#     for k in sort(collect(keys(M)))
#         data[i, 1] = k
#         n_real_solsᵢ = count(Mᵢ -> is_real(result.pathresults[Mᵢ[1]]), M[k])
#         n_solsᵢ = length(M[k])
#         data[i, 2] = n_solsᵢ
#         data[i, 3] = n_real_solsᵢ
#         data[i, 4] = n_solsᵢ - n_real_solsᵢ
#         n_higher_real += n_real_solsᵢ
#         i += 1
#     end
#
#     if mult_1_exists
#         data[1, 1] = 1
#         data[1, 2] = stats.singular - n_higher_mults_total
#         data[1, 3] = stats.real_singular - n_higher_real
#         data[1, 4] = data[1, 2] - data[1, 3]
#     end
#
#     PrettyTables.pretty_table(
#         io,
#         data,
#         headers;
#         tf = PrettyTables.unicode,
#         alignment = :c,
#         header_crayon = PrettyTables.Crayon(bold = false),
#         border_crayon = PrettyTables.Crayon(faint = true),
#     )
# end
