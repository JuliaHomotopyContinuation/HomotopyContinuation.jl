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

import IterTools: imap

abstract type AbstractResult end
abstract type AbstractSolver end

######################
## MultiplicityInfo ##
######################

"""
    MultiplicityInfo

This contains information about the multiplicities of the solutions.
"""
struct MultiplicityInfo
    # indices of multiple solutions grouped by multiplicity
    multiplicities::Dict{Int,Vector{Vector{Int}}}
    # This set stores all path numbers which we want to ignore
    # if we don't show multiple results
    multiple_indicator::Set{Int}
end

function MultiplicityInfo(pathresults::Vector{PathResult})
    multiple_indicator = Set{Int32}()
    multiplicities = compute_multiplicities(pathresults)
    for clusters in values(multiplicities), cluster in clusters
        for i = 2:length(cluster)
            push!(multiple_indicator, path_number(pathresults[cluster[i]]))
        end
    end
    # the multiplicities are currently with respect to the indices of the pathresults
    # but we want to have them with respect to path numbers in case pathresults
    # is just a subset of all paths
    for k in keys(multiplicities)
        multiplicities[k] = map(multiplicities[k]) do cluster
            map(i -> path_number(pathresults[i]), cluster)
        end
    end
    MultiplicityInfo(multiplicities, multiple_indicator)
end

function compute_multiplicities(result::Vector{<:PathResult}; kwargs...)
    D = Dict{Int,Vector{Vector{Int}}}()
    for m in multiplicities(solution, result)
        if haskey(D, length(m))
            push!(D[length(m)], m)
        else
            D[length(m)] = [m]
        end
    end
    D
end

function assign_multiplicities!(path_results::Vector{<:PathResult}, I::MultiplicityInfo)
    #assign multiplicities
    for (k, clusters) in I.multiplicities, cluster in clusters, vᵢ in cluster
        path_results[vᵢ].multiplicity =
            max(k, something(path_results[vᵢ].winding_number, 1))
    end
    path_results
end



"""
    Result

The result of [`solve`](@ref). This is a wrapper around the results of each single path
([`PathResult`](@ref)) and it contains some additional information like a random seed to
replicate the result.
"""
struct Result <: AbstractResult
    path_results::Vector{PathResult}
    tracked_paths::Int
    seed::Union{Nothing,UInt32}
    start_system::Union{Nothing,Symbol}
    multiplicity_info::MultiplicityInfo
end
function Result(
    path_results::Vector{PathResult};
    seed::Union{Nothing,UInt32} = nothing,
    tracked_paths = length(path_results),
    start_system = nothing,
)
    multiplicity_info = MultiplicityInfo(filter(is_singular, path_results))
    assign_multiplicities!(path_results, multiplicity_info)
    Result(path_results, tracked_paths, seed, start_system, multiplicity_info)
end

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

multiple_indicator(::AbstractResult) =
    error("[`multiple_indicator`] Not defined for abstract results yet")
multiple_indicator(R::Result) = R.multiplicity_info.multiple_indicator
is_multiple_result(r::PathResult, R::AbstractResult) =
    path_number(r) ∈ multiple_indicator(R)
is_multiple_result(r::PathResult, R::AbstractVector{PathResult}) = false


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

function ResultStatistics(result::Result; real_tol::Float64 = 1e-6)
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
        elseif is_singular(r)
            if is_real(r, real_tol)
                real_singular += 1
                real_singular_with_multiplicity += something(multiplicity(r), 1)
            end
            singular += 1
            singular_with_multiplicity += something(multiplicity(r), 1)
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
    statistics(R::Result; real_tol = 1e-6)

Statistic about the number of (real) singular and non-singular solutions etc.
"""
statistics(r; kwargs...) = ResultStatistics(r; kwargs...)

const Results = Union{Result,AbstractVector{<:PathResult}}
const AbstractResults = Union{AbstractResult,AbstractVector{<:PathResult}}

"""
    results(
        result;
        only_real = false,
        real_tol = 1e-6,
        only_nonsingular = false,
        only_singular = false,
        only_finite = true,
        multiple_results = false,
    )
    results(f, result; options...)

Return all [`PathResult`](@ref)s for which satisfy the given conditions and apply,
if provided, the function `f`.
"""
results(R::AbstractResults; kwargs...) = results(identity, R; kwargs...)
function results(
    f::Function,
    R::AbstractResults;
    only_real::Bool = false,
    real_tol::Float64 = 1e-6,
    only_nonsingular::Bool = false,
    only_singular::Bool = false,
    only_finite::Bool = true,
    multiple_results::Bool = false,
)
    if multiple_results == false && !(typeof(R)<:Results)
        println("Warning: Since result is a ResultIterator, counting multiple results")
        multiple_results = true
    end
    filter_function =
        r ->
            (!only_real || is_real(r, real_tol)) &&
            (!only_nonsingular || is_nonsingular(r)) &&
            (!only_singular || is_singular(r)) &&
            (!only_finite || is_finite(r)) &&
            (multiple_results || !is_multiple_result(r, R))
    return_iter = IterTools.imap(f, Iterators.filter(filter_function, R))
    if typeof(R) <: Results
        return (collect(return_iter))
    else
        return (return_iter)
    end
end



"""
    ResultIterator

"""
struct ResultIterator{Iter} <: AbstractResult
    starts::Iter                       # The start solution iterator
    S::AbstractSolver
    bitmask::Union{BitVector,Nothing}  # `nothing` means no filtering
end

function ResultIterator(starts::Iter, S::AbstractSolver; predicate = nothing) where {Iter}
    bitmask = isnothing(predicate) ? nothing : BitVector([predicate(S(x)) for x in starts])
    return (ResultIterator{Iter}(starts, S, bitmask))
end


function Base.show(io::IO, ri::ResultIterator{Iter}) where {Iter}
    print(io, "ResultIterator induced by $Iter")
    !isnothing(ri.bitmask) && print(" and a filtering bitmask")
end

function Base.IteratorSize(ri::ResultIterator)
    ri.bitmask === nothing ? Base.IteratorSize(ri.starts) : Base.HasLength()
end
Base.IteratorEltype(::ResultIterator) = Base.HasEltype()
Base.eltype(::ResultIterator) = PathResult



function Base.iterate(ri::ResultIterator, state = nothing) #States of the induced iterator are pairs (i::Int,state)
    native_state = state === nothing ? 0 : state[1]
    next_ss = state === nothing ? iterate(ri.starts) : iterate(ri.starts, state[2])
    next_ss === nothing && return nothing  # End of iteration

    start_value, new_ss_state = next_ss
    new_state = (native_state+1, new_ss_state)

    if ri.bitmask === nothing
        return (track(ri.S.trackers[1], start_value), new_state)
    else
        # Apply the filter by checking the bitmask
        while !ri.bitmask[new_state[1]] && next_ss !== nothing
            next_ss = iterate(ri.starts, new_state[2])
            next_ss === nothing && return nothing  # End of iteration
            start_value, new_ss_state = next_ss
            new_state = (new_state[1]+1, new_ss_state)
        end

        return (track(ri.S.trackers[1], start_value), new_state)
    end
end


function Base.length(ri::ResultIterator)
    if Base.IteratorSize(ri) == Base.SizeUnknown()
        k = 0
        for _ in ri.starts
            k += 1
        end
        k
    elseif Base.IteratorSize(ri) == Base.HasLength()
        if ri.bitmask !== nothing
            return (sum(ri.bitmask))
        else
            return (length(ri.starts))
        end
    end
end


function bitmask(f::Function, ri::ResultIterator)
    BitVector(map(f, ri))
end

function bitmask_filter(f::Function, ri::ResultIterator)
    bm = bitmask(f, ri)
    return (ResultIterator(ri.starts, ri.S, bm))
end

function trace(iter::ResultIterator)
    s = solution(first(iter))
    mapreduce(x->solution(x), +, iter, init = 0.0 .* s)
end

function Result(ri::ResultIterator)
    C = collect(ri)
    for i = 1:length(C)
        C[i].path_number = i
    end
    Result(C; seed = ri.S.seed, start_system = ri.S.start_system)
end

######################
## Helper functions ##
######################

"""
    nresults(
        result;
        only_real = false,
        real_tol = 1e-6,
        only_nonsingular = false,
        only_singular = false,
        only_finite = true,
        multiple_results = false,
    )

Count the number of results which satisfy the corresponding conditions. See also
[`results`](@ref).
"""
function nresults(
    R::AbstractResults;
    only_real::Bool = false,
    real_tol::Float64 = 1e-6,
    only_nonsingular::Bool = false,
    only_singular::Bool = false,
    onlyfinite::Bool = true, # deprecated
    only_finite::Bool = onlyfinite,
    multiple_results::Bool = false,
)
    if multiple_results == false && !(typeof(R)<:Results)
        println("Warning: Since result is a ResultIterator, counting multiple results")
        multiple_results = true
    end
    count(R) do r
        (!only_real || is_real(r, real_tol)) &&
            (!only_nonsingular || is_nonsingular(r)) &&
            (!only_singular || is_singular(r)) &&
            (!only_finite || isfinite(r)) &&
            (multiple_results || !is_multiple_result(r, R))
    end
end


"""
    solutions(result; only_nonsingular = true, conditions...)

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
solutions(result::AbstractResults; only_nonsingular = true, kwargs...) =
    results(solution, result; only_nonsingular = only_nonsingular, kwargs...)

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
function real_solutions(result::AbstractResults; tol::Float64 = 1e-6, kwargs...)
    results(real ∘ solution, result; only_real = true, real_tol = tol, kwargs...)
end


"""
    nonsingular(result; conditions...)

Return all [`PathResult`](@ref)s for which the solution is non-singular.
This is just a shorthand for `results(R; only_nonsingular=true, conditions...)`.
For the possible `conditions` see [`results`](@ref).
"""
nonsingular(R::AbstractResults; kwargs...) = results(R; only_nonsingular = true, kwargs...)

"""
    singular(result; multiple_results=false, kwargs...)

Return all [`PathResult`]s for which the solution is singular.
If `multiple_results=false` only one point from each cluster of multiple solutions is
returned. If `multiple_results = true` all singular solutions in `R` are returned.
For the possible `kwargs` see [`results`](@ref).
"""
function singular(R::AbstractResults; kwargs...)
    results(R; only_singular = true, kwargs...)
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

"""
    nsolutions(result; only_nonsingular = true, options...)

The number of solutions. See [`results`](@ref) for the possible options.
"""
nsolutions(R::AbstractResults; only_nonsingular = true, options...) =
    nresults(R; only_nonsingular = only_nonsingular, options...)

"""
    nsingular(
        result;
        counting_multiplicities = false,
        kwargs...,
    )

The number of singular solutions. A solution is considered singular if its winding number is
larger than 1 or the condition number is larger than `tol`.
If `counting_multiplicities=true` the number of singular solutions times their
multiplicities is returned.
"""
function nsingular(R::AbstractResult; counting_multiplicities::Bool = false, kwargs...)
    if counting_multiplicities
        return (nresults(R; only_singular = true, multiple_results = true, kwargs...))
    else
        return (nresults(R; only_singular = true, multiple_results = false, kwargs...))
    end
end


"""
    nat_infinity(result)

The number of solutions at infinity.
"""
nat_infinity(R::AbstractResults) = count(is_at_infinity, R)

"""
    nexcess_solutions(result)

The number of exess solutions. See also [`excess_solution_check`](@ref).
"""
nexcess_solutions(R::AbstractResults) = count(is_excess_solution, R)

"""
    nfailed(result)

The number of failed paths.
"""
nfailed(R::AbstractResults) = count(is_failed, R)

"""
    nnonsingular(result)

The number of non-singular solutions. See also [`is_singular`](@ref).
"""
nnonsingular(R::AbstractResults) = count(is_nonsingular, R)

"""
    nreal(result; tol=1e-6)

The number of real solutions. See also [`is_real`](@ref).
"""
nreal(R::Result; tol = 1e-6) = nresults(R, only_real = true, real_tol = tol)
nreal(R::AbstractResult; tol = 1e-6) =
    nresults(R, only_real = true, real_tol = tol, multiple_results = true)

"""
    ntracked(result)

Returns the total number of paths tracked.
"""
ntracked(R::AbstractResult) = length(R) #error("ntracked not implemented for abstract results")
ntracked(R::Result) = R.tracked_paths

###
### Show
####
function Base.show(io::IO, x::Result)
    s = statistics(x)
    total = s.nonsingular + s.singular
    header = "Result with $total $(plural("solution", total))"
    println(io, header)
    println(io, "="^(length(header)))
    println(io, "• $(ntracked(x)) $(plural("path", ntracked(x))) tracked")
    println(
        io,
        "• $(s.nonsingular) non-singular $(plural("solution", s.nonsingular)) ",
        "($(s.real_nonsingular) real)",
    )
    if s.singular > 0
        println(
            io,
            "• $(s.singular) singular $(plural("solution", s.singular)) ",
            "($(s.real_singular) real)",
        )
    end
    # s.at_infinity > 0 &&
    #     println(io, "• $(s.at_infinity) $(plural("solution", s.at_infinity)) at infinity")
    s.excess_solution > 0 && println(
        io,
        "• $(s.excess_solution) excess $(plural("solution", s.excess_solution))",
    )
    println(io, "• random_seed: ", sprint(show, seed(x)))
    if !isnothing(x.start_system)
        println(io, "• start_system: :", x.start_system)
    end
    if s.singular > 0
        println(io, "• multiplicity table of singular solutions:")
        singular_multiplicities_table(io, x, s)
    end
end

function TreeViews.treelabel(io::IO, x::Result, ::MIME"application/prs.juno.inline")
    s = statistics(x)
    total = s.nonsingular + s.singular
    print(
        io,
        "<span><span class=\"syntax--support syntax--type syntax--julia\">Result</span>" *
        " with $total $(plural("solution", total))<span>",
    )
end
TreeViews.hastreeview(::Result) = true
TreeViews.numberofnodes(::Result) = 10
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
    elseif i == 5 && s.excess_solution > 0
        print(io, "$(s.excess_solution) excess $(plural("solution", s.excess_solution))")
    elseif i == 6 && s.at_infinity > 0
        print(io, "$(s.at_infinity) at_infinity")
    elseif i == 7 && s.failed > 0
        print(io, "$(s.failed) failed")
    elseif i == 8
        print(io, "Random seed used")
    elseif i == 9 && !isnothing(x.start_system)
        print(io, "start_system")
    elseif i == 10 && s.singular > 0
        print(io, "  multiplicity table of singular solutions: \n")
        singular_multiplicities_table(io, x, s)
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
    elseif i == 5 && s.excess_solution > 0
        return filter(is_excess_solution, path_results(r))
    elseif i == 6 && s.at_infinity > 0
        return at_infinity(r)
    elseif i == 7 && s.failed > 0
        return failed(r)
    elseif i == 8
        return seed(r)
    elseif i == 9 && !isnothing(r.start_system)
        return r.start_system
    elseif i == 10 && s.singular > 0
        return missing
    end
    missing
end


function singular_multiplicities_table(io, result::Result, stats = statistics(result))
    M = result.multiplicity_info.multiplicities
    if isempty(M)
        n_higher_mults_total = 0
    else
        n_higher_mults_total = sum(length, values(M))
    end

    headers = ["mult.", "total", "# real", "# non-real"]
    mult_1_exists = n_higher_mults_total < stats.singular
    data = Matrix{Int}(undef, mult_1_exists + length(M), 4)

    n_higher_real = 0
    i = mult_1_exists + 1
    for k in sort(collect(keys(M)))
        data[i, 1] = k
        n_real_solsᵢ = count(Mᵢ -> is_real(result.path_results[Mᵢ[1]]), M[k])
        n_solsᵢ = length(M[k])
        data[i, 2] = n_solsᵢ
        data[i, 3] = n_real_solsᵢ
        data[i, 4] = n_solsᵢ - n_real_solsᵢ
        n_higher_real += n_real_solsᵢ
        i += 1
    end

    if mult_1_exists
        data[1, 1] = 1
        data[1, 2] = stats.singular - n_higher_mults_total
        data[1, 3] = stats.real_singular - n_higher_real
        data[1, 4] = data[1, 2] - data[1, 3]
    end

    PrettyTables.pretty_table(
        io,
        data;
        header = headers,
        tf = PrettyTables.tf_unicode_rounded,
        alignment = :c,
        header_crayon = PrettyTables.Crayon(bold = false),
        border_crayon = PrettyTables.Crayon(faint = true),
    )
end
