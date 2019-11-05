export Result,
       nresults,
       nsolutions,
       nfinite,
       nsingular,
       nat_infinity,
       nfailed,
       nnonsingular,
       nreal,
       ntracked,
       finite,
       results,
       mapresults,
       failed,
       at_infinity,
       singular,
       nonsingular,
       seed,
       solutions,
       real_solutions,
       multiplicities,
       statistics,
       multiplicities!

######################
## MultiplicityInfo ##
######################

"""
    MultiplicityInfo

This contains informations about the multiplicities of the solutions.
"""
struct MultiplicityInfo
    multiplicities::Dict{Int,Vector{Vector{Int}}}
    # This set stores all path numbers which we want to ignore
    # if we don't show multiple results
    multiple_indicator::Set{Int32}
end

function MultiplicityInfo(pathresults::Vector{<:PathResult}; tol = 1e-6)
    multiple_indicator = Set{Int32}()
    multiplicities = compute_multiplicities(pathresults; tol = float(tol))
    for clusters in values(multiplicities), cluster in clusters
        for i = 2:length(cluster)
            push!(multiple_indicator, path_number(pathresults[cluster[i]]))
        end
    end
    MultiplicityInfo(multiplicities, multiple_indicator)
end

function compute_multiplicities(result::Vector{<:PathResult}; kwargs...) where {T}
    D = Dict{Int,Vector{Vector{Int}}}()
    compute_multiplicities!(D, result; kwargs...)
end
function compute_multiplicities!(
    D::Dict,
    result::Vector{<:PathResult};
    tol::Float64 = 1e-6,
) where {T}
    for m in multiplicities(solution, result, tol = tol)
        if haskey(D, length(m))
            push!(D[length(m)], m)
        else
            D[length(m)] = [m]
        end
    end
    D
end

"""
    assign_multiplicities!(path_results::Vector{<:PathResult}, I::MultiplicityInfo)

Assign the path results their multiplicity according to the multiplicity info `I`.
"""
function assign_multiplicities!(path_results::Vector{<:PathResult}, I::MultiplicityInfo)
    #assign multiplicities
    for r in path_results
        r.multiplicity[] = 1
    end
    for (k, clusters) in I.multiplicities, cluster in clusters, vᵢ in cluster
        path_results[vᵢ].multiplicity[] = k
    end
    path_results
end

is_multiple_result(r::PathResult, I::MultiplicityInfo) =
    path_number(r) ∈ I.multiple_indicator


"""
    Result{V<:AbstractVector}

The result of `solve`. This is a wrapper around the results of each single path
([`PathResult`](@ref)) and it contains some additional informations like a random seed to
replicate the result.
"""
struct Result{V}
    pathresults::Vector{PathResult{V}}
    tracked_paths::Int
    seed::Int
    retracked_paths::Int
    multiplicity_info::Base.RefValue{MultiplicityInfo}
end

function Result(
    pathresults::Vector{PathResult{V}},
    tracked_paths::Int,
    seed::Int;
    retracked_paths::Int = 0,
    multiplicity_tol::Float64 = 1e-6,
) where {V}
    multiplicity_info = MultiplicityInfo(pathresults; tol = multiplicity_tol)
    assign_multiplicities!(pathresults, multiplicity_info)
    Result(pathresults, tracked_paths, seed, retracked_paths, Ref(multiplicity_info))
end

Base.length(r::Result) = length(r.pathresults)
Base.getindex(r::Result, I) = getindex(r.pathresults, I)

Base.iterate(r::Result) = iterate(r.pathresults)
Base.iterate(r::Result, state) = iterate(r.pathresults, state)
Base.lastindex(r::Result) = lastindex(r.pathresults)
Base.eltype(r::Type{Result{V}}) where {V} = PathResult{V}
solution_type(::Result{V}) where {V} = V

is_multiple_result(r::PathResult, R::Result) = is_multiple_result(r, R.multiplicity_info[])
is_multiple_result(r::PathResult, R::Vector{<:PathResult}) = false



"""
    multiplicities!(result::Result; tol=1e-6)

Compute the multiplicities of the solutions in `result` with respect to the given tolerance.
"""
function multiplicities!(result::Result; tol = 1e-6)
    multiplicity_info = MultiplicityInfo(result.pathresults; tol = tol)
    result.multiplicity_info[] = multiplicity_info
    assign_multiplicities!(result.pathresults, multiplicity_info)
    result
end

const Results = Union{Result,Vector{<:PathResult}}
const ProjectiveResult = Result{<:PVector}

"""
    nresults(
        result;
        only_real = false,
        real_tol = 1e-6,
        only_nonsingular = false,
        singular_tol = 1e10,
        onlyfinite = true,
    )

The number of solutions which satisfy the corresponding predicates.

## Example
```julia
result = solve(F)
# Get all non-singular results where all imaginary parts are smaller than 1e-8
nresults(result, only_real=true, real_tol=1e-8, only_nonsingular=true)
```
"""
function nresults(
    R::Results;
    only_real = false,
    real_tol = 1e-6,
    only_nonsingular = false,
    only_singular = false,
    singular_tol = 1e10,
    onlyfinite = true,
    multiple_results = false,
)
    count(R) do r
        (!only_real || is_real(r, real_tol)) &&
        (!only_nonsingular || is_nonsingular(r, singular_tol)) &&
        (!only_singular || is_singular(r, singular_tol)) &&
        (!onlyfinite || isfinite(r) || is_projective(r)) &&
        (multiple_results || !is_multiple_result(r, R))
    end
end

"""
    statistics(
        R::Result;
        only_real = false,
        real_tol = 1e-6,
        only_nonsingular = false,
        only_singular = false,
        singular_tol = 1e10,
    )

Statistic about the number of (real) singular and non-singular solutions etc.
Returns a named tuple with the statistics.

## Example
```julia
julia> statistics(solve([x^2+y^2-2, 2x+3y-1]))
    (
     nonsingular = 2,
     singular = 0,
     real_nonsingular = 2,
     real_singular = 0,
     real = 2,
     at_infinity = 0,
     failed = 0,
     total = 2,
    )
"""
function statistics(
    R::Results,
    only_real = false,
    real_tol = 1e-6,
    only_nonsingular = false,
    only_singular = false,
    singular_tol = 1e10,
)

    failed = at_infinity = nonsingular = singular = real_nonsingular = real_singular = 0
    singular_with_multiplicity = real_singular_with_multiplicity = 0
    for r in R
        is_multiple_result(r, R) && continue

        if is_failed(r)
            failed += 1
        elseif is_singular(r, singular_tol)
            if is_real(r, real_tol)
                real_singular += 1
                real_singular_with_multiplicity += unpack(multiplicity(r), 1)
            end
            singular += 1
            singular_with_multiplicity += unpack(multiplicity(r), 1)
        elseif !is_projective(r) && !isfinite(r)
            at_infinity += 1
        else # finite, nonsingular
            if is_real(r, real_tol)
                real_nonsingular += 1
            end
            nonsingular += 1
        end
    end
    (
     nonsingular = nonsingular,
     singular = singular,
     singular_with_multiplicity = singular_with_multiplicity,
     real_nonsingular = real_nonsingular,
     real_singular = real_singular,
     real_singular_with_multiplicity = real_singular_with_multiplicity,
     real = real_nonsingular + real_singular,
     at_infinity = at_infinity,
     failed = failed,
     total = R.tracked_paths,
    )
end

"""
    nfinite(result)

The number of finite solutions.
"""
nfinite(R::Results) = count(isfinite, R)

"""
    nsolutions(result)

The number of solutions.
"""
nsolutions(R::Results) = nresults(R)

"""
    nsingular(
        result;
        singular_tol = 1e10,
        multiplicitytol = 1e-5,
        counting_multiplicities = false,
        kwargs...,
    )

The number of singular solutions. A solution is considered singular if its windingnumber is
larger than 1 or the condition number is larger than `tol`.
If `counting_multiplicities=true` the number of singular solutions times their
multiplicities is returned.
"""
function nsingular(
    R::Results;
    singular_tol = 1e10,
    counting_multiplicities = false,
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
    nafailed(result)

The number of failed paths.
"""
nfailed(R::Results) = count(is_failed, R)

"""
    nnonsingular(result; tol=1e-10)

The number of non-singular solutions.
"""
nnonsingular(R::Result; tol = 1e10) = count(r -> is_nonsingular(r, tol), R)

"""
    nreal(result; tol=1e-6)

The number of real solutions where all imaginary parts of each solution
are smaller than `tol`.
"""
nreal(R::Results; tol = 1e-6) = count(r -> is_real(r, tol), R)


"""
    ntracked(R::Result)

Returns the total number of paths tracked.
"""
ntracked(R::Result) = R.tracked_paths

"""
    seed(result)

The random seed used in the computation.
"""
seed(result::Result) = result.seed



# Filtering
"""
    results(result; only_real=false, real_tol=1e-6, only_nonsingular=false,
                onlysigular=false, singular_tol=1e10, onlyfinite=true, multiple_results=false)

Return all `PathResult`s for which the given conditions apply.

## Example

```julia
R = solve(F)

# This gives us all PathResults considered non-singular and real (but still as a complex vector).
real_solutions = results(R, only_real=true, only_nonsingular=true)
```
"""
results(R::Results; kwargs...) = mapresults(identity, R; kwargs...)
results(f::Function, R::Results; kwargs...) = mapresults(f, R; kwargs...)

"""
    mapresults(f::Function, result; conditions...)

Apply the function `f` to all `PathResult`s for which the given conditions apply. For the
possible conditions see [`results`](@ref).

## Example
```julia
# This gives us all solutions considered real (but still as a complex vector).
real_solutions = mapresults(solution, R, only_real=true)
```
"""
function mapresults(
    f::Function,
    R::Results;
    only_real = false,
    real_tol = 1e-6,
    only_nonsingular = false,
    only_singular = false,
    singular_tol = 1e10,
    onlyfinite = true,
    multiple_results = false,
)
    [f(r)
        for r in R if (!only_real || is_real(r, real_tol)) &&
                      (!only_nonsingular || is_nonsingular(r, singular_tol)) &&
                      (!only_singular || is_singular(r, singular_tol)) &&
                      (!onlyfinite || isfinite(r) || is_projective(r)) &&
                      (multiple_results || !is_multiple_result(r, R))]
end

"""
    solutions(result; conditions...)

Return all solution (as `Vector`s) for which the given conditions apply.
For the possible `conditions` see [`results`](@ref).

## Example
```julia
julia> @polyvar x y
julia> result = solve([(x-2)y, y+x+3]);
julia> solutions(result)
[[2.0+0.0im, -5.0+0.0im], [-3.0+0.0im, 0.0+0.0im]]
```
"""
function solutions(result::Results; kwargs...)
    mapresults(solution, result; kwargs...)
end

"""
    real_solutions(result; tol=1e-6, conditions...)

Return all real solution (as `Vector`s of reals) for which the given conditions apply.
For the possible `conditions` see [`results`](@ref). Note that `only_real` is always `true`
and `real_tol` is now `tol`.

## Example
```julia
julia> @polyvar x y
julia> result = solve([(x-2)y, y+x+3]);
julia> real_solutions(result)
[[2.0, -5.0], [-3.0, 0.0]]
```
"""
function real_solutions(result::Results; only_real = true, tol = 1e-6, kwargs...)
    mapresults(real_vector ∘ solution, result; only_real = true, real_tol = tol, kwargs...)
end
@deprecate realsolutions(result; kwargs...) real_solutions(result; kwargs...)

"""
    nonsingular(result::Results; conditions...)

Return all `PathResult`s for which the solution is non-singular. This is just a shorthand
for `results(R; only_nonsingular=true, conditions...)`. For the possible `conditions` see
[`results`](@ref).
"""
nonsingular(R::Results; kwargs...) = results(R; only_nonsingular = true, kwargs...)

"""
    singular(R::Results; tol=1e10, multiple_results=false, kwargs...)

Return all `PathResult`s for which the solution is singular. A solution is labeled singular
if the condition number is greater than `singular_tol`, or if the winding number is > 1.
If `multiple_results=false` only one point from each cluster of multiple solutions is returned.
If If `multiple_results=true` all singular solutions in `R` are returned.
For the possible `kwargs` see [`results`](@ref).
"""
function singular(R::Results; tol = 1e10, kwargs...)
    results(R; only_singular = true, singular_tol = tol, kwargs...)
end


"""
    finite(result::AffineResults; conditions...)

Return all `PathResult`s for which the solution is finite. This is just a shorthand
for `results(R; onlyfinite=true, conditions...)`. For the possible `conditions` see
[`results`](@ref).
"""
finite(R::Results; kwargs...) = results(R; onlyfinite = true, kwargs...)

"""
    real(result, tol=1e-6)

Get all results where the solutions are real with the given tolerance `tol`.
See [`is_real`](@ref) for details regarding the determination of 'realness'.
"""
Base.real(R::Results; tol = 1e-6) = [r for r in R if is_real(r, tol)]

"""
    failed(result)

Get all results where the path tracking failed.
"""
failed(R::Results) = [r for r in R if is_failed(r)]

"""
    at_infinity(result::AffineResult)

Get all results where the solutions is at infinity.
"""
at_infinity(R::Results) = [r for r in R if is_at_infinity(r)]

"""
    multiplicities(V::Results; tol=1e-6)

Returns a `Vector` of `Vector{PathResult}`s grouping the `PathResult`s whose solutions
appear with multiplicities *greater* 1 in 'V'.
Two solutions are regarded as equal, when their pairwise distance is less than 'tol'.
"""
function multiplicities(results::Results; tol = 1e-6)
    map(i -> results[i], multiplicities(solution, results; tol = tol))
end

####Show function

function Base.show(io::IO, x::Result)
    s = statistics(x)
    header = "Result{$(solution_type(x))} with $(s.nonsingular + s.singular) solutions"
    println(io, header)
    println(io, "="^(length(header)))
    println(
        io,
        "• $(s.nonsingular) non-singular $(plural("solution", s.nonsingular)) ($(s.real_nonsingular) real)",
    )
    println(
        io,
        "• $(s.singular) singular $(plural("solution", s.singular)) ($(s.real_singular) real)",
    )
    s.at_infinity > 0 && println(
        io,
        "• $(s.at_infinity) $(plural("solution", s.at_infinity)) at infinity",
    )
    s.failed > 0 && println(io, "• $(s.failed) failed $(plural("path", s.failed))")
    println(io, "• $(ntracked(x)) paths tracked")
    println(io, "• random seed: $(seed(x))")
    if s.singular > 0
        println(io, "• multiplicity table of singular solutions:")
        singular_multiplicities_table(io, x, s)
    end
end

function Base.show(io::IO, x::ProjectiveResult)
    s = statistics(x)
    header = "Result{$(solution_type(x))} with $(s.nonsingular + s.singular) solutions"
    println(io, header)
    println(io, "="^(length(header)))
    println(
        io,
        "• $(s.nonsingular) non-singular $(plural("solution", s.nonsingular)) ($(s.real_nonsingular) real)",
    )
    println(
        io,
        "• $(s.singular) singular $(plural("solution", s.singular)) ($(s.real_singular) real)",
    )
    s.failed > 0 && println(io, "• $(s.failed) failed $(plural("path", s.failed))")
    println(io, "• $(ntracked(x)) paths tracked")
    println(io, "• random seed: $(seed(x))")
    if s.singular > 0
        println(io, "• multiplicity table of singular solutions:")
        singular_multiplicities_table(io, x, s)
    end
end

TreeViews.hastreeview(::Result) = true
TreeViews.hastreeview(::ProjectiveResult) = true
TreeViews.numberofnodes(::Result) = 8
TreeViews.numberofnodes(::ProjectiveResult) = 7
TreeViews.treelabel(io::IO, x::Result, ::MIME"application/prs.juno.inline") = print(
    io,
    "<span class=\"syntax--support syntax--type syntax--julia\">" *
    "Result{$(solution_type(x))}" * "</span>",
)

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
        singular_multiplicities_table(io, x, s)
    end
end

function TreeViews.treenode(r::Result, i::Integer)
    s = statistics(r)
    if i == 1
        return ntracked(r)
    elseif i == 2 && s.nonsingular > 0
        return finite(r, only_nonsingular = true)
    elseif i == 3 && s.singular > 0
        return finite(r, only_singular = true)
    elseif i == 4 && (s.real_nonsingular + s.real_singular) > 0
        return finite(r, only_real = true)
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

function TreeViews.nodelabel(
    io::IO,
    x::ProjectiveResult,
    i::Int,
    ::MIME"application/prs.juno.inline",
)
    s = statistics(x)
    if i == 1
        print(io, "Paths tracked")
    elseif i == 2 && s.nonsingular > 0
        print(io, "$(s.nonsingular) non-singular ($(s.real_nonsingular) real)")
    elseif i == 3 && s.singular > 0
        print(io, "$(s.singular) singular ($(s.real_singular) real)")
    elseif i == 4 && (s.real_nonsingular + s.real_singular) > 0
        print(io, "$(s.real_nonsingular + s.real_singular) real solutions")
    elseif i == 5 && s.failed > 0
        print(io, "$(s.failed) failed")
    elseif i == 6
        print(io, "Random seed used")
    elseif i == 7 && s.singular > 0
        print(io, "  multiplicity table of singular solutions: \n")
        singular_multiplicities_table(io, x, s)
    end
end

function TreeViews.treenode(r::ProjectiveResult, i::Integer)
    s = statistics(r)
    if i == 1
        return length(r)
    elseif i == 2 && s.nonsingular > 0
        return finite(r, only_nonsingular = true)
    elseif i == 3 && s.singular > 0
        return finite(r, only_singular = true)
    elseif i == 4 && (s.real_nonsingular + s.real_singular) > 0
        return finite(r, only_real = true)
    elseif i == 5 && s.failed > 0
        return failed(r)
    elseif i == 6
        return seed(r)
    elseif i == 7 && s.singular > 0
        return missing
    end
    missing
end

plural(singularstr, n) = n == 1 ? singularstr : singularstr * "s"

function singular_multiplicities_table(io, result::Result, stats = statistics(result))
    M = result.multiplicity_info[].multiplicities
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
        n_real_solsᵢ = count(Mᵢ -> is_real(result.pathresults[Mᵢ[1]]), M[k])
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
        data,
        headers,
        PrettyTables.unicode;
        alignment = :c,
        header_crayon = PrettyTables.Crayon(bold = false),
        border_crayon = PrettyTables.Crayon(faint = true),
    )
end
