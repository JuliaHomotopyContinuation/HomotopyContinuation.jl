export AffineResult, ProjectiveResult,
    nresults, nfinite, nsingular, natinfinity, nfailed, nnonsingular, nreal,
    finite, results, mapresults, failed, atinfinity, singular, nonsingular, seed,
    solutions, realsolutions, multiplicities, uniquesolutions, statistics

"""
    Result

Abstract super type of [`ProjectiveResult`](@ref) and [`AffineResult`](@ref).
"""
abstract type Result end

"""
    ProjectiveResult <: Result

The result of a homogenous system of polynomials.
"""
struct ProjectiveResult{T1, T2, T3} <: Result
    pathresults::Vector{PathResult{T1, T2, T3}}
    seed::Int
end

"""
    AffineResult <: Result

The result of an (non-homogenous) system of polynomials.
"""
struct AffineResult{T1, T2, T3} <: Result
    pathresults::Vector{PathResult{T1, T2, T3}}
    seed::Int
end


Base.length(r::Result) = length(r.pathresults)
Base.getindex(r::Result, I) = getindex(r.pathresults, I)

Base.iterate(r::Result) = iterate(r.pathresults)
Base.iterate(r::Result, state) = iterate(r.pathresults, state)
Base.lastindex(r::Result) = lastindex(r.pathresults)
Base.eltype(r::Type{AffineResult{T1, T2, T3}}) where {T1, T2, T3} = PathResult{T1, T2, T3}
Base.eltype(r::Type{ProjectiveResult{T1, T2, T3}}) where {T1, T2, T3} = PathResult{T1,T2,T3}

const AffineResults = Union{AffineResult, Vector{<:PathResult}}
const ProjectiveResults = Union{ProjectiveResult, Vector{<:PathResult}}
const Results = Union{Result, Vector{<:PathResult}}

"""
    nresults(result; onlyreal=false, realtol=1e-6, onlynonsingular=false, singulartol=1e14, onlyfinite=true)

The number of solutions which satisfy the corresponding predicates.

## Example
```julia
result = solve(F)
# Get all non-singular results where all imaginary parts are smaller than 1e-8
nresults(result, onlyreal=true, realtol=1e-8, onlynonsingular=true)
```
"""
function nresults(R::Result; onlyreal=false, realtol=1e-6,
    onlynonsingular=false, onlysingular=false, singulartol=1e14, onlyfinite=true)
    count(R) do r
        (!onlyreal || isreal(r, realtol)) &&
        (!onlynonsingular || isnonsingular(r, singulartol)) &&
        (!onlysingular || issingular(r, singulartol)) &&
        (!onlyfinite || isfinite(r) || isprojective(r))
    end
end

"""
    statistics(R::Result; onlyreal=false, realtol=1e-6,
        onlynonsingular=false, onlysingular=false, singulartol=1e14)

Statistic about the number of (real) singular and non-singular solutions etc. Returns a named tuple with the statistics.

## Example
```julia
julia> statistics(solve([x^2+y^2-2, 2x+3y-1]))
(nonsingular = 2, singular = 0, real_nonsingular = 2, real_singular = 0, real = 2, atinfinity = 0, failed = 0, total = 2)
"""
statistics(R::Result; kwargs...) = statistics(R.pathresults; kwargs...)
function statistics(R::Vector{<:PathResult}; onlyreal=false, realtol=1e-6,
    onlynonsingular=false, onlysingular=false, singulartol=1e14)

    failed = atinfinity = nonsingular = singular = real_nonsingular = real_singular = 0

    for r in R
        if isfailed(r)
            failed += 1
        elseif issingular(r, singulartol)
            if isreal(r, realtol)
                real_singular += 1
            end
            singular += 1
        elseif !isprojective(r) && !isfinite(r)
            atinfinity += 1
        else # finite, nonsingular
            if isreal(r, realtol)
                real_nonsingular += 1
            end
            nonsingular += 1
        end
    end
    (nonsingular = nonsingular,
    singular = singular,
    real_nonsingular = real_nonsingular,
    real_singular = real_singular,
    real = real_nonsingular + real_singular,
    atinfinity = atinfinity,
    failed = failed,
    total = length(R))
end

"""
    nfinite(affineresult)

The number of finite solutions.
"""
nfinite(R::AffineResult) = count(isfinite, R)

"""
    nsingular(result; tol=1e10)

The number of singular solutions. A solution is considered singular
if its windingnumber is larger than 1 or the condition number is larger than `tol`.
"""
nsingular(R::Results; tol=1e10) = count(R) do r
    issingular(r, tol)
end

"""
    natinfinity(affineresult)

The number of solutions at infinity.
"""
natinfinity(R::AffineResult) = count(isatinfinity, R)

"""
    nafailed(result)

The number of failed paths.
"""
nfailed(R::Result) = count(isfailed, R)

"""
    nnonsingular(result; tol=1e-10)

The number of non-singular solutions.
"""
nnonsingular(R::Result; tol = 1e10) = count(r -> isnonsingular(r, tol), R)

"""
    nreal(result; tol=1e-6)

The number of real solutions where all imaginary parts of each solution
are smaller than `tol`.
"""
nreal(R::Result; tol = 1e-6) = count(r -> isreal(r, tol), R)


"""
    seed(result)

The random seed used in the computation.
"""
seed(result::Result) = result.seed



# Filtering
"""
    results(result; onlyreal=false, realtol=1e-6, onlynonsingular=false, onlysigular=false, singulartol=1e14, onlyfinite=true)

Return all `PathResult`s for which the given conditions apply.

## Example

```julia
R = solve(F)

# This gives us all PathResults considered non-singular and real (but still as a complex vector).
realsolutions = results(R, onlyreal=true, onlynonsingular=true)
```
"""
results(R::Results; kwargs...) = mapresults(identity, R; kwargs...)
# fallback since we introduced mapresults only in 0.2.1
results(f::Function, R::Results; kwargs...) = mapresults(f, R; kwargs...)

"""
    mapresults(f::Function, result; conditions...)

Apply the function `f` to all `PathResult`s for which the given conditions apply. For the possible
conditions see [`results`](@ref).

## Example
```julia
# This gives us all solutions considered real (but still as a complex vector).
realsolutions = results(solution, R, onlyreal=true)
```
"""
function mapresults(f::Function, R::Results;
    onlyreal=false, realtol=1e-6, onlynonsingular=false, onlysingular=false, singulartol=1e14,
    onlyfinite=true)
    [f(r) for r in R if
        (!onlyreal || isreal(r, realtol)) &&
        (!onlynonsingular || isnonsingular(r, singulartol)) &&
        (!onlysingular || issingular(r, singulartol)) &&
        (!onlyfinite || isfinite(r) || isprojective(r))]
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
    realsolutions(result; tol=1e-6, conditions...)

Return all real solution (as `Vector`s of reals) for which the given conditions apply.
For the possible `conditions` see [`results`](@ref). Note that `onlyreal` is always `true`
and `realtol` is now `tol`.

## Example
```julia
julia> @polyvar x y
julia> result = solve([(x-2)y, y+x+3]);
julia> realsolutions(result)
[[2.0, -5.0], [-3.0, 0.0]]
```
"""
function realsolutions(result::Results; onlyreal=true, tol=1e-6, kwargs...)
    mapresults(r -> real.(solution(r)), result; onlyreal=true, realtol=tol, kwargs...)
end

"""
    nonsingular(result::Results; conditions...)

Return all `PathResult`s for which the solution is non-singular. This is just a shorthand
for `results(R; onlynonsingular=true, conditions...)`. For the possible `conditions` see [`results`](@ref).
"""
nonsingular(R::Results; kwargs...) = results(R; onlynonsingular=true, kwargs...)

"""
    singular(result::Results; conditions...)

Return all `PathResult`s for which the solution is singular. This is just a shorthand
for `results(R; onlysingular=true, conditions...)`. For the possible `conditions` see [`results`](@ref).
"""
function singular(R::Results; singulartol=1e14, tol=singulartol, kwargs...)
    results(R; onlysingular=true, singulartol=tol, kwargs...)
end


"""
    finite(result::AffineResults; conditions...)

Return all `PathResult`s for which the solution is finite. This is just a shorthand
for `results(R; onlyfinite=true, conditions...)`. For the possible `conditions` see [`results`](@ref).
"""
finite(R::Results; kwargs...) = results(R; onlyfinite=true, kwargs...)

"""
    real(result, tol=1e-6)

Get all results where the solutions are real with the given tolerance `tol`.
See [`isreal`](@ref) for details regarding the determination of 'realness'.
"""
Base.real(R::Results; tol=1e-6) = [r for r in R if isreal(r, tol)]

"""
    failed(result)

Get all results where the path tracking failed.
"""
failed(R::Results) = [r for r in R if isfailed(r)]

"""
    atinfinity(result::AffineResult)

Get all results where the solutions is at infinity.
"""
function atinfinity(R::AffineResults)
    [r for r in R if isatinfinity(r)]
end



"""
    multiplicities(V::Results; tol=1e-6)

Returns a `Vector` of `Vector{PathResult}`s grouping the `PathResult`s whose solutions appear with multiplicities *greater* 1 in 'V'.
Two solutions are regarded as equal, when their pairwise distance is less than 'tol'.
"""
function multiplicities(V::Results; tol=1e-6)
    output = Vector{Vector{PathResult}}()
    if all(v -> v.solution_type == :affine, V)
        M = multiplicities(solutions(V), infinity_norm, tol = tol)
        for m in M
            push!(output, V[m])
        end
    elseif all(v -> v.solution_type == :projective, V)
        M = multiplicities(LinearAlgebra.normalize.(solutions(V)), fubini_study, tol = tol)
        for m in M
            push!(output, V[m])
        end
    else
        @warn("Input contains both affine and projective data. Empty vector is returned.")
    end
    output
end

"""
    uniquesolutions(R::Result; tol=1e-6, multiplicities=false, conditions...)

Return all *unique* solutions. If `multiplicities` is `true`, then
all *unique* solutions with their correspnding multiplicities as pairs `(s, m)`
where `s` is the solution and `m` the multiplicity are returned.
For the possible `conditions` see [`results`](@ref).

## Example
```julia-repl
julia> @polyvar x;
julia> uniquesolutions([(x-3)^3*(x+2)], multiplicities=true)
[([3.0+0.0im], 3), ([-2.0+0.0im], 1)]
julia> uniquesolutions([(x-3)^3*(x+2)])
[[3.0+0.0im], [-2.0+0.0im]]
```
"""
function uniquesolutions(R::Results; tol=1e-6, multiplicities=false, conditions...)
    uniquesolutions(R, Val{multiplicities}; tol=tol, conditions...)
end

function uniquesolutions(R::Results, ::Type{Val{B}}; tol=1e-6, conditions...) where B
    sols = solutions(R; conditions...)
    if R isa AffineResult
        M = multiplicities(sols, infinity_norm, tol = tol)
    elseif R isa ProjectiveResult
        M = multiplicities(LinearAlgebra.normalize.(sols), fubini_study, tol = tol)
    end
    _uniquesolutions(sols, M, Val{B})
end

function _uniquesolutions(solutions, multiplicities, T::Type{Val{B}}) where B
    indicator = trues(length(solutions))
    uniques = map(multiplicities) do m
        for k in m
            indicator[k] = false
        end
        if T === Val{true}
            (solutions[m[1]], length(m))
        else
            solutions[m[1]]
        end
    end
    for (k, s) in enumerate(solutions)
        if indicator[k]
            if T == Val{true}
                push!(uniques, (s, 1))
            else
                push!(uniques, s)
            end
        end
    end
    uniques
end



####Show functions
function Base.show(io::IO, x::AffineResult)
    s = statistics(x)
    println(io, "AffineResult with $(length(x)) tracked paths")
    println(io, "==================================")
    println(io, "• $(s.nonsingular) non-singular finite $(plural("solution", s.nonsingular)) ($(s.real_nonsingular) real)")
    println(io, "• $(s.singular) singular finite $(plural("solution", s.singular)) ($(s.real_singular) real)")
    println(io, "• $(s.atinfinity) $(plural("solution", s.atinfinity)) at infinity")
    println(io, "• $(s.failed) failed $(plural("path", s.failed))")
    println(io, "• random seed: $(seed(x))")
end

function Base.show(io::IO, x::ProjectiveResult)
    s = statistics(x)
    println(io, "ProjectiveResult with $(length(x)) tracked paths")
    println(io, "==================================")
    println(io, "• $(s.nonsingular) non-singular $(plural("solution", s.nonsingular)) ($(s.real_nonsingular) real)")
    println(io, "• $(s.singular) singular finite $(plural("solution", s.singular)) ($(s.real_singular) real)")
    println(io, "• $(s.failed) failed $(plural("path", s.failed))")
    println(io, "• random seed: $(seed(x))")
end

TreeViews.hastreeview(::AffineResult) = true
TreeViews.hastreeview(::ProjectiveResult) = true
TreeViews.numberofnodes(::AffineResult) = 7
TreeViews.numberofnodes(::ProjectiveResult) = 6
TreeViews.treelabel(io::IO, x::AffineResult, ::MIME"application/prs.juno.inline") =
    print(io, "<span class=\"syntax--support syntax--type syntax--julia\">AffineResult</span>")
TreeViews.treelabel(io::IO, x::ProjectiveResult, ::MIME"application/prs.juno.inline") =
    print(io, "<span class=\"syntax--support syntax--type syntax--julia\">ProjectiveResult</span>")

function TreeViews.nodelabel(io::IO, x::AffineResult, i::Int, ::MIME"application/prs.juno.inline")
    s = statistics(x)
    if i == 1
        print(io, "Paths tracked")
    elseif i == 2 && s.nonsingular > 0
        print(io, "$(s.nonsingular) finite non-singular ($(s.real_nonsingular) real)")
    elseif i == 3 && s.singular > 0
        print(io, "$(s.singular) finite singular ($(s.real_singular) real)")
    elseif i == 4 && (s.real_nonsingular+s.real_singular) > 0
        print(io, "$(s.real_nonsingular+s.real_singular) finite real")
    elseif i == 5 && s.atinfinity > 0
        print(io, "$(s.atinfinity) atinfinity")
    elseif i == 6 && s.failed > 0
        print(io, "$(s.failed) failed")
    elseif i == 7
        print(io, "Random seed used")
    end
end

function TreeViews.treenode(r::AffineResult, i::Integer)
    s = statistics(r)
    if i == 1
        return length(r)
    elseif i == 2 && s.nonsingular > 0
        return finite(r, onlynonsingular=true)
    elseif i == 3 && s.singular > 0
        return finite(r, onlysingular=true)
    elseif i == 4 && (s.real_nonsingular+s.real_singular) > 0
        return finite(r, onlyreal = true)
    elseif i == 5 && s.atinfinity > 0
        return atinfinity(r)
    elseif i == 6 && s.failed > 0
        return failed(r)
    elseif i == 7
        return seed(r)
    end
    missing
end

function TreeViews.nodelabel(io::IO, x::ProjectiveResult, i::Int, ::MIME"application/prs.juno.inline")
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
    end
end

function TreeViews.treenode(r::ProjectiveResult, i::Integer)
    s = statistics(r)
    if i == 1
        return length(r)
    elseif i == 2 && s.nonsingular > 0
        return finite(r, onlynonsingular=true)
    elseif i == 3 && s.singular > 0
        return finite(r, onlysingular=true)
    elseif i == 4 && (s.real_nonsingular + s.real_singular) > 0
        return finite(r, onlyreal=true)
    elseif i == 5 && s.failed > 0
        return failed(r)
    elseif i == 6
        return seed(r)
    end
    missing
end

plural(singularstr, n) = n == 1 ? singularstr :  singularstr * "s"
