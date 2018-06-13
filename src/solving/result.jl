export AffineResult, ProjectiveResult, nresults, nfinite, nsingular, natinfinity, nfailed, nsmooth,
    finite, results, failed, atinfinity, singular, smooth, seed

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

Base.start(r::Result) = start(r.pathresults)
Base.next(r::Result, state) = next(r.pathresults, state)
Base.done(r::Result, state) = done(r.pathresults, state)
Base.endof(r::Result) = endof(r.pathresults)
Base.iteratorsize(::Type{<:Result}) = Base.HasLength()
Base.eltype(r::Result) = eltype(r.pathresults)

"""
    nresults(result)

The number of proper solutions.
"""
nresults(R::Result) = count(issuccess, R)

"""
    nresults(affineresult)

The number of finite solutions.
"""
nfinite(R::AffineResult) = count(isfinite, R)

"""
    nsingular(result; tol=1e10)

The number of singular solutions. A solution is considered singular
if its windingnumber is larger than 1 or the condition number is larger than `tol`.
"""
nsingular(R::Result; tol=1e10) = count(R) do r
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
    smooth(result)

The number of smooth solutions.
"""
nsmooth(R::Result; tol = 1e10) = count(r -> issmooth(r, tol), R)

"""
    seed(result)

The random seed used in the computation.
"""
seed(result::Result) = result.seed

const AffineResults = Union{AffineResult, Vector{<:PathResult}}
const ProjectiveResults = Union{ProjectiveResult, Vector{<:PathResult}}
const Results = Union{Result, Vector{<:PathResult}}

# Filtering
"""
    results(result; only_real=false, realtol=1e-6, onlysmooth=true, singulartol=1e10, onlyfinite=true)

Return all `PathResult`s for which the given conditions apply.

    results(f::Function, result; kwargs...)

Additionally you can apply a transformation `f` on each result.

## Example

```julia
R = solve(F)

# This gives us all solutions considered real (but still as a complex vector).
realsolutions = results(solution, R, only_real=true)
"""
results(R::Result; kwargs...) = results(identity, R; kwargs...)
function results(f::Function, R::Result;
    onlyreal=false, realtol=1e-6, onlysmooth=false, singulartol=1e10,
    onlyfinite=true)
    [f(r) for r in R if
        (!onlyreal || isreal(r, realtol)) &&
        (!onlysmooth || issmooth(r, singulartol)) &&
        (!onlyfinite || isfinite(r) || isprojective(r))]
end

"""
    smooth(result::AffineResult)

Return all `PathResult`s for which the solution is smooth.
"""
smooth(R::Results; kwargs...) = smooth(identity, R; kwargs...)
function smooth(f::Function, R::Results; tol=1e10)
    [f(r) for r in R if issmooth(r, tol)]
end

"""
    finite(result::AffineResult)

Return all `PathResult`s for which the result is successfull and
the contained `solution` is indeed a solution of the system.

    finite(f::Function, result)

Additionally you can apply a transformation `f` on each result.
"""
finite(R::Results; kwargs...) = finite(identity, R; kwargs...)
function finite(f::Function, R::AffineResults; onlysmooth=false, tol=1e10)
    if !onlysmooth
        [f(r) for r in R if isfinite(r)]
    else
        [f(r) for r in R if isfinite(r) && issmooth(r, tol)]
    end
end


"""
    singular(result; tol=1e10)

Get all singular solutions. A solution is considered singular
if its windingnumber is larger than 1 or the condition number is larger than `tol`.
"""
singular(R::Result; tol=1e10) = [r for r in R if (issuccess(r) && issingular(r, tol))]


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
failed(R::Result) = [r for r in R if isfailed(r)]

"""
    atinfinity(result::AffineResult)

Get all results where the solutions is at infinity.
"""
function atinfinity(R::AffineResults)
    [r for r in R if isatinfinity(r)]
end

function Base.show(io::IO, r::AffineResult)
    println(io, "-----------------------------------------------")
    println(io, "Paths tracked: $(length(r))")
    println(io, "# smooth finite solutions:  $(nsmooth(r))")
    println(io, "# singular finite solutions:  $(nsingular(r))")
    println(io, "# solutions at infinity:  $(natinfinity(r))")
    println(io, "# failed paths:  $(nfailed(r))")
    println(io, "Random seed used: $(seed(r))")
    println(io, "-----------------------------------------------")
end

function Base.show(io::IO, r::ProjectiveResult)
    println(io, "-----------------------------------------------")
    println(io, "Paths tracked: $(length(r))")
    println(io, "# smooth solutions:  $(nsmooth(r))")
    println(io, "# singular solutions:  $(nsingular(r))")
    println(io, "# failed paths:  $(nfailed(r))")
    println(io, "Random seed used: $(seed(r))")
    println(io, "-----------------------------------------------")
end


function Juno.render(i::Juno.Inline, r::AffineResult)
        t = Juno.render(i, Juno.defaultrepr(r))
        t[:children] = [
            Juno.render(i, Text("Paths tracked → $(length(r))"))]
        #
        n = nsmooth(r)
        if n > 0
            t_result = Juno.render(i, finite(r))
            t_result[:head] = Juno.render(i, Text("$n finite smooth solutions"))
            push!(t[:children], t_result)
        end

        n = nsingular(r)
        if n > 0
            t_result = Juno.render(i, finite(r))
            t_result[:head] = Juno.render(i, Text("$n finite singular solutions"))
            push!(t[:children], t_result)
        end

        n = natinfinity(r)
        if n > 0
            t_result = Juno.render(i, atinfinity(r))
            t_result[:head] = Juno.render(i, Text("$n solutions at ∞"))
            push!(t[:children], t_result)
        end

        n = nfailed(r)
        if n > 0
            t_result = Juno.render(i, failed(r))
            t_result[:head] = Juno.render(i, Text("$n paths failed"))
            push!(t[:children], t_result)
        end
        push!(t[:children], Juno.render(i, Text("Random seed used → $(seed(r))")))

        return t
end

function Juno.render(i::Juno.Inline, r::ProjectiveResult)
        t = Juno.render(i, Juno.defaultrepr(r))
        t[:children] = [
            Juno.render(i, Text("Paths tracked → $(length(r))"))]
        #
        n = nsmooth(r)
        if n > 0
            t_result = Juno.render(i, smooth(r))
            t_result[:head] = Juno.render(i, Text("$n smooth solutions"))
            push!(t[:children], t_result)
        end

        n = nsingular(r)
        if n > 0
            t_result = Juno.render(i, singular(r))
            t_result[:head] = Juno.render(i, Text("$n singular solutions"))
            push!(t[:children], t_result)
        end

        n = nfailed(r)
        if n > 0
            t_result = Juno.render(i, failed(r))
            t_result[:head] = Juno.render(i, Text("$n paths failed"))
            push!(t[:children], t_result)
        end
        push!(t[:children], Juno.render(i, Text("Random seed used → $(seed(r))")))

        return t
end
