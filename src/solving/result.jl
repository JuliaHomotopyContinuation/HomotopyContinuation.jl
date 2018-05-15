export nresults, nfinite, nsingular, natinfinity, nfailed,
    finite, results, failed, atinfinity

abstract type Result end

"""
    Result(Vector{PathResult})

Constructs a summary of `PathResult`s.  A `Result` contains:
* `PathResults` the vector of PathResults
* `tracked`: length of `PathResults`
* `finite`: Number of finite_results.
* `at_infinity`: Number of results at infinity.
* `singular`: Number of singular results.
* `failed`: Number of failed paths.
"""
struct AffineResult{T1, T2, T3} <: Result
    pathresults::Vector{PathResult{T1, T2, T3}}
end

Base.length(r::Result) = length(r.pathresults)
Base.getindex(r::Result, I) = getindex(r.pathresults, I)

Base.start(r::Result) = start(r.pathresults)
Base.next(r::Result, state) = next(r.pathresults, state)
Base.done(r::Result, state) = done(r.pathresults, state)
function Base.iteratorsize(::Type{AffineResult{T1,T2,T3}}) where {T1, T2, T3}
    Base.HasLength()
end
Base.eltype(r::Result) = eltype(r.pathresults)

"""
    nresults(result)

The number of proper solutions.
"""
nresults(R::AffineResult) = count(isfinite, R)

"""
    nresults(result)

The number of finite (isolated) solutions.
"""
nfinite(R::AffineResult) = count(isfinite, R)

"""
    nsingular(result; tol=1e10)

The number of singular solutions. A solution is considered singular
if its windingnumber is larger than 1 or the condition number is larger than `tol`.
"""
nsingular(R::Result; tol=1e10) = count(R) do r
    isfinite(r) && issingular(r, tol)
end

"""
    natinfinity(result)

The number of solutions at infinity.
"""
natinfinity(R::AffineResult) = count(isatinfinity, R)

"""
    natinfinity(result)

The number of failed paths.
"""
nfailed(R::Result) = count(isfailed, R)

const AffineResults = Union{AffineResult, Vector{<:PathResult}}
const Results = Union{Result, Vector{<:PathResult}}

# Filtering
"""
    results(result; only_real=false, realtol=1e-6,
        include_singular=true, singulartol=1e10,
        include_atinfinity=false)

Return all `PathResult`s for which the given conditions apply.

    results(f::Function, result; kwargs...)

Additionally you can apply a transformation `f` on each result.

## Example

```julia
R = solve(F)

# This gives us all solutions considered real (but still as a complex vector).
realsolutions = results(solution, R, only_real=true)
"""
results(R::Results; kwargs...) = results(identity, R; kwargs...)
function results(f::Function, R::AffineResult;
    only_real=false, realtol=1e-6, include_singular=true, singulartol=1e10,
    include_atinfinity=false)
    [f(r) for r in R if
        (!only_real || isreal(r, realtol)) &&
        (include_singular || !issingular(r, singulartol)) &&
        (include_atinfinity || isfinite(r))]
end


"""
    finite(result::AffineResult)

Return all `PathResult`s for which the result is successfull and
the contained `solution` is indeed a solution of the system.

    finite(f::Function, result)

Additionally you can apply a transformation `f` on each result.

## Example

```julia
R = solve(F)

# This gives us the actual solutions.
solutions = results(solution, R)

# We also want to get an overview over the condition numbers
condition_numbers = results(condition_number, R)
"""
finite(R::AffineResults; kwargs...) = finite(identity, R; kwargs...)
function finite(f::Function, R::AffineResults; include_singular=true, tol=1e10)
    if include_singular
        [f(r) for r in R if isfinite(r)]
    else
        [f(r) for r in R if isfiniter(r) && !issingular(r, tol)]
    end
end

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
atinfinity(R::Result) = [r for r in R if isatinfinity(r)]

function Base.show(io::IO, r::AffineResult)
    println(io, "-----------------------------------------------")
    println(io, "Paths tracked: $(length(r))")
    println(io, "# finite (isolated) solutions:  $(nfinite(r))")
    println(io, "# singular finite solutions:  $(nsingular(r))")
    println(io, "# solutions at infinity:  $(natinfinity(r))")
    println(io, "# failed paths:  $(nfailed(r))")
    println(io, "-----------------------------------------------")
end


function Juno.render(i::Juno.Inline, r::AffineResult)
        t = Juno.render(i, Juno.defaultrepr(r))
        t[:children] = [
            Juno.render(i, Text("Paths tracked → $(length(r))"))]
        #
        n = nfinite(r)
        if n > 0
            t_result = Juno.render(i, finite(r))
            t_result[:head] = Juno.render(i, Text("$n finite solutions"))
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

        return t
    end

#
#
