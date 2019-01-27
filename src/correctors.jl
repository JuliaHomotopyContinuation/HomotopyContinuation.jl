import ..AffinePatches
using ..Utilities

export AbstractCorrector,
    AbstractCorrectorCache,
    CorrectorResult,
    cache,
    correct!

# Interface
"""
Each corrector has to be a subtype of `AbstractCorrector`.
"""
abstract type AbstractCorrector end

"""
Each corrector has to construct a cache via `[cache](@ref)` which has
to be a subtype of `AbstractCorrectorCache`.
"""
abstract type AbstractCorrectorCache end

@enum ReturnCode begin
    converged
    terminated
    terminated_no_approximate
    maximal_iterations
end

"""
    Result{T, Corrector}

Structure holding information about a `correct!` step. The fields are
* `converged::Bool` Indicating whether the correction was successfull
* `res::T` The residual of the final result.
* `iters::Int` The number of iterations used.
"""
struct CorrectorResult{T}
    retcode::ReturnCode
    accuracy::T
    iters::Int
    ω₀::Float64
    ω::Float64
    norm_Δx₀::T
    cond::Float64
end

Base.show(io::IO, ::MIME"application/prs.juno.inline", r::CorrectorResult) = r
Base.show(io::IO, result::CorrectorResult) = print_fieldnames(io, result)


"""
    isconverged(result::Result)

Returns whether the correction was successfull.
"""
isconverged(result::CorrectorResult) = result.retcode == converged

"""
    cache(::AbstractCorrector, ::HomotopyWithCache{M, N}, x, t)::AbstractCorrectorCache

Construct a cache to avoid allocations.
"""
cache(c::AbstractCorrector, args...) = throw(MethodError(cache, tuple(c, args...)))


"""
    correct!(xnext, ::AbstractCorrector, ::AbstractCorrectorCache, H::HomotopyWithCache, x, t, tol)::Result

Perform a correction step such that in the end `H(xnext,t) < tol`.
Returns a [`Result`](@ref).
"""
function correct! end

include("correctors/newton.jl")
