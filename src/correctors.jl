export AbstractCorrector,
    AbstractCorrectorCache,
    CorrectorResult,
    cache,
    correct!,
    NewtonCorrector

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

module CorrectorReturnCode
    import ..NewtonReturnCode

    @enum codes begin
        converged
        terminated
        terminated_no_approximate
        maximal_iterations
    end

    function code(c::NewtonReturnCode.codes)
        if c == NewtonReturnCode.converged
            converged
        elseif c == NewtonReturnCode.terminated
            terminated
        elseif c == NewtonReturnCode.terminated_no_approximate
            terminated_no_approximate
        else # maximal_iterations
            maximal_iterations
        end
    end
end


"""
    CorrectorResult{T}

Structure holding information about a `correct!` step. The fields are
* `converged::Bool` Indicating whether the correction was successfull
* `res::T` The residual of the final result.
* `iters::Int` The number of iterations used.
"""
struct CorrectorResult{T}
    return_code::CorrectorReturnCode.codes
    accuracy::T
    iters::Int
    ω₀::Float64
    ω::Float64
    norm_Δx₀::T
end

Base.show(io::IO, ::MIME"application/prs.juno.inline", r::CorrectorResult) = r
Base.show(io::IO, result::CorrectorResult) = print_fieldnames(io, result)


"""
    isconverged(result::Result)

Returns whether the correction was successfull.
"""
isconverged(result::CorrectorResult) = result.return_code == CorrectorReturnCode.converged

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

# Newton Corrector
"""
    NewtonCorrector(;simplified_last_step=true)

An ordinary Newton's method. If `simplified_last_step` is `true`, then for the last iteration
the previously Jacobian will be used. This uses an LU-factorization for square systems
and a QR-factorization for overdetermined.
"""
struct NewtonCorrector <: AbstractCorrector
    simplified_last_step::Bool
end

NewtonCorrector(;simplified_last_step=true) = NewtonCorrector(simplified_last_step)


struct NewtonCorrectorCache{FH<:FixedHomotopy, NC<:AbstractNewtonCache} <: AbstractCorrectorCache
    F::FH
    C::NC
end

function cache(::NewtonCorrector, H::HomotopyWithCache, x, t)
    F = FixedHomotopy(H, t)
    C = newton_cache(F, x)

    NewtonCorrectorCache(F, C)
end


@inline function correct!(out, alg::NewtonCorrector, cache::NewtonCorrectorCache,
                  H::HomotopyWithCache, x, t, norm, jacobian::Jacobian, tol::Float64,
                  maxiters::Int, step_size_model;
                  precision::PrecisionOption=PRECISION_FIXED_64,
                  update_jacobian_infos::Bool=false,
                  use_qr::Bool=false)
    cache.F.t = t
    result = newton!(out, cache.F, x, norm, cache.C, jacobian, tol, 1, maxiters,
        alg.simplified_last_step, update_jacobian_infos, use_qr, step_size_model.ω,
        step_size_model.expected_Δx₀, precision)
    CorrectorResult(result)
end

function CorrectorResult(R::NewtonResult)
    CorrectorResult(CorrectorReturnCode.code(R.return_code),
                    R.accuracy, R.iters, R.ω₀, R.ω, R.norm_Δx₀,)
end
