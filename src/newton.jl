export NewtonResult, newton, newton!, NewtonCache, NewtonReturnCode


module NewtonReturnCode
    @doc """
        NewtonReturnCode.codes

    The possible return codes of Newton's method

    * `NewtonReturnCode.converged`
    * `NewtonReturnCode.terminated`
    * `NewtonReturnCode.terminated_no_approximate`
    * `NewtonReturnCode.maximal_iterations`
    """
    @enum codes begin
        converged
        terminated
        terminated_no_approximate
        maximal_iterations
    end
end


"""
    NewtonResult{T}

Structure holding information about the outcome of the `newton` function. The fields are.
* `return_code::NewtonReturnCode.codes`: The return code of computation. `NewtonReturnCode.converged` means that `accuracy ≤ tol`.
* `accuracy::T`: |xᵢ-xᵢ₋₁| for i = iters and x₀,x₁,…,xᵢ₋₁,xᵢ are the Newton iterates.
* `iters::Int`: The number of iterations used.
* `digits_lost::Float64` Estimate of the (relative) lost digits in the linear algebra.
"""
struct NewtonResult{T}
    return_code::NewtonReturnCode.codes
    accuracy::T
    iters::Int
    ω₀::Float64
    ω::Float64
    norm_Δx₀::T
end

Base.show(io::IO, ::MIME"application/prs.juno.inline", r::NewtonResult) = r
Base.show(io::IO, result::NewtonResult) = print_fieldnames(io, result)


"""
    isconverged(result::Result)

Returns whether the correction was successfull.
"""
isconverged(result::NewtonResult) = result.return_code == NewtonReturnCode.converged


abstract type AbstractNewtonCache end
"""
    NewtonCache(F::AbstractSystem, x)

Cache for the [`newton`](@ref) function.
"""
struct NewtonCache{T, SC<:AbstractSystemCache} <: AbstractNewtonCache
    Jac::Jacobian{T}
    rᵢ::Vector{T}
    Δxᵢ::Vector{T}
    system_cache::SC
end

struct NewtonCacheF64{SC<:AbstractSystemCache, SC2<:AbstractSystemCache, AV<:AbstractVector{Complex{Double64}}} <: AbstractNewtonCache
    Jac::Jacobian{ComplexF64}
    xᵢ_D64::AV
    rᵢ::Vector{ComplexF64}
    rᵢ_D64::Vector{Complex{Double64}}
    Δxᵢ::Vector{ComplexF64}
    system_cache::SC
    system_cache_D64::SC2
end

function newton_cache(F::AbstractSystem, x)
    system_cache = cache(F, x)
    Jac = Jacobian(Random.randn!(jacobian(F, x, system_cache)))
    rᵢ = evaluate(F, x, system_cache)
    Δxᵢ = zeros(eltype(x), length(x))

    NewtonCache(Jac, rᵢ, Δxᵢ, system_cache)
end

function newton_cache(F::AbstractSystem, x::AbstractVector{ComplexF64})
    system_cache = cache(F, x)
    Jac = Jacobian(Random.randn!(complex.(float.(jacobian(F, x, system_cache)))))
    rᵢ = evaluate(F, x, system_cache)
    xD64 = similar(x, Complex{Double64})
    system_cache_D64 = cache(F, xD64)
    rᵢ_D64 = similar(rᵢ, Complex{Double64})
    Δxᵢ = zeros(eltype(x), length(x))

    NewtonCacheF64(Jac, xD64, rᵢ, rᵢ_D64, Δxᵢ, system_cache, system_cache_D64)
end

"""
    newton(F::AbstractSystem, x₀, norm=euclidean_norm, cache=NewtonCache(F, x₀); tol=1e-6, miniters=1, maxiters=3, simplified_last_step=true)

An ordinary Newton's method. If `simplified_last_step` is `true`, then for the last iteration
the previously Jacobian will be used. This uses an LU-factorization for square systems
and a QR-factorization for overdetermined.
"""
function newton(F::AbstractSystem, x₀, norm=euclidean_norm, cache=NewtonCache(F, x₀); kwargs...)
    x = copy(x₀)
    x, newton!(x, F, x₀, norm, cache; kwargs...)
end

"""
    newton!(out, F::AbstractSystem, x₀, norm, cache::AbstractNewtonCache; tol=1e-6, miniters=1, maxiters=3, simplified_last_step=true)

In-place version of [`newton`](@ref). Needs a [`NewtonCache`](@ref) and `norm` as input.
"""
@inline function newton!(out, F::AbstractSystem, x₀, norm, cache::AbstractNewtonCache;
                 tol::Float64=1e-6, miniters::Int=1, maxiters::Int=3,
                 simplified_last_step::Bool=true,
                 update_jacobian_infos::Bool=false,
                 use_qr::Bool=false,
                 ω::Float64=0.0, expected_Δx₀::Float64=0.0,
                 precision::PrecisionOption=PRECISION_ADAPTIVE)
    newton!(out, F, x₀, norm, cache, cache.Jac, tol,
            miniters, maxiters, simplified_last_step,
            update_jacobian_infos, use_qr, ω, expected_Δx₀, precision)
end

function newton!(out, F::AbstractSystem, x₀, norm, cache::AbstractNewtonCache, Jac::Jacobian,
                 tol::Float64, miniters::Int, maxiters::Int, simplified_last_step::Bool,
                 update_jacobian_infos::Bool, use_qr::Bool, ω::Float64,
                 expected_Δx₀::Float64, precision::PrecisionOption)
    rᵢ, Δxᵢ = cache.rᵢ, cache.Δxᵢ
    Jᵢ = Jac.J
    copyto!(out, x₀)
    xᵢ₊₁ = xᵢ = out # just alias to make logic easier
    rᵢ₊₁ = rᵢ

    # initialize variables
    T = real(eltype(xᵢ))
    Θ₀ = Θᵢ₋₁ = norm_Δxᵢ₋₁ = norm_Δxᵢ = norm_Δx₀ = zero(T)
    accuracy = T(Inf)
    ω₀ = ω
    for i ∈ 0:maxiters
        curr_accuracy = i == 0 ? expected_Δx₀ : accuracy
        newton_step!(cache, Jac, F, xᵢ;
                        tol=tol,
                        curr_accuracy=curr_accuracy,
                        ω = ω,
                        simple_step=(i == maxiters && simplified_last_step),
                        use_qr=use_qr,
                        update_infos=(update_jacobian_infos && i == 0),
                        precision=precision)
        norm_Δxᵢ₋₁ = norm_Δxᵢ
        norm_Δxᵢ = Float64(norm(Δxᵢ))
        @inbounds for k in eachindex(xᵢ)
            xᵢ₊₁[k] = xᵢ[k] - Δxᵢ[k]
        end

        if i == 0
            accuracy = norm_Δx₀ = norm_Δxᵢ₋₁ = norm_Δxᵢ
            if norm_Δx₀ ≤ tol && i + 1 ≥ miniters
                return NewtonResult(NewtonReturnCode.converged, norm_Δx₀, i + 1, ω₀, ω, norm_Δx₀)
            end

        else
            Θᵢ₋₁ = norm_Δxᵢ / norm_Δxᵢ₋₁
            ωᵢ₋₁ = 2Θᵢ₋₁ / norm_Δxᵢ₋₁
            if i == 1
                ω = ω₀ = ωᵢ₋₁
            else
                ω = @fastmath max(ω, ωᵢ₋₁)
            end

            if Θᵢ₋₁ > 0.5
                return NewtonResult(NewtonReturnCode.terminated, accuracy, i + 1, ω₀, ω, norm_Δx₀)
            end

            accuracy = norm_Δxᵢ / (1 - 2Θᵢ₋₁^2)
            if accuracy ≤ tol && i + 1 ≥ miniters
                return NewtonResult(NewtonReturnCode.converged, accuracy, i + 1, ω₀, ω, norm_Δx₀)
            end
        end
    end

    return NewtonResult(NewtonReturnCode.maximal_iterations, accuracy, maxiters+1, ω₀, ω, norm_Δx₀)
end

@inline function newton_step!(cache::NewtonCache, Jac::Jacobian, F, xᵢ;
                            tol::Float64=1e-7, ω::Float64=NaN, curr_accuracy::Float64=NaN,
                            simple_step::Bool=false, update_infos::Bool=true,
                            use_qr::Bool=false,
                            precision::PrecisionOption=PRECISION_ADAPTIVE)
    rᵢ, Δxᵢ = cache.rᵢ, cache.Δxᵢ
    if simple_step
        evaluate!(rᵢ, F, xᵢ, cache.system_cache)
    else
        evaluate_and_jacobian!(rᵢ, Jᵢ, F, xᵢ, cache.system_cache)
        updated_jacobian!(Jac, update_infos=(update_infos || use_qr))
    end
    solve!(Δxᵢ, Jac, rᵢ, update_digits_lost=update_infos)
end

@inline function newton_step!(cache::NewtonCacheF64, Jac::Jacobian, F, xᵢ;
                    tol::Float64=1e-7, ω::Float64=NaN, curr_accuracy::Float64=NaN,
                    simple_step::Bool=false, update_infos::Bool=true,
                    use_qr::Bool=false,
                    precision::PrecisionOption=PRECISION_ADAPTIVE)
    Jᵢ = Jac.J
    @unpack rᵢ, rᵢ_D64, Δxᵢ, xᵢ_D64 = cache
    if simple_step
        if higher_precision(Jac, tol, curr_accuracy, ω, precision)
            xᵢ_D64 .= xᵢ
            evaluate!(rᵢ, F, xᵢ_D64, cache.system_cache_D64)
        else
            evaluate!(rᵢ, F, xᵢ, cache.system_cache)
        end
    else
        if higher_precision(Jac, tol, curr_accuracy, ω, precision)
            xᵢ_D64 .= xᵢ
            evaluate!(rᵢ, F, xᵢ_D64, cache.system_cache_D64)
            jacobian!(Jᵢ, F, xᵢ, cache.system_cache)
        else
            evaluate_and_jacobian!(rᵢ, Jᵢ, F, xᵢ, cache.system_cache)
        end
        updated_jacobian!(Jac, update_infos=(update_infos || use_qr))
    end
    solve!(Δxᵢ, Jac, rᵢ, update_digits_lost=update_infos)
end

@inline function higher_precision(Jac::Jacobian{ComplexF64}, tol::Float64,
                            curr_accuracy::Float64, ω::Float64,
                            precision::PrecisionOption)
    if precision == PRECISION_FIXED_64
        return false
    elseif precision == PRECISION_FIXED_128
        return true
    end
    digits_lost = unpack(Jac.digits_lost, 0.0)
    if !isfinite(curr_accuracy) || !isfinite(ω) || iszero(ω)
        digits_lost > 12 + log₁₀(tol)
    else
        digits_lost > 10 + log₁₀(max(0.5ω * curr_accuracy^2, tol))
    end
end
