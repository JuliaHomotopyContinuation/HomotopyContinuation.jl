export NewtonResult, newton, newton!, NewtonReturnCode

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
* `ω₀::Float64`: Lipschitz constant estimate for the first newton update.
* `ω::Float64`: Maximal lipschitz constant estimate over all newton updates.
* `norm_Δx₀::T`: Norm of the first newton update.
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


abstract type AbstractNewtonCache{T<:Number} end

"""
    NewtonCacheF64(F::AbstractSystem, x) <: AbstractNewtonCache

Cache for the [`newton`](@ref) function using `Float64`.
"""
struct NewtonCacheF64{SC<:AbstractSystemCache, SC2<:AbstractSystemCache, AV<:AbstractVector{Complex{Double64}}} <: AbstractNewtonCache{ComplexF64}
    JM::JacobianMonitor{Float64}
    xᵢ_D64::AV
    rᵢ::Vector{ComplexF64}
    rᵢ_D64::Vector{Complex{Double64}}
    Δxᵢ::Vector{ComplexF64}
    system_cache::SC
    system_cache_D64::SC2
end

newton_cache(F::AbstractSystem, x::AbstractVector) =  newton_cache(F, complex.(float.(x)))
function newton_cache(F::AbstractSystem, x::AbstractVector{ComplexF64})
    system_cache = cache(F, x)
    JM = JacobianMonitor(complex.(float.(jacobian(F, x, system_cache))))
    rᵢ = evaluate(F, x, system_cache)
    xD64 = similar(x, Complex{Double64})
    system_cache_D64 = cache(F, xD64)
    rᵢ_D64 = similar(rᵢ, Complex{Double64})
    Δxᵢ = zeros(eltype(x), length(x))

    NewtonCacheF64(JM, xD64, rᵢ, rᵢ_D64, Δxᵢ, system_cache, system_cache_D64)
end

"""
    newton(F::AbstractSystem, x₀, norm=InfNorm(), cache=NewtonCache(F, x₀); tol=1e-6, miniters=1, maxiters=3, simplified_last_step=true)

An ordinary Newton's method. If `simplified_last_step` is `true`, then for the last iteration
the previously Jacobian will be used. This uses an LU-factorization for square systems
and a QR-factorization for overdetermined.
"""
function newton(F::AbstractSystem, x₀, norm::AbstractNorm=InfNorm(), cache::AbstractNewtonCache{T}=newton_cache(F, x₀); kwargs...) where {T}
    x = Vector{T}(undef, length(x₀))
    x, newton!(x, F, x₀, norm, cache; kwargs...)
end

"""
    newton!(out, F::AbstractSystem, x₀, norm, cache::AbstractNewtonCache; tol=1e-6, miniters=1, maxiters=3, simplified_last_step=true)

In-place version of [`newton`](@ref). Needs a [`NewtonCache`](@ref) and `norm` as input.
"""
@inline function newton!(out, F::AbstractSystem, x₀, norm::AbstractNorm, cache::AbstractNewtonCache;
                 tol::Float64=1e-6, miniters::Int=1, maxiters::Int=3,
                 simplified_last_step::Bool=true,
                 jacobian_monitor_update::JacobianMonitorUpdates=JAC_MONITOR_UPDATE_NOTHING,
                 update_all_steps::Bool=false,
                 use_qr::Bool=false,
                 ω::Float64=0.0, expected_Δx₀::Float64=0.0,
                 precision::PrecisionOption=PRECISION_ADAPTIVE)
    newton!(out, F, x₀, norm, cache, cache.JM, tol,
            miniters, maxiters, simplified_last_step,
            jacobian_monitor_update, update_all_steps, ω, expected_Δx₀, precision)
end

function newton!(out, F::AbstractSystem, x₀, norm::AbstractNorm, cache::AbstractNewtonCache, JM::JacobianMonitor,
                 tol::Float64, miniters::Int, maxiters::Int, simplified_last_step::Bool,
                 jacobian_monitor_update::JacobianMonitorUpdates,
                 update_all_steps::Bool, ω::Float64,
                 expected_Δx₀::Float64, precision::PrecisionOption)
    copyto!(out, x₀)
    xᵢ₊₁ = xᵢ = out # just alias to make logic easier

    # initialize variables
    T = real(eltype(xᵢ))
    Θ₀ = Θᵢ₋₁ = norm_Δxᵢ₋₁ = norm_Δxᵢ = norm_Δx₀ = zero(T)
    accuracy = T(Inf)
    ω₀ = ω
    for i ∈ 0:maxiters
        curr_accuracy = i == 0 ? expected_Δx₀ : accuracy
        simple_step = i == maxiters && simplified_last_step
        Δxᵢ = newton_step!(cache, JM, F, xᵢ, norm, tol, curr_accuracy, ω, simple_step,
                     jacobian_monitor_update, precision)
        @show JM
        norm_Δxᵢ₋₁ = norm_Δxᵢ
        norm_Δxᵢ = Float64(norm(Δxᵢ))
        @inbounds for k in eachindex(xᵢ)
            xᵢ₊₁[k] = xᵢ[k] - Δxᵢ[k]
        end

        println("||⋅||: ", xᵢ₊₁ - [√2, -√2])
        if i == 0
            accuracy = norm_Δx₀ = norm_Δxᵢ₋₁ = norm_Δxᵢ
            @show norm_Δx₀
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

            @show norm_Δxᵢ Θᵢ₋₁
            if Θᵢ₋₁ > 0.25
                return NewtonResult(NewtonReturnCode.terminated, accuracy, i + 1, ω₀, ω, norm_Δx₀)
            end

            accuracy = norm_Δxᵢ / (1 - 2Θᵢ₋₁^2)
            @show accuracy
            if accuracy ≤ tol && i + 1 ≥ miniters
                return NewtonResult(NewtonReturnCode.converged, accuracy, i + 1, ω₀, ω, norm_Δx₀)
            end
        end
        # only update things in the first step
        if !update_all_steps
            jacobian_monitor_update = JAC_MONITOR_UPDATE_NOTHING
        end
    end

    return NewtonResult(NewtonReturnCode.maximal_iterations, accuracy, maxiters+1, ω₀, ω, norm_Δx₀)
end

@inline function newton_step!(cache::NewtonCacheF64, JM::JacobianMonitor, F, xᵢ,
                    norm::AbstractNorm, tol::Float64, ω::Float64, curr_accuracy::Float64,
                    simple_step::Bool, jac_monitor_update::JacobianMonitorUpdates,
                    precision::PrecisionOption)
    Jᵢ = jacobian(JM)
    @unpack rᵢ, rᵢ_D64, Δxᵢ, xᵢ_D64 = cache

    if simple_step
        if higher_precision(JM, tol, curr_accuracy, ω, precision)
            xᵢ_D64 .= xᵢ
            evaluate!(rᵢ, F, xᵢ_D64, cache.system_cache_D64)
        else
            evaluate!(rᵢ, F, xᵢ, cache.system_cache)
        end
    else
        if higher_precision(JM, tol, curr_accuracy, ω, precision)
            xᵢ_D64 .= xᵢ
            evaluate!(rᵢ, F, xᵢ_D64, cache.system_cache_D64)
            jacobian!(Jᵢ, F, xᵢ, cache.system_cache)
        else
            evaluate_and_jacobian!(rᵢ, Jᵢ, F, xᵢ, cache.system_cache)
        end
        updated!(JM)
    end
    @show rᵢ
    LA.ldiv!(Δxᵢ, JM, rᵢ, norm, jac_monitor_update)
    Δxᵢ
end

@inline function higher_precision(JM::JacobianMonitor{Float64}, tol::Float64,
                            curr_accuracy::Float64, ω::Float64,
                            precision::PrecisionOption)
    if precision == PRECISION_FIXED_64
        return false
    elseif precision == PRECISION_FIXED_128
        return true
    end

    if !isfinite(curr_accuracy) || !isfinite(ω) || iszero(ω)
        tol < forward_err(JM) * 1e4
    else
        max(0.5ω * curr_accuracy^2, tol) < forward_err(JM) * 1e4
    end
end
