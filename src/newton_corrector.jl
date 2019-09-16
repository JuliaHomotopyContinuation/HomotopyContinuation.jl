struct NewtonCorrector end

struct NewtonCorrectorCache{T}
    Δx::Vector{Complex{T}}
    r::Vector{Complex{T}}
end

function NewtonCorrectorCache(H::HomotopyWithCache, x::AbstractVector, t::Number)
    rx = evaluate(H, x, t)
    T = complex(float(eltype(rx)))
    Δx = Vector{T}(undef, length(x))
    r = Vector{T}(undef, length(rx))
    NewtonCorrectorCache(Δx, r)
end

cache(::NewtonCorrector, H, x, t) = NewtonCorrectorCache(H, x, t)

@doc """
    NewtonCorrectorReturnCode

The possible return codes of Newton's method

* `NEWT_CONVERGED`
* `NEWT_TERMINATED`
* `NEWT_MAX_ITERS`
""" @enum NewtonCorrectorCodes begin
    NEWT_CONVERGED
    NEWT_TERMINATED
    NEWT_MAX_ITERS
end

struct NewtonCorrectorResult{T}
    return_code::NewtonCorrectorCodes
    accuracy::T
    iters::Int
    ω₀::Float64
    ω::Float64
    θ₀::Float64
    θ::Float64
    norm_Δx₀::T
end

Base.show(io::IO, ::MIME"application/prs.juno.inline", r::NewtonCorrectorResult) = r
Base.show(io::IO, result::NewtonCorrectorResult) = print_fieldnames(io, result)
is_converged(R::NewtonCorrectorResult) = R.return_code == NEWT_CONVERGED

function newton!(
    x̄::AbstractVector,
    H::HomotopyWithCache,
    x₀::AbstractVector,
    t::Number,
    JM::JacobianMonitor,
    norm::AbstractNorm,
    cache::NewtonCorrectorCache;
    tol::Float64 = 1e-6, max_iters::Int = 3, debug::Bool = false
)

    # Setup values
    x̄ .= x₀
    acc = limit_acc = norm_Δxᵢ = norm_Δxᵢ₋₁ = norm_Δx₀ = Inf
    ω = ω₀ = θ₀ = θ = NaN
    @unpack Δx, r = cache
    # alias to make logic easier
    xᵢ₊₁ = xᵢ = x̄

    for i ∈ 1:max_iters
        debug && println("i = ", i)

        compute_jacobian = i == 1 || i < max_iters
        if compute_jacobian
            evaluate_and_jacobian!(r, jacobian(JM), H, xᵢ, t)
            updated!(JM)
        else
            evaluate!(r, H, xᵢ, t)
        end

        # Update cond info etc?
        LA.ldiv!(Δx, JM, r, norm, JAC_MONITOR_UPDATE_NOTHING)

        norm_Δxᵢ₋₁ = norm_Δxᵢ
        norm_Δxᵢ = Float64(norm(Δx))

        debug && println("||Δxᵢ|| = ", norm_Δxᵢ)

        xᵢ₊₁ .= xᵢ .- Δx

        if i == 1
            acc = norm_Δx₀ = norm_Δxᵢ
            if acc ≤ tol
                return NewtonCorrectorResult(NEWT_CONVERGED, acc, i, ω₀, ω, θ₀, θ, norm_Δx₀)
            end
            continue
        end

        if i == 2
            θ = θ₀ = norm_Δxᵢ / norm_Δxᵢ₋₁
            ω = ω₀ = 2θ / norm_Δxᵢ₋₁
        elseif i > 2
            θ = norm_Δxᵢ / norm_Δxᵢ₋₁
            ω = max(ω, 2θ / norm_Δxᵢ₋₁)
        end
        acc = norm_Δxᵢ / (1.0 - min(0.5, 2θ^2))
        if acc ≤ tol
            return NewtonCorrectorResult(NEWT_CONVERGED, acc, i, ω₀, ω, θ₀, θ, norm_Δx₀)
        end
        if θ ≥ 0.5
            return NewtonCorrectorResult(NEWT_TERMINATED, acc, i, ω₀, ω, θ₀, θ, norm_Δx₀)
        end
    end

    return return NewtonCorrectorResult(
        NEWT_MAX_ITERS,
        acc,
        max_iters,
        ω₀,
        ω,
        θ₀,
        θ,
        norm_Δx₀
    )
end


"""
    limit_accuracy!(
        res::AbstractVector{<:Real},
        H::HomotopyWithCache, x t,
        JM::JacobianMonitor, norm, cache::NewtonCorrectorCache;
        accuracy::Float64 = 1e-6,
        safety_margin::Float64 = 1e2,
        update_jacobian::Bool = false)

Obtain an estimate of the limiting accuracy of Newton's method
and of the limiting residual.
If `update_jacobian == true` a new Jacobian is computed, otherwise the currently stored
Jacobian is used.
The method first produces a single (simplified) Newton update `Δx`.
If `accuracy / norm(Δx) > safety_margin` then `norm(Δx)` is obtained.
Otherwise a second simplified Newton update is produced and the norm of this is returned.
The reason for the two step procedure is that the first update does not necessarily yield
the limiting accuracy, but if it is sufficiently good then there is no need to obtain the
actual limiting accuracy.
"""
function limit_accuracy!(
    limit_res::AbstractVector{<:Real},
    H::HomotopyWithCache,
    x::AbstractVector,
    t::Number,
    JM::JacobianMonitor,
    norm::AbstractNorm,
    cache::NewtonCorrectorCache;
    accuracy::Float64 = 1e-6, safety_margin::Float64 = 1e2, update_jacobian::Bool = false
)
    @unpack Δx, r = cache

    if update_jacobian
        evaluate_and_jacobian!(r, jacobian(JM), H, x, t)
        updated!(JM)
    else
        evaluate!(r, H, x, t)
    end
    # limit residual ≈ abs.(r) is
    limit_res .= abs.(r)
    
    # Update cond info etc?
    LA.ldiv!(Δx, JM, r, norm, JAC_MONITOR_UPDATE_NOTHING)
    x .-= Δx

    limit_acc = Float64(norm(Δx))

    # accuracy is not yet sufficient, maybe this is just an artifact
    # of the Newton step. So let's do another simplified Newton step
    if accuracy / limit_acc < safety_margin
        evaluate!(r, H, x, t)
        # Update cond info etc?
        LA.ldiv!(Δx, JM, r, norm, JAC_MONITOR_UPDATE_NOTHING)
        x .-= Δx
        limit_acc = Float64(norm(Δx))
    end


    return limit_acc
end
