struct NewtonCorrector end

struct NewtonCorrectorCache{HC<:AbstractHomotopyCache,AV<:AbstractVector{ComplexDF64}}
    Δx::Vector{ComplexF64}
    r::Vector{ComplexF64}
    abs_r::Vector{Float64}
    # D64 evaluation
    x_D64::AV
    cache_D64::HC
end

function NewtonCorrectorCache(H::HomotopyWithCache, x::AbstractVector, t::Number)
    m, n = size(H)
    Δx = Vector{ComplexF64}(undef, n)
    r = Vector{ComplexF64}(undef, m)
    abs_r = Vector{Float64}(undef, m)
    x_D64 = similar(x, ComplexDF64)
    cache_D64 = cache(H.homotopy, x_D64, t)
    NewtonCorrectorCache(Δx, r, abs_r, x_D64, cache_D64)
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

struct NewtonCorrectorResult
    return_code::NewtonCorrectorCodes
    accuracy::Float64
    iters::Int
    ω₀::Float64
    ω::Float64
    θ₀::Float64
    θ::Float64
    norm_Δx₀::Float64
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
    tol::Float64 = 1e-6,
    min_iters::Int = 1,
    max_iters::Int = 3,
    debug::Bool = false,
    simplified_last_step::Bool = true,
    full_steps::Int = max_iters - simplified_last_step,
    double_64_evaluation::Bool = false,
)

    # Setup values
    x̄ .= x₀
    acc = limit_acc = norm_Δxᵢ = norm_Δxᵢ₋₁ = norm_Δx₀ = Inf
    ω = ω₀ = θ₀ = θ = NaN
    @unpack Δx, r, x_D64, cache_D64 = cache
    # alias to make logic easier
    xᵢ₊₁ = xᵢ = x̄
    is_converged = false
    for i ∈ 1:max_iters
        debug && println("i = ", i)

        compute_jacobian = i == 1 || i ≤ full_steps

        if compute_jacobian && double_64_evaluation
            x_D64 .= xᵢ
            evaluate!(r, H.homotopy, x_D64, t, cache_D64)
            jacobian!(jacobian(JM), H, xᵢ, t)
            updated!(JM)
        elseif compute_jacobian
            evaluate_and_jacobian!(r, jacobian(JM), H, xᵢ, t)
            updated!(JM)
        elseif double_64_evaluation
            x_D64 .= xᵢ
            evaluate!(r, H.homotopy, x_D64, t, cache_D64)
        else
            evaluate!(r, H, xᵢ, t)
        end

        LA.ldiv!(Δx, JM, r, norm)
        norm_Δxᵢ₋₁ = norm_Δxᵢ
        norm_Δxᵢ = Float64(norm(Δx))

        if isnan(norm_Δxᵢ)
            acc = norm_Δx₀ = norm_Δxᵢ
            return NewtonCorrectorResult(NEWT_TERMINATED, acc, i, ω₀, ω, θ₀, θ, norm_Δx₀)
        end

        debug && println("||Δxᵢ|| = ", norm_Δxᵢ)

        if i == 1
            xᵢ₊₁ .= xᵢ .- Δx
            acc = norm_Δx₀ = norm_Δxᵢ
            is_converged = acc ≤ tol
            if is_converged && i ≥ min_iters
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
        if isnan(θ) || θ ≥ 0.5
            code = is_converged ? NEWT_CONVERGED : NEWT_TERMINATED
            return NewtonCorrectorResult(code, acc, i, ω₀, ω, θ₀, θ, norm_Δx₀)
        end

        xᵢ₊₁ .= xᵢ .- Δx
        acc = norm_Δxᵢ / (1.0 - min(0.5, 2 * θ^2))
        is_converged = acc ≤ tol
        if is_converged && i ≥ min_iters
            return NewtonCorrectorResult(NEWT_CONVERGED, acc, i, ω₀, ω, θ₀, θ, norm_Δx₀)
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
        norm_Δx₀,
    )
end


"""
    limit_accuracy!(
        H::HomotopyWithCache,
        x::AbstractVector,
        t::Number,
        JM::JacobianMonitor,
        norm::AbstractNorm,
        cache::NewtonCorrectorCache;
        update_jacobian::Bool = false,
        compute_cond::Bool = false,
        compute_eval_err::Bool = false,
        x_accuracy::Float64 = 1e-8,
    )

Obtain an estimate of the limiting accuracy of Newton's method
and of the limiting residual.
If `update_jacobian == true` a new Jacobian is computed, otherwise the currently stored
Jacobian is used.
"""
function limit_accuracy!(
    H::HomotopyWithCache,
    x::AbstractVector,
    t::Number,
    JM::JacobianMonitor,
    norm::AbstractNorm,
    cache::NewtonCorrectorCache;
    update_jacobian::Bool = false,
    compute_cond::Bool = false,
    compute_eval_err::Bool = false,
    x_accuracy::Float64 = 1e-8,
)
    @unpack Δx, r = cache

    if update_jacobian
        evaluate_and_jacobian!(r, jacobian(JM), H, x, t)
        updated!(JM)
    else
        evaluate!(r, H, x, t)
    end
    LA.ldiv!(Δx, JM, r, norm)
    limit_acc = Float64(norm(Δx))

    if compute_cond || compute_eval_err
        cache.abs_r .= abs.(r)
    end

    compute_cond && cond!(JM, norm, cache.abs_r)

    eval_err = eps()
    if compute_eval_err
        row_scaling!(jacobian(JM), norm, cache.abs_r)
        apply_row_scaling!(cache.abs_r, JM)
        eval_err = max(eps(), maximum(cache.abs_r))
    end

    return (limit_acc = limit_acc, eval_err = eval_err)
end
