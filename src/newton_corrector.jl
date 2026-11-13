# This file contains the Newton corrector used for the path tracking in tracker.jl
# For a derivation see:
# Mixed Precision Path Tracking for Polynomial Homotopy Continuation,
# Sascha Timme (2020), arXiv:1902.02968

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
    NEWT_SINGULARITY
end

struct NewtonCorrectorResult
    return_code::NewtonCorrectorCodes
    accuracy::Float64
    iters::Int
    ω::Float64
    θ::Float64
    μ_low::Float64
    norm_Δx₀::Float64
end

Base.show(io::IO, ::MIME"application/prs.juno.inline", r::NewtonCorrectorResult) = r
Base.show(io::IO, result::NewtonCorrectorResult) = print_fieldnames(io, result)
is_converged(R::NewtonCorrectorResult) = R.return_code == NEWT_CONVERGED

struct NewtonCorrector
    a::Float64
    # derived
    h_a::Float64
    Δx::Vector{ComplexF64}
    r::Vector{ComplexF64}
    r̄::Vector{ComplexF64}
    x_extended::Vector{ComplexDF64}
    extended_workspace::Base.RefValue{Union{Nothing,ExtendedLinearWorkspace}}
    # Shared with the predictor through Tracker; reset when a new path starts.
    full_extended_precision::Base.RefValue{Bool}
end

function NewtonCorrector(a::Float64, x::AbstractVector{ComplexF64}, m::Int)
    n = length(x)
    h_a = 2a * (√(4 * a^2 + 1) - 2a)
    Δx = zeros(ComplexF64, n)
    r = zeros(ComplexF64, m)
    r̄ = zero(r)
    x_extended = ComplexDF64.(x)
    NewtonCorrector(
        a,
        h_a,
        Δx,
        r,
        r̄,
        x_extended,
        Ref{Union{Nothing,ExtendedLinearWorkspace}}(nothing),
        Ref(false),
    )
end

function extended_newton_correction!(NC, H, x, t)
    # Never leave a correction from a failed lower-precision solve available.
    fill!(NC.Δx, NaN)
    n = length(x)
    size(H) == (n, n) && extended_jacobian_supported(H) || return false
    W = NC.extended_workspace[]
    if W === nothing
        W = ExtendedLinearWorkspace(n)
        NC.extended_workspace[] = W
    end
    NC.x_extended .= x
    evaluate_and_jacobian!(W.residual, W.A, H, NC.x_extended, t)
    extended_factorize!(W)
    extended_ldiv!(W) || return false
    NC.Δx .= W.correction
    NC.full_extended_precision[] = true
    true
end

function full_precision_correction!(NC, H, x, t, J)
    ok = extended_newton_correction!(NC, H, x, t)
    J.factorizations[] += 1
    J.ldivs[] += 1
    if ok
        W = NC.extended_workspace[]
        matrix(J) .= W.A
        NC.r .= W.residual
        updated!(J)
    else
        fill!(NC.Δx, NaN)
    end
    ok
end

function refine_newton_correction!(NC, H, x, t, J, norm; tol)
    if NC.full_extended_precision[]
        full_precision_correction!(NC, H, x, t, J)
        return nothing
    end
    result = iterative_refinement!(NC.Δx, J, NC.r, norm; tol, max_iters = 3)
    if (result.diverged || !isfinite(result.accuracy) || result.accuracy > tol) &&
       (!all(iszero, NC.r) || !all(isfinite, NC.Δx))
        full_precision_correction!(NC, H, x, t, J)
    end
    nothing
end

function extended_prec_refinement_step!(
    x̄::AbstractVector,
    NC::NewtonCorrector,
    H::AbstractHomotopy,
    x::AbstractVector,
    t::Number,
    J::Jacobian,
    norm::AbstractNorm;
    simple_newton_step::Bool = true,
)
    @unpack Δx, r, x_extended = NC
    if NC.full_extended_precision[]
        full_precision_correction!(NC, H, x, t, J)
    else
        evaluate_and_jacobian!(r, matrix(J), H, x, t)
        x_extended .= x
        evaluate!(r, H, x_extended, t)
        LA.ldiv!(Δx, updated!(J), r, norm)
        refine_newton_correction!(NC, H, x, t, J, norm; tol = 1e-8)
    end
    x̄ .= x .- Δx
    if simple_newton_step
        if NC.full_extended_precision[]
            full_precision_correction!(NC, H, x̄, t, J)
        else
            x_extended .= x̄
            evaluate!(r, H, x_extended, t)
            LA.ldiv!(Δx, J, r, norm)
            NC.extended_workspace[] === nothing ||
                refine_newton_correction!(NC, H, x̄, t, J, norm; tol = 1e-8)
        end
    end
    norm(Δx)
end

function newton!(
    x̄::AbstractVector,
    NC::NewtonCorrector,
    H::AbstractHomotopy,
    x₀::AbstractVector,
    t::Number,
    J::Jacobian,
    norm::AbstractNorm;
    μ::Float64 = throw(UndefKeywordError(:μ)),
    ω::Float64 = throw(UndefKeywordError(:ω)),
    extended_precision::Bool = false,
    accurate_μ::Bool = false,
    first_correction::Bool = false,
)
    @unpack a, h_a, Δx, r, r̄, x_extended = NC

    x̄ .= x₀
    xᵢ₊₃ = xᵢ₊₂ = xᵢ₊₁ = xᵢ = x̄
    Δxᵢ = Δx
    μ_low = θ = norm_Δxᵢ = norm_Δxᵢ₋₁ = norm_Δx₀ = NaN
    ā = a
    for i = 0:10
        if extended_precision && NC.full_extended_precision[]
            full_precision_correction!(NC, H, xᵢ, t, J)
        else
            evaluate_and_jacobian!(r, matrix(J), H, xᵢ, t)
            if extended_precision
                copyto!(r̄, r)
                x_extended .= xᵢ
                evaluate!(r, H, x_extended, t)
            end
            LA.ldiv!(Δxᵢ, updated!(J), r, norm)
            if extended_precision
                refine_newton_correction!(NC, H, xᵢ, t, J, norm; tol = ā^2)
            end
        end
        norm_Δxᵢ = norm(Δxᵢ)

        if isnan(norm_Δxᵢ)
            @label return_singular
            return NewtonCorrectorResult(
                NEWT_SINGULARITY,
                norm(Δxᵢ),
                i + 1,
                ω,
                θ,
                μ_low,
                norm_Δx₀,
            )
        end

        # Detect a correction lost to rounding while applying it. If the point
        # is unchanged, its residual and Jacobian remain valid for the final
        # accuracy check; no new evaluation or factorization is necessary.
        unchanged = true
        @inbounds for k in eachindex(xᵢ, Δxᵢ)
            old = xᵢ[k]
            new = old - Δxᵢ[k]
            unchanged &= isequal(old, new)
            xᵢ₊₁[k] = new
        end

        i == 0 && (norm_Δx₀ = norm_Δxᵢ)
        i == 1 && (ω = 2 * norm_Δxᵢ / norm_Δxᵢ₋₁^2)
        i >= 1 && (θ = norm_Δxᵢ / norm_Δxᵢ₋₁)
        if i >= 1 && θ > ā || (i == 0 && !first_correction && 0.125 * norm_Δx₀ * ω > h_a)
            @label return_terminated
            return NewtonCorrectorResult(
                NEWT_TERMINATED,
                norm(Δxᵢ),
                i + 1,
                ω,
                θ,
                μ_low,
                norm_Δx₀,
            )
        elseif ω * norm_Δxᵢ^2 < 2 * μ * sqrt(1 - 2 * h_a)
            if extended_precision && NC.full_extended_precision[]
                μ_low = Inf
                unchanged || full_precision_correction!(NC, H, xᵢ, t, J)
            else
                if !unchanged
                    evaluate_and_jacobian!(r, matrix(J), H, xᵢ, t)
                    updated!(J)
                end
                if extended_precision
                    LA.ldiv!(Δxᵢ, J, unchanged ? r̄ : r)
                    μ_low = norm(Δxᵢ)
                    isfinite(μ_low) || (μ_low = Inf)
                    if !unchanged
                        x_extended .= xᵢ
                        evaluate!(r, H, x_extended, t)
                    end
                end
                LA.ldiv!(Δxᵢ, J, r)
                if extended_precision
                    refine_newton_correction!(NC, H, xᵢ, t, J, norm; tol = ā^2)
                end
            end
            xᵢ₊₂ .= xᵢ₊₁ .- Δxᵢ
            norm_Δxᵢ₊₁ = norm(Δxᵢ)
            if isnan(norm_Δxᵢ₊₁)
                @goto return_singular
                # it seems that the solution is singular, so we are diverging again when
                # we are getting tooo close
            elseif norm_Δxᵢ₊₁ > √norm_Δxᵢ
                θ = norm_Δxᵢ₊₁ / norm_Δxᵢ
                return NewtonCorrectorResult(
                    NEWT_TERMINATED,
                    norm_Δxᵢ₊₁,
                    i + 2,
                    ω,
                    θ,
                    μ_low,
                    norm_Δx₀,
                )
            end

            if norm_Δxᵢ₊₁ > 2μ && extended_precision
                if NC.full_extended_precision[]
                    full_precision_correction!(NC, H, xᵢ, t, J)
                else
                    x_extended .= xᵢ
                    evaluate!(r, H, x_extended, t)
                    LA.ldiv!(Δxᵢ, J, r)
                    refine_newton_correction!(NC, H, xᵢ, t, J, norm; tol = ā^2)
                end
                norm_Δxᵢ = norm_Δxᵢ₊₁
                μ = norm_Δxᵢ₊₁ = norm(Δxᵢ)
            elseif norm_Δxᵢ₊₁ > 2μ || accurate_μ
                if extended_precision && NC.full_extended_precision[]
                    full_precision_correction!(NC, H, xᵢ, t, J)
                else
                    evaluate!(r, H, xᵢ, t)
                    LA.ldiv!(Δxᵢ, J, r)
                end
                μ = norm(Δxᵢ)
            else
                μ = norm_Δxᵢ₊₁
            end

            if i == 0
                ω̄ = 2 * norm_Δxᵢ / norm_Δxᵢ₊₁^2
                if ω̄ < ω
                    ω = ω̄
                else
                    ω *= 0.25
                end
            end

            return NewtonCorrectorResult(NEWT_CONVERGED, μ, i + 2, ω, θ, μ_low, norm_Δx₀)
        end

        norm_Δxᵢ₋₁ = norm_Δxᵢ
        i >= 1 && (ā *= ā)
    end

    return NewtonCorrectorResult(NEWT_MAX_ITERS, μ, 11, ω, θ, μ_low, norm_Δx₀)
end

function init_newton!(
    x̄::AbstractVector,
    NC::NewtonCorrector,
    H::AbstractHomotopy,
    x₀::AbstractVector,
    t::Number,
    J::Jacobian,
    norm::WeightedNorm;
    a::Float64,
    extended_precision::Bool = true,
)
    x₂ = x₁ = x̄ # alias to make logic easier
    @unpack a, Δx, r, x_extended = NC

    if extended_precision && NC.full_extended_precision[]
        full_precision_correction!(NC, H, x₀, t, J)
    else
        evaluate_and_jacobian!(r, matrix(J), H, x₀, t)
        if extended_precision
            x_extended .= x₀
            evaluate!(r, H, x_extended, t)
        end
        LA.ldiv!(Δx, updated!(J), r, norm)
    end
    v = norm(Δx) + eps()
    if extended_precision && !isfinite(v)
        refine_newton_correction!(NC, H, x₀, t, J, norm; tol = 1e-8)
        v = norm(Δx) + eps()
    end
    valid = false
    ω = μ = NaN
    ε = sqrt(v)
    for k = 1:3
        x̄ .= x₀ .+ ε .* weights(norm)

        if extended_precision && NC.full_extended_precision[]
            full_precision_correction!(NC, H, x̄, t, J)
        else
            evaluate_and_jacobian!(r, matrix(J), H, x̄, t)
            if extended_precision
                x_extended .= x̄
                evaluate!(r, H, x_extended, t)
            end
            LA.ldiv!(Δx, updated!(J), r, norm)
        end

        if extended_precision && !all(isfinite, Δx)
            refine_newton_correction!(NC, H, x̄, t, J, norm; tol = 1e-8)
        end

        x₁ .= x̄ .- Δx
        norm_Δx₀ = norm(Δx)
        if extended_precision && NC.full_extended_precision[]
            full_precision_correction!(NC, H, x₁, t, J)
        else
            if extended_precision
                x_extended .= x₁
                evaluate!(r, H, x_extended, t)
            else
                evaluate!(r, H, x₁, t)
            end
            LA.ldiv!(Δx, J, r, norm)
        end
        if extended_precision && !all(isfinite, Δx)
            refine_newton_correction!(NC, H, x₁, t, J, norm; tol = 1e-8)
        end
        x₂ .= x₁ .- Δx
        norm_Δx₁ = norm(Δx) + eps()
        if norm_Δx₁ < a * norm_Δx₀
            ω = 2 * norm_Δx₁ / norm_Δx₀^2
            μ = norm_Δx₁
            if ω * μ > a^7
                refined_res = newton!(
                    x̄,
                    NC,
                    H,
                    x̄,
                    t,
                    J,
                    norm;
                    ω = ω,
                    μ = a^7 / ω,
                    accurate_μ = true,
                    extended_precision = extended_precision,
                )
                if is_converged(refined_res)
                    valid = true
                    ω = refined_res.ω
                    μ = refined_res.accuracy
                else
                    valid = false
                end
            else
                valid = true
                break
            end
        else
            ε *= sqrt(ε)
        end
    end

    (valid = valid, ω = ω, μ = μ)
end
