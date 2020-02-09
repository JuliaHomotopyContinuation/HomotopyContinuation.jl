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
    x_high::Vector{ComplexDF64}
end

function NewtonCorrector(a::Float64, (m, n)::NTuple{2,Int})
    h_a = 2a * (√(4 * a^2 + 1) - 2a)
    Δx = zeros(ComplexF64, n)
    r = zeros(ComplexF64, m)
    r̄ = zero(r)
    x_high = zeros(ComplexDF64, n)
    NewtonCorrector(a, h_a, Δx, r, r̄, x_high)
end

function high_prec_refinement_step!(
    x̄::AbstractVector,
    NC::NewtonCorrector,
    H::AbstractHomotopy,
    x::AbstractVector,
    t::Number,
    JM::Jacobian,
    norm::AbstractNorm;
    high_precision::Bool = false,
)
    @unpack Δx, r, x_high = NC
    evaluate_and_jacobian!(r, jacobian(JM), H, x, t)
    x_high .= x
    evaluate!(r, H, x_high, ComplexDF64(t))
    LA.ldiv!(Δx, updated!(JM), r, norm)
    x̄ .= x .- Δx
    x_high .= x̄
    evaluate!(r, H, x_high, ComplexDF64(t))
    LA.ldiv!(Δx, JM, r, norm)
    norm(Δx)
end

function newton!(
    x̄::AbstractVector,
    NC::NewtonCorrector,
    H::AbstractHomotopy,
    x₀::AbstractVector,
    t::Number,
    JM::Jacobian,
    norm::AbstractNorm;
    μ::Float64 = throw(UndefKeywordError(:μ)),
    ω::Float64 = throw(UndefKeywordError(:ω)),
    high_precision::Bool = false,
)
    @unpack a, h_a, Δx, r, r̄, x_high = NC

    x̄ .= x₀
    xᵢ₊₃ = xᵢ₊₂ = xᵢ₊₁ = xᵢ = x̄
    Δxᵢ = Δx
    μ_low = θ = norm_Δxᵢ = norm_Δxᵢ₋₁ = norm_Δx₀ = NaN
    ā = a
    for i = 0:10
        evaluate_and_jacobian!(r, jacobian(JM), H, xᵢ, t)
        if high_precision
            x_high .= xᵢ
            evaluate!(r, H, x_high, ComplexDF64(t))
        end
        LA.ldiv!(Δxᵢ, updated!(JM), r, norm)
        if i > 0 && high_precision
            iterative_refinement!(Δxᵢ, JM, r, norm)
        end
        xᵢ₊₁ .= xᵢ .- Δxᵢ
        norm_Δxᵢ = norm(Δxᵢ)

        i == 0 && (norm_Δx₀ = norm_Δxᵢ)
        i == 1 && (ω = 2 * norm_Δxᵢ / norm_Δxᵢ₋₁^2)
        i >= 1 && (θ = norm_Δxᵢ / norm_Δxᵢ₋₁)

        if i >= 1 && θ > ā
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
            evaluate_and_jacobian!(r, jacobian(JM), H, xᵢ, t)
            updated!(JM)
            if high_precision
                LA.ldiv!(Δxᵢ, JM, r)
                μ_low = norm(Δxᵢ)
                x_high .= xᵢ
                evaluate!(r, H, x_high, ComplexDF64(t))
            end
            LA.ldiv!(Δxᵢ, JM, r)
            if high_precision
                iterative_refinement!(Δxᵢ, JM, r)
            end
            xᵢ₊₂ .= xᵢ₊₁ .- Δxᵢ
            μ = norm_Δxᵢ₊₁ = norm(Δxᵢ)

            # TODO: Sometimes a second iteration?
            if i == 0 && !high_precision
                evaluate!(r, H, xᵢ, t)
                LA.ldiv!(Δxᵢ, JM, r)
                μ = norm(Δxᵢ)
            end

            # TODO: NOT IN PAPER SO FAR
            if i == 0
                ω̄ = 2 * norm_Δxᵢ / norm_Δxᵢ₊₁^2
                if ω̄ < ω
                    ω = ω̄
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
    JM::Jacobian,
    norm::WeightedNorm;
    a::Float64 = throw(UndefKeywordError(:a)),
)
    x₂ = x₁ = x̄ # alias to make logic easier
    @unpack a, Δx, r = NC

    evaluate_and_jacobian!(r, jacobian(JM), H, x₀, t)
    LA.ldiv!(Δx, updated!(JM), r, norm)

    v = norm(Δx) + eps()
    valid = false
    ω = μ = NaN
    ε = sqrt(v)
    for k = 1:3
        x̄ .= x₀ .+ ε .* weights(norm)

        evaluate_and_jacobian!(r, jacobian(JM), H, x̄, t)
        LA.ldiv!(Δx, updated!(JM), r, norm)
        x₁ .= x̄ .- Δx
        norm_Δx₀ = norm(Δx)
        evaluate!(r, H, x₁, t)
        LA.ldiv!(Δx, JM, r, norm)
        x₂ .= x₁ .- Δx
        norm_Δx₁ = norm(Δx) + eps()
        if norm_Δx₁ < a * norm_Δx₀
            ω = 2 * norm_Δx₁ / norm_Δx₀^2
            μ = norm_Δx₁
            valid = true
            break
        else
            ε *= sqrt(ε)
        end
    end

    (valid = valid, ω = ω, μ = μ)
end
