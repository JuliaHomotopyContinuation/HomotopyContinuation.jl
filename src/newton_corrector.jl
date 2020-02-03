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
    norm_Δx₀::Float64
end

Base.show(
    io::IO,
    ::MIME"application/prs.juno.inline",
    r::NewtonCorrectorResult,
) = r
Base.show(io::IO, result::NewtonCorrectorResult) = print_fieldnames(io, result)
is_converged(R::NewtonCorrectorResult) = R.return_code == NEWT_CONVERGED

struct NewtonCorrector
    a::Float64
    # derived
    h_a::Float64
    Δx::Vector{ComplexF64}
    r::Vector{ComplexF64}
end

function NewtonCorrector(a::Float64, (m, n)::NTuple{2,Int})
    h_a = 2a * (√(4 * a^2 + 1) - 2a)
    Δx = zeros(ComplexF64, n)
    r = zeros(ComplexF64, m)
    NewtonCorrector(a, h_a, Δx, r)
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
)
    @unpack a, h_a, Δx, r = NC

    x̄ .= x₀
    xᵢ₊₃ = xᵢ₊₂ = xᵢ₊₁ = xᵢ = x̄
    Δxᵢ = Δx
    θ = norm_Δxᵢ = norm_Δxᵢ₋₁ = norm_Δx₀ = NaN
    ā = a
    for i = 0:10
        evaluate_and_jacobian!(r, jacobian(JM), H, xᵢ, t)
        LA.ldiv!(Δxᵢ, updated!(JM), r, norm)
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
                norm_Δx₀,
            )
        elseif ω * norm_Δxᵢ^2 < 2 * μ * sqrt(1 - 2 * h_a)
            # TODO: Can we sometimes skip this?

            evaluate_and_jacobian!(r, jacobian(JM), H, xᵢ, t)
            LA.ldiv!(Δxᵢ, updated!(JM), r)
            xᵢ₊₂ .= xᵢ₊₁ .- Δxᵢ
            norm_Δxᵢ₊₁ = norm(Δxᵢ)

            # TODO: Sometimes a second iteration?
            xᵢ₊₂ .= xᵢ₊₁ .- Δxᵢ
            μ = norm(Δxᵢ)

            # TODO: NOT IN PAPER SO FAR
            if i == 0
                ω̄ = 2norm_Δxᵢ / norm_Δxᵢ₊₁^2
                if !isnan(ω̄) && ω̄ < ω
                    ω = ω̄
                else
                    ω = 0.5 * ω
                end
            end

            return NewtonCorrectorResult(
                NEWT_CONVERGED,
                μ,
                i + 2,
                ω,
                θ,
                norm_Δx₀,
            )
        end

        norm_Δxᵢ₋₁ = norm_Δxᵢ
        i >= 1 && (ā *= ā)
    end

    return NewtonCorrectorResult(NEWT_MAX_ITERS, μ, 11, ω, θ, norm_Δx₀)
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
        norm_Δx₁ = norm(Δx)

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
