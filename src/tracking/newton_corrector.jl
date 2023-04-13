@doc """
    NewtonCorrectorReturnCode

The possible return codes of Newton's method

* `NEWT_CONVERGED`
* `NEWT_TERMINATED`
* `NEWT_MAX_ITERS`
""" @enum NewtonCorrectorResultCodes begin
    NEWT_CONVERGED
    NEWT_TERMINATED
    NEWT_MAX_ITERS
    NEWT_SINGULARITY
end

struct NewtonCorrectorResult
    return_code::NewtonCorrectorResultCodes
    accuracy::Float64
    iters::Int
    ω::Float64
    θ::Float64
    μ_low::Float64
    norm_Δx₀::Float64
end

Base.show(io::IO, result::NewtonCorrectorResult) = print_fieldnames(io, result)
is_converged(R::NewtonCorrectorResult) = R.return_code == NEWT_CONVERGED

Base.@kwdef struct NewtonCorrector
    ws::LinearSolveWorkspace
    norm::WeightedNorm{InfNorm}
    r::Vector{ComplexF64}
    x_extended::Vector{ComplexDF64}
end
Base.show(io::IO, NC::NewtonCorrector) = print(io, "NewtonCorrector($(NC.ws)))")

function NewtonCorrector(ws::LinearSolveWorkspace, norm::WeightedNorm{InfNorm})
    m, n = size(ws)
    r = zeros(ComplexF64, m)
    x_extended = zeros(ComplexDF64, n)
    NewtonCorrector(ws, norm, r, x_extended)
end


function newton!(
    x̄::AbstractVector,
    NC::NewtonCorrector,
    H::AbstractHomotopy,
    x₀::AbstractVector,
    t::Number;
    max_iters::Int,
    tol::Float64,
    extended_precision::Bool = false,
)
    @unpack r, x_extended, ws, norm = NC
    x̄ .= x₀
    xᵢ₊₁ = xᵢ = x̄
    μ_low = θ = norm_Δxᵢ = norm_Δxᵢ₋₁ = norm_Δx₀ = NaN
    J = get_A(ws)
    N = max_iters
    μ = ω = NaN

    for iter = 0:(N-1)
        try
            evaluate_and_jacobian!(r, J, H, xᵢ, t)
            set_A!(ws, J)

            if extended_precision
                x_extended .= xᵢ
                evaluate!(r, H, x_extended, t)
            end
            set_b!(ws, r)
            solve!(ws)
            Δxᵢ = solution(ws)

            norm_Δxᵢ = norm(Δxᵢ)
            iter == 0 && (norm_Δx₀ = norm_Δxᵢ)

            if isnan(norm_Δxᵢ)
                return NewtonCorrectorResult(
                    NEWT_SINGULARITY,
                    norm_Δxᵢ,
                    iter + 1,
                    ω,
                    θ,
                    μ_low,
                    norm_Δx₀,
                )
            end

            xᵢ₊₁ .= xᵢ .- Δxᵢ

            iter == 1 && (ω = 2 * norm_Δxᵢ / norm_Δxᵢ₋₁^2)
            iter >= 1 && (θ = norm_Δxᵢ / norm_Δxᵢ₋₁)

            if norm_Δxᵢ ≤ tol
                evaluate_and_jacobian!(r, J, H, xᵢ, t)
                set_A!(ws, J)
                set_b!(ws, r)
                solve!(ws)
                Δxᵢ₊₁ = solution(ws)
                μ_low = norm(Δxᵢ₊₁)
                if extended_precision
                    x_extended .= xᵢ
                    evaluate!(r, H, x_extended, t)
                    set_b!(ws, r)
                    Δxᵢ₊₁ = solve!(ws)
                    μ = norm(Δxᵢ₊₁)
                else
                    μ = μ_low
                end
                xᵢ₊₁ .= xᵢ .- Δxᵢ₊₁

                return NewtonCorrectorResult(NEWT_CONVERGED, μ, iter, ω, θ, μ_low, norm_Δx₀)
            end
            norm_Δxᵢ₋₁ = norm_Δxᵢ
        catch
            return NewtonCorrectorResult(
                NEWT_SINGULARITY,
                norm_Δxᵢ,
                iter + 1,
                ω,
                θ,
                μ_low,
                norm_Δx₀,
            )

        end
    end

    return NewtonCorrectorResult(NEWT_MAX_ITERS, norm_Δxᵢ, N, ω, θ, μ_low, norm_Δx₀)
end
