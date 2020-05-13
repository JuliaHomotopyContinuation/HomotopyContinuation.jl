"""
    Predict(H::AbstractHomotopy, ::AD{N}) where {N}

Construct a predictor. The predicor uses a Padé approximant of order (2,1) to produce
an initial guess for the Newton corrector.
For this, the predictor tries to compute the local derivatives up to order 4.
This is done using automatic differentiation for the derivatives up to order `N`.
Otherwise numerical differentiation is used.
If the numerical derivatives up to order 3 cannot be computed with at least some accuracy
the predictor switches to a Ralston's 2nd order Runge Kutta method.
"""
mutable struct Predictor{T<:AD}
    AD::T
    ND::NumericalDifferentiation
    trust_region::Float64
    local_error::Float64
    pade::Bool
    taylor::BitVector
    steps_since_pade_attempt::Int
    cond_H_ẋ::Float64
    # buffer for the computation of derivatives
    tx⁰::TaylorVector{1,ComplexF64}
    tx¹::TaylorVector{2,ComplexF64}
    tx²::TaylorVector{3,ComplexF64}
    tx³::TaylorVector{4,ComplexF64}
    tx⁴::TaylorVector{5,ComplexF64}
    trust_tx::Vector{Bool}
    u::Vector{ComplexF64}
    m₁::Vector{ComplexF64}
    m₂::Vector{ComplexF64}
    m₃::Vector{ComplexF64}
    y⁰::TaylorVector{1,ComplexF64}
end

function Predictor(@nospecialize(H::AbstractHomotopy), ad::AD)
    m, n = size(H)
    ND = NumericalDifferentiation(m, n)
    tx⁴ = TaylorVector{5}(ComplexF64, n)
    tx⁰ = TaylorVector{1}(tx⁴)
    tx¹ = TaylorVector{2}(tx⁴)
    tx² = TaylorVector{3}(tx⁴)
    tx³ = TaylorVector{4}(tx⁴)
    trust_tx = [true, true, true, true]
    u = zeros(ComplexF64, m)
    m₁ = zeros(ComplexF64, n)
    m₂ = zeros(ComplexF64, n)
    m₃ = zeros(ComplexF64, n)
    y⁰ = TaylorVector{1}(ComplexF64, n)

    Predictor(
        ad,
        ND,
        Inf,
        Inf,
        true,
        falses(n),
        0,
        1.0,
        tx⁰,
        tx¹,
        tx²,
        tx³,
        tx⁴,
        trust_tx,
        u,
        m₁,
        m₂,
        m₃,
        y⁰,
    )
end

function init!(predictor::Predictor)
    predictor.steps_since_pade_attempt = 0
    predictor.pade = true
    predictor.cond_H_ẋ = 1.0
    predictor
end

function update!(
    P::Predictor,
    H::AbstractHomotopy,
    x,
    t,
    jacobian::Jacobian,
    norm::AbstractNorm,
    x̂ = nothing, # the predicted value
    h = nothing; # the last step size
    dist_to_target::Float64,
)
    # If we cannot use a Pade approximant then don't try it in every step again, but
    # rather only every k steps
    if P.pade || P.steps_since_pade_attempt ≥ 3
        P.steps_since_pade_attempt = 0
        trusted = compute_derivatives!(
            P,
            H,
            x,
            t,
            jacobian,
            norm;
            dist_to_target = dist_to_target,
        )
    else
        P.steps_since_pade_attempt += 1
        trusted = compute_derivatives!(
            P,
            H,
            x,
            t,
            jacobian,
            norm;
            dist_to_target = dist_to_target,
            only_first = true,
        )
    end
    prev_order = P.pade ? 4 : 3
    if trusted ≥ 3 && update_pade!(P, norm)
        P.pade = true
    else
        P.local_error = Inf
        P.trust_region = Inf
        P.pade = false
    end

    if trusted < 4 && (isnothing(x̂) || isnothing(h))
        # just make a pessimistic guess
        P.local_error = 100.0
    elseif !(isnothing(x̂) || isnothing(h))
        P.local_error = min(P.local_error, norm(x, x̂) / h^prev_order)
        # x⁴ can be completely wrong, as fallback extrapolate the norm of x⁴ from x³ and τ
    elseif P.pade
        P.local_error = min(P.local_error, norm(last(vectors(P.tx³))) / P.trust_region)
    end

    P
end

function compute_derivatives!(
    predictor::Predictor,
    H::AbstractHomotopy,
    x,
    t,
    J::Jacobian,
    norm::AbstractNorm;
    dist_to_target::Float64,
    only_first::Bool = false,
)
    @unpack AD, ND, u, tx⁰, tx¹, tx², tx³, tx⁴, trust_tx = predictor
    x⁰, x¹, x², x³, x⁴ = vectors(tx⁴)

    τ = local_error(predictor)
    x⁰ .= x
    taylor!(u, Val(1), H, tx⁰, t, AD, ND; dist_to_target = dist_to_target)
    u .= .-u
    LA.ldiv!(x¹, J, u)
    # Check if we have to do iterative refinment for all the others as well
    δ = fixed_precision_iterative_refinement!(x¹, workspace(J), u, InfNorm())
    predictor.cond_H_ẋ = δ / eps()
    if δ > 1e-12
        δ̂ = iterative_refinement!(x¹, J, u, InfNorm(); tol = 1e-12, max_iters = 5)
        # @show δ δ̂
    end
    only_first && return 1

    trust_tx[2] = taylor!(
        u,
        Val(2),
        H,
        tx¹,
        t,
        AD,
        ND;
        dist_to_target = dist_to_target,
        incremental = true,
    )
    trust_tx[2] || (trust_tx[3] = trust_tx[4] = false; return 1)
    u .= .-u
    LA.ldiv!(x², J, u)
    δ > 1e-6 && iterative_refinement!(x², J, u, norm; tol = 1e-6, max_iters = 4)

    trust_tx[3] = taylor!(
        u,
        Val(3),
        H,
        tx²,
        t,
        AD,
        ND;
        dist_to_target = dist_to_target,
        incremental = true,
    )
    trust_tx[3] || (trust_tx[4] = false; return 2)

    u .= .-u
    LA.ldiv!(x³, J, u)
    δ > 1e-4 && iterative_refinement!(x³, J, u, norm; tol = 1e-4, max_iters = 3)

    trust_tx[4] = taylor!(
        u,
        Val(4),
        H,
        tx³,
        t,
        AD,
        ND;
        dist_to_target = dist_to_target,
        incremental = true,
    )
    trust_tx[4] || return 3
    u .= .-u
    LA.ldiv!(x⁴, J, u)

    return 4
end

function update_pade!(pred::Predictor{AD{M}}, norm::AbstractNorm) where {M}
    tx = pred.tx⁴
    local_err_vec = pred.m₁
    trust_tx = pred.trust_tx
    # This is an adaption of the algorithm outlined in
    # Robust Pad{\'e} approximation via SVD (Trefethen et al)
    τ = Inf
    local_error = Inf
    if M ≥ 3
        λ_min = 3.814697265625e-6 # exp2(-18)
        for i = 1:length(tx)
            x, x¹, x², x³, x⁴ = tx[i]
            # @show tx[i]
            c, c¹, c², c³ = abs.((x, x¹, x², x³))
            λ = exp(0.1 * (3 * log(c) + log(c¹) - log(c²) - 3 * log(c³)))
            tol = 1e-14 * c
            # @show i, c, λ
            if !isfinite(λ) || λ * c¹ ≤ tol || λ^2 * c² ≤ tol || λ^3 * c³ ≤ tol
                pred.taylor[i] = true
            else
                pred.taylor[i] = false
                τᵢ = c² / c³
                τᵢ < τ && (τ = τᵢ)
            end
            #
            # @show λ
            # @show c, λ * c¹, λ^2 * c¹, λ^3 * c³
            # # c¹ = abs(x¹)
            # c¹ = abs(x¹)
            # λ = max(λ_min, c¹)
            # c¹ /= λ
            # c² = abs(x²) / λ^2
            # absx³ = abs(x³)
            # c³ = absx³ / λ^3
            # tol = 1e-14 * max(c¹, c², c³)
            # if (c¹ ≤ tol && c² ≤ tol && c³ ≤ tol) || c² ≤ tol
            #     pred.taylor[i] = true
            # else
            #     pred.taylor[i] = false
            #     τᵢ = (c² / c³) / λ
            #     τᵢ < τ && (τ = τᵢ)
            # end
        end
    elseif trust_tx[3]
        # check that the ratio (x¹ / x²) / (x² / x³) doesn't differ by more than a
        # factor of 10 ≈ exp(2.3)
        d = pred.ND.logabs_norm[2] - 2 * pred.ND.logabs_norm[3] + pred.ND.logabs_norm[4]
        abs(d) > 2.3 && return false
        τ = exp(pred.ND.logabs_norm[3] - pred.ND.logabs_norm[4])
        τ2 = τ^2
        τ3 = τ2 * τ
        for i = 1:length(tx)
            x, x¹, x², x³, x⁴ = tx[i]
            tol = 1e-14 * abs(x)
            if abs(x²) * τ2 ≤ tol || abs(x³) * τ3 ≤ tol
                pred.taylor[i] = true
            else
                pred.taylor[i] = false
            end
        end
    end
    isnan(τ) && return false

    if trust_tx[4]
        for i = 1:length(tx)
            x, x¹, x², x³, x⁴ = tx[i]
            if pred.taylor[i]
                local_err_vec[i] = x⁴
            else
                e = x⁴ - x³ * (x³ / x²)
                if !isnan(e)
                    local_err_vec[i] = e
                end
            end
        end

        local_error = norm(local_err_vec)
    end

    pred.local_error = nanmin(local_error, Inf)
    pred.trust_region = nanmin(τ, Inf)

    true
end


function predict!(x̂, pred::Predictor, H, x, t, Δt, J::Jacobian)
    if pred.pade
        predict_pade!(x̂, pred, H, Δt)
    else
        predict_ralston2!(x̂, pred, H, x, t, Δt, J)
    end
    nothing
end

function predict_pade!(x̂, pred::Predictor, H, Δt)
    tx = pred.tx³
    for i = 1:length(tx)
        x, x¹, x², x³ = tx[i]
        if pred.taylor[i]
            x̂[i] = x + Δt * (x¹ + Δt * x²)
        else
            δᵢ = 1 - Δt * x³ / x²
            x̂[i] = x + Δt * (x¹ + Δt * x² / δᵢ)
        end
    end

    x̂
end

function predict_ralston2!(x̂, pred::Predictor, H, x, t, Δt, J::Jacobian)
    @unpack u, y⁰, m₁, m₂, m₃ = pred
    # 0 | 0 0 0
    m₁ .= last(vectors(pred.tx¹))

    # 2/3 | 2/3 0 0
    h = 2Δt / 3
    @. y⁰ = x + h * m₁
    x̂ .= first.(y⁰)
    evaluate_and_jacobian!(u, matrix(J), H, x̂, t + h)
    taylor!(u, Val(1), H, y⁰, t + h, pred.AD, pred.ND; dist_to_target = abs(Δt))
    u .= .-u
    LA.ldiv!(m₂, updated!(J), u)

    #   | 1/4 3/4
    @. x̂ = x + 0.25 * Δt * (m₁ + 3m₂)
    x̂
end

function predict_ralston4!(x̂, pred::Predictor, H, x, t, Δt, J::Jacobian)
    @unpack u, y⁰, m₁, m₂, m₃ = pred
    m₁ .= last(vectors(pred.tx¹))

    @. y⁰ = x + 0.4 * Δt * m₁
    x̂ .= first.(y⁰)
    evaluate_and_jacobian!(u, matrix(J), H, x̂, t + 0.4Δt)
    taylor!(u, Val(1), H, y⁰, t + 0.4Δt, pred.AD, pred.ND; dist_to_target = abs(Δt))
    u .= .-u
    LA.ldiv!(m₂, updated!(J), u)

    @. y⁰ = x + Δt * (0.29697761 * m₁ + 0.15875964 * m₂)
    x̂ .= first.(y⁰)
    h = 0.45573725 * Δt
    evaluate_and_jacobian!(u, matrix(J), H, x̂, t + h)
    taylor!(u, Val(1), H, y⁰, t + h, pred.AD, pred.ND; dist_to_target = abs(Δt))
    u .= .-u
    LA.ldiv!(m₃, updated!(J), u)

    m₄ = x̂
    @. y⁰ = x + Δt * (0.21810040 * m₁ - 3.05096516 * m₂ + 3.83286476m₃)
    x̂ .= first.(y⁰)
    evaluate_and_jacobian!(u, matrix(J), H, x̂, t + Δt)
    taylor!(u, Val(1), H, y⁰, t + Δt, pred.AD, pred.ND; dist_to_target = abs(Δt))
    u .= .-u
    LA.ldiv!(m₄, updated!(J), u)

    for i in eachindex(x̂)
        x̂[i] =
            x[i] +
            Δt * (
                0.17476028 * m₁[i] - 0.55148066 * m₂[i] +
                1.20553560 * m₃[i] +
                0.17118478 * m₄[i]
            )
    end

    x̂
end

order(P::Predictor) = P.pade ? 4 : 3
local_error(predictor::Predictor) = predictor.local_error
trust_region(predictor::Predictor) = predictor.trust_region
