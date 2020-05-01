# This file contains a Padé predictor of order (2,1) used for the tracking in tracker.jl
# For a motivation of this see:
# Mixed Precision Path Tracking for Polynomial Homotopy Continuation,
# Sascha Timme (2020), arXiv:1902.02968

"""
    Predict(H::AbstractHomotopy, ::AD{N}) where {N}

Construct a predictor. The predicor uses a Padé approximant of order (2,1) to produce
an initial guess for the Newton corrector.
For this, the predictor tries to compute the local derivatives up to order 4.
This is done using automatic differentiation for the derivatives up to order `N`.
Otherwise numerical differentiation is used.
If the numerical derivatives up to order 3 cannot be computed with at least some accuracy
the predictor switches to a Raltson's 4th order Runge Kutta method.
"""
mutable struct Predictor{T<:AD}
    AD::T
    ND::NumericalDifferentiation
    trust_region::Float64
    local_error::Float64
    pade::Bool
    taylor::BitVector
    steps_since_pade_attempt::Int
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
    predictor
end

function update!(
    P::Predictor,
    H::AbstractHomotopy,
    x,
    t,
    jacobian::Jacobian,
    norm::AbstractNorm;
    dist_to_target::Float64,
)
    # If we cannot use a Pade approximant then don't try it in every step again, but
    # rather only every k steps
    if P.steps_since_pade_attempt < 5
        ord = compute_derivatives!(
            P,
            H,
            x,
            t,
            jacobian,
            norm;
            dist_to_target = dist_to_target,
        )
        P.steps_since_pade_attempt = 0
        if ord ≥ 3
            update_pade!(P, norm)
            P.pade = true
            return P
        end
    end
    P.steps_since_pade_attempt += 1
    P.local_error = Inf
    P.trust_region = Inf
    P.pade = false
    P
end

function compute_derivatives!(
    predictor::Predictor,
    H::AbstractHomotopy,
    x,
    t,
    jacobian::Jacobian,
    norm::AbstractNorm;
    dist_to_target::Float64,
    min_acc::Float64 = sqrt(eps()),
    max_iters::Int = 5,
)
    @unpack AD, ND, u, tx⁰, tx¹, tx², tx³, tx⁴, trust_tx = predictor
    x⁰, x¹, x², x³, x⁴ = vectors(tx⁴)

    x⁰ .= x
    taylor!(u, Val(1), H, tx⁰, t, AD, ND; cond = 1.0, dist_to_target = dist_to_target)
    u .= .-u
    LA.ldiv!(x¹, jacobian, u)

    # Check if we have to do iterative refinment for all the others as well
    δ = fixed_precision_iterative_refinement!(x¹, workspace(jacobian), u, norm)
    cond = δ / eps()
    iterative_refinement = δ > min_acc
    iterative_refinement &&
    iterative_refinement!(x¹, jacobian, u, norm; tol = min_acc, max_iters = max_iters)
    # if diverged
    #     state.code = TrackerCode.terminated_ill_conditioned
    # end

    trust_tx[2] = taylor!(
        u,
        Val(2),
        H,
        tx¹,
        t,
        AD,
        ND;
        cond = cond,
        dist_to_target = dist_to_target,
        incremental = true,
    )

    trust_tx[2] || return 1
    u .= .-u
    LA.ldiv!(x², jacobian, u)
    iterative_refinement &&
    iterative_refinement!(x², jacobian, u, norm; tol = min_acc, max_iters = max_iters)

    trust_tx[3] = taylor!(
        u,
        Val(3),
        H,
        tx²,
        t,
        AD,
        ND;
        cond = cond,
        dist_to_target = dist_to_target,
        incremental = true,
    )
    trust_tx[3] || return 2

    u .= .-u
    LA.ldiv!(x³, jacobian, u)
    iterative_refinement &&
    iterative_refinement!(x³, jacobian, u, norm; tol = min_acc, max_iters = max_iters)

    trust_tx[4] = taylor!(
        u,
        Val(4),
        H,
        tx³,
        t,
        AD,
        ND;
        dist_to_target = dist_to_target,
        cond = cond,
        incremental = true,
    )
    trust_tx[4] || return 3
    u .= .-u
    LA.ldiv!(x⁴, jacobian, u)

    return 4
end

function update_pade!(pred::Predictor, norm::AbstractNorm)
    tx = pred.tx⁴
    trust_tx = pred.trust_tx
    # This is an adaption of the algorithm outlined in
    # Robust Pad{\'e} approximation via SVD (Trefethen et al)
    τ = Inf
    λ_min = 3.814697265625e-6 # exp2(-18)
    for i = 1:length(tx)
        x, x¹, x², x³, x⁴ = tx[i]
        c¹ = abs(x¹)
        λ = max(λ_min, c¹)
        c¹ /= λ
        c² = abs(x²) / λ^2
        absx³ = abs(x³)
        c³ = absx³ / λ^3
        tol = 1e-14 * max(c¹, c², c³)
        if (c¹ ≤ tol && c² ≤ tol && c³ ≤ tol) || c² ≤ tol
            pred.taylor[i] = true
            pred.u[i] = x⁴
        else
            pred.taylor[i] = false
            τᵢ = (c² / c³) / λ
            τᵢ < τ && (τ = τᵢ)
            pred.u[i] = x⁴ - x³ * (x³ / x²)
        end
    end
    pred.local_error = norm(pred.u)
    pred.trust_region = τ

    pred
end


function predict!(x̂, pred::Predictor, H, x, t, Δt, J::Jacobian)
    if pred.pade
        predict_pade!(x̂, pred, H, Δt)
    else
        predict_ralston4!(x̂, pred, H, x, t, Δt, J)
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

function predict_ralston3!(x̂, pred::Predictor, H, x, t, Δt, J::Jacobian)
    @unpack u, y⁰, m₁, m₂, m₃ = pred
    # 0 | 0 0 0
    m₁ .= last(vectors(pred.tx¹))

    # 1/2 | 1/2 0 0
    @. y⁰ = x + 0.5 * Δt * m₁
    x̂ .= first.(y⁰)
    evaluate_and_jacobian!(u, matrix(J), H, x̂, t + 0.5Δt)
    taylor!(u, Val(1), H, y⁰, t + 0.5Δt, pred.AD, pred.ND; dist_to_target = abs(Δt))
    u .= .-u
    LA.ldiv!(m₂, updated!(J), u)

    # 3/4 | 0 3/4 0
    @. y⁰ = x + 0.75 * Δt * m₂
    x̂ .= first.(y⁰)
    evaluate_and_jacobian!(u, matrix(J), H, x̂, t + 0.75Δt)
    taylor!(u, Val(1), H, y⁰, t + 0.75Δt, pred.AD, pred.ND; dist_to_target = abs(Δt))
    u .= .-u
    LA.ldiv!(m₃, updated!(J), u)

    #   | 2/9 3/9 4/9
    @. x̂ = x + Δt * (2m₁ + 3m₂ + 4m₃) / 9
    x̂
end

function predict_ralston4!(x̂, pred::Predictor, H, x, t, Δt, J::Jacobian)
    @unpack u, y⁰, m₁, m₂, m₃ = pred
    # 0 | 0 0 0
    m₁ .= last(vectors(pred.tx¹))

    # 1/2 | 1/2 0 0
    @. y⁰ = x + 0.4 * Δt * m₁
    x̂ .= first.(y⁰)
    evaluate_and_jacobian!(u, matrix(J), H, x̂, t + 0.4Δt)
    taylor!(u, Val(1), H, y⁰, t + 0.4Δt, pred.AD, pred.ND; dist_to_target = abs(Δt))
    u .= .-u
    LA.ldiv!(m₂, updated!(J), u)
    # 1/1 | 0 1/2 0

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
    #   | 2/9 3/9 4/9
    @. x̂ = x + Δt * (0.17476028 * m₁ - 0.55148066 * m₂ + 1.20553560m₃ + 0.17118478m₄)
    x̂
end

order(P::Predictor) = 4
local_error(predictor::Predictor) = predictor.local_error
trust_region(predictor::Predictor) = predictor.trust_region
