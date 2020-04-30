# This file contains a Padé predictor of order (2,1) used for the tracking in tracker.jl
# For a motivation of this see:
# Mixed Precision Path Tracking for Polynomial Homotopy Continuation,
# Sascha Timme (2020), arXiv:1902.02968
mutable struct Predictor{N}
    AD::AD{N}
    ND::NumericalDifferentiation
    taylor::BitVector
    τ::Float64
    local_err::Vector{ComplexF64}
    ord::Int
    # buffer for the computation of derivatives
    tx⁰::TaylorVector{1,ComplexF64}
    tx¹::TaylorVector{2,ComplexF64}
    tx²::TaylorVector{3,ComplexF64}
    tx³::TaylorVector{4,ComplexF64}
    tx⁴::TaylorVector{5,ComplexF64}
    trust_tx::Vector{Bool}
    u::Vector{ComplexF64}
end

function Predictor(@nospecialize(H::AbstractHomotopy), ad::AD)
    m, n = size(H)
    ND = NumericalDifferentiation(n)
    local_err = zeros(ComplexF64, n)
    tx⁴ = TaylorVector{5}(ComplexF64, n)
    tx⁰ = TaylorVector{1}(tx⁴)
    tx¹ = TaylorVector{2}(tx⁴)
    tx² = TaylorVector{3}(tx⁴)
    tx³ = TaylorVector{4}(tx⁴)
    trust_tx = [true, true, true, true]
    u = zeros(ComplexF64, m)

    Predictor(ad, ND, falses(n), Inf, local_err, 4, tx⁰, tx¹, tx², tx³, tx⁴, trust_tx, u)
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
    compute_derivatives!(P, H, x, t, jacobian, norm; dist_to_target = dist_to_target)
    update_pade!(P)
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

    trust_tx[2] || return predictor
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
    trust_tx[3] || return predictor

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
    trust_tx[4] || return predictor
    u .= .-u
    LA.ldiv!(x⁴, jacobian, u)

    predictor
end

function update_pade!(pred::Predictor)
    tx = pred.tx⁴
    trust_tx = pred.trust_tx
    # This is an adaption of the algorithm outlined in
    # Robust Pad{\'e} approximation via SVD (Trefethen et al)
    τ = Inf
    λ_min = 3.814697265625e-6 # exp2(-18)
    if trust_tx[1] && trust_tx[2] && trust_tx[3]
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
                pred.local_err[i] = x⁴
            else
                pred.taylor[i] = false
                τᵢ = (c² / c³) / λ
                τᵢ < τ && (τ = τᵢ)
                pred.local_err[i] = x⁴ - x³ * (x³ / x²)
            end
        end
        pred.ord = 4
        pred.τ = τ
    else
        _, x¹, x², x³, = vectors(tx)
        pred.taylor .= true
        if trust_tx[2]
            pred.ord = 3
            pred.local_err .= x³
        else
            pred.ord = 2
            pred.local_err .= x²
        end
    end

    pred.τ = τ

    pred
end

function predict!(x̂, predictor::Predictor, H, Δt)
    tx = predictor.tx³
    for i = 1:length(tx)
        x, x¹, x², x³ = tx[i]
        if predictor.taylor[i]
            if predictor.ord ≥ 3
                x̂[i] = x + Δt * (x¹ + Δt * x²)
            else
                x̂[i] = x + Δt * x¹
            end
        else
            δᵢ = 1 - Δt * x³ / x²
            x̂[i] = x + Δt * (x¹ + Δt * x² / δᵢ)
        end
    end

    nothing
end

order(P::Predictor) = P.ord
local_error(predictor::Predictor) = predictor.local_err
trust_region(predictor::Predictor) = predictor.τ
