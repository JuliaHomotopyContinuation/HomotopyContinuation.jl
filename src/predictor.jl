#=

Predictor

## Introduction

In an ideal world the predictor only uses an order (2,1) Padé approximant to generate
an initial guess for the Newton corrector.
The choice of Padé approximants as predictors is motived by the results in [1].
In particular, the ability to provide a trust region of the prediction method to
significantly decrease the likelihood of path jumping.
additionally, sice Padé approximants are rational functions, these perform very well
for diverging paths and close to singularities.

Now, for an Padé approximant of order (2,1) the first 3 derivatives are needed.
The most robust way to achieve this it to use automatic differentiation.
However, due to significant compile times, we also have to consider the case
that the Padé approximants are computed by numerical differentiation using the scheme
derived in [2]. The numerical differentiation scheme only uses O(h^2) approximations
since we are mostly constrained by the available precision.

However, the smaller the trust region is (i.e. the closer we are to a singularity) the
smaller we need to take the discretization in the numerical differentiation.
Thus, the error in the numerical differentiation is inversely corrolated to the trust
region.
Additionally, the third derivative needs to be sufficiently accurate. Otherwise
the predictor doesn't generate sufficiently accurate intial guesses resulting
in very small step sizes.
Therefore it can necessary to drop down to an Padé approximant of order (1,1).
Usually the arithmetic is sufficient to obtain for the second derivative a sufficiently
accurate approximation.

If also the second derivative is not available we use cubic hermite interpolation.
For diverging paths this works much much better than (explicit) Runge-Kutta methods.
The reason that we do not use expl. RK methods is that, then when numerical differentiation
has to give up, we are in a regime where the ODE of the path is *stiff*. Thus, expl.
RK just don't yield satisfying results.

## Endgame and coordinates

Additionally, the predictor computes in s-coordinates, that is for the path x(t)
the predictor internally computes with y(s) = x(s^m) where m is the *winding number*.
If m = 1 then this is the same as computing with x(t). However, during the endgame
`m` is possibly set.
All APIs expect input in t-coordinates and return results in t-coordinates.


[1]: A Robust Numerical Path Tracking Algorithm for Polynomial Homotopy Continuation,
 Telen, Van Barel, Verschelde https://arxiv.org/abs/1909.04984

[2]:Mackens, Wolfgang. "Numerical differentiation of implicitly defined space curves."
   Computing 41.3 (1989): 237-260.
=#

module PredictionMethod
@enum methods begin
    Pade21
    Pade11
    Hermite
    Euler
end
end
using .PredictionMethod: PredictionMethod


## Type dispatch on automatic differentiation or numerical differentiation
struct AD{N} end
function AD(N::Int)
    if !(0 ≤ N ≤ 3)
        throw(ArgumentError("`automatic_differentiation` has to be between 0 and 3."))
    end
    AD{N}()
end


"""
    Predict(H::AbstractHomotopy, ::AD{N}) where {N}

Construct a predictor. The predicor uses a Padé approximant of order (2,1) to produce
an initial guess for the Newton corrector.
For this, the predictor tries to compute the local derivatives up to order 3.
This is done using automatic differentiation for the derivatives up to order `3`.
Otherwise numerical differentiation is used.
"""
Base.@kwdef mutable struct Predictor{T<:AD}
    AD::T

    method::PredictionMethod.methods = PredictionMethod.Pade21
    order::Int = 4
    use_hermite::Bool = true
    trust_region::Float64 = Inf
    local_error::Float64 = Inf
    cond_H_ẋ::Float64 = Inf

    tx⁰::TaylorVector{1,ComplexF64}
    tx¹::TaylorVector{2,ComplexF64}
    tx²::TaylorVector{3,ComplexF64}
    tx³::TaylorVector{4,ComplexF64}
    t::ComplexF64 = complex(NaN)
    tx_norm::Vector{Float64} = zeros(4)

    # finite diff
    xtemp::Vector{ComplexF64}
    u::Vector{ComplexF64}
    u₁::Vector{ComplexF64}
    u₂::Vector{ComplexF64}

    # for interpolation
    prev_tx¹::TaylorVector{2,ComplexF64}
    prev_t::ComplexF64 = complex(NaN)

    # endgame predictor
    winding_number::Int = 1
    s::ComplexF64 = complex(NaN)
    prev_s::ComplexF64 = complex(NaN)
    ty¹::TaylorVector{2,ComplexF64}
    prev_ty¹::TaylorVector{2,ComplexF64}
end


function Predictor(H::AbstractHomotopy, autodiff::AD)
    m, n = size(H)
    tx³ = TaylorVector{4}(ComplexF64, n)
    Predictor(
        AD = autodiff,
        tx⁰ = TaylorVector{1}(tx³),
        tx¹ = TaylorVector{2}(tx³),
        tx² = TaylorVector{3}(tx³),
        tx³ = tx³,
        xtemp = zeros(ComplexF64, n),
        u = zeros(ComplexF64, m),
        u₁ = zeros(ComplexF64, m),
        u₂ = zeros(ComplexF64, m),
        prev_tx¹ = TaylorVector{2}(ComplexF64, n),
        ty¹ = TaylorVector{2}(ComplexF64, n),
        prev_ty¹ = TaylorVector{2}(ComplexF64, n),
    )
end

function init!(predictor::Predictor)
    predictor.cond_H_ẋ = 1.0
    predictor.winding_number = 1

    predictor.t = predictor.prev_t = NaN
    predictor.s = predictor.prev_s = NaN
    predictor.trust_region = predictor.local_error = NaN
    predictor.use_hermite = true
    predictor
end

order(predictor::Predictor) = predictor.order
local_error(predictor::Predictor) = predictor.local_error
trust_region(predictor::Predictor) = predictor.trust_region

###################
### FINITE DIFF ###
###################

function g!(u, H, x::Vector, t, h, xtemp)
    evaluate!(u, H, x, t + h)
end
function g!(u, H, tx::TaylorVector{2}, t, h, xtemp)
    @inbounds for i = 1:length(tx)
        x, x¹ = tx[i]
        xtemp[i] = muladd(h, x¹, x)
    end
    evaluate!(u, H, xtemp, t + h)
end
function g!(u, H, tx::TaylorVector{3}, t, h, xtemp)
    @inbounds for i = 1:length(tx)
        x, x¹, x² = tx[i]
        xtemp[i] = muladd(h, muladd(h, x², x¹), x)
    end
    evaluate!(u, H, xtemp, t + h)
end

function g!(u, H, tx::TaylorVector{4}, t, h, xtemp)
    @inbounds for i = 1:length(tx)
        x, x¹, x², x³ = tx[i]
        xtemp[i] = muladd(h, muladd(h, muladd(h, x³, x²), x¹), x)
    end
    evaluate!(u, H, xtemp, t + h)
end

function finite_diff!(u, predictor::Predictor, H::AbstractHomotopy, x, t, h; order::Int)
    @unpack u₁, u₂, xtemp = predictor

    g!(u₁, H, x, t, h, xtemp)
    g!(u₂, H, x, t, -h, xtemp)

    hk = h^order
    if iseven(order)
        u .= 0.5 .* (u₁ .+ u₂) ./ hk
    else
        u .= 0.5 .* (u₁ .- u₂) ./ hk
    end
end


function finite_diff_taylor!(
    u,
    ::Val{N},
    predictor::Predictor,
    H::AbstractHomotopy,
    x,
    t;
    prev_λ::Float64 = NaN,
) where {N}
    # Use not machine precision to account for some error in the evaluation
    ε = 1.4210854715202004e-14
    # λ should be a scaling factor such that the derivatives of x(t/λ) are of of the same
    # order of magnitude. We can either compute this by balancing the derivatives or
    # re-using the trust region size of the previous step
    if N == 1
        λ = isfinite(prev_λ) ? prev_λ : 1.0
    elseif N == 2
        λ = isfinite(prev_λ) ? prev_λ : predictor.tx_norm[1] / predictor.tx_norm[2]
    elseif N == 3
        λ = isfinite(prev_λ) ? prev_λ : predictor.tx_norm[2] / predictor.tx_norm[3]
    end
    # truncation err = (1/λ)^2*h^2
    # round-off err   = ε/h^N
    # -->
    # Compute optimal ĥ by
    #      round-off err = trunc err
    # <=>  ελ^-N/ĥ^3 = λ^(-N-2)*ĥ^2
    # <=>  (ε λ^2)^(1/(N+2)) = ĥ
    # To be a little bit pessimistic, reduce λ by a quarter
    λ *= 0.25

    h = (ε * λ^2)^(1 / (N + 2))
    # h should be at most half λ
    h = min(h, 0.5λ)
    # Check that truncation and round-off error are both acceptable
    trunc_err = h^2 / λ^2
    rnd_err = ε / h^N
    err = trunc_err + rnd_err
    # @show N
    # @show λ, h, trunc_err, rnd_err
    if trunc_err + rnd_err > 0.01 || !isfinite(h)
        return false, err
    end

    finite_diff!(u, predictor, H, x, t, h; order = N)

    return true, err^2
end

## Default handling ignores incremental
ModelKit.taylor!(u, v::Val, H::AbstractHomotopy, tx, t, incremental::Bool) =
    taylor!(u, v, H, tx, t)
@generated function ModelKit.taylor!(
    u,
    v::Val{M},
    predictor::Predictor{AD{N}},
    H,
    tx,
    t,
) where {M,N}
    if M ≤ N
        quote
            taylor!(u, v, H, tx, t, true)
            true, eps()
        end
    else
        quote
            finite_diff_taylor!(u, v, predictor, H, tx, t; prev_λ = predictor.trust_region)
        end
    end
end

function winding_number!(P::Predictor, m)
    P.winding_number = m
    P
end

"""
    update!(
        P::Predictor,
        H::AbstractHomotopy,
        x,
        t,
        jacobian::Jacobian,
        norm::AbstractNorm,
        x̂ = nothing,
    )

Upadte the predictor with the new solution `(x,t)` of `H(x,t) = 0`.
This computes new local derivatives and chooses an appriopriate prediction method.
"""
function update!(
    predictor::Predictor,
    H::AbstractHomotopy,
    x,
    t,
    J::Jacobian,
    norm::AbstractNorm,
    x̂ = nothing, # the predicted value
)
    @unpack u, tx⁰, tx¹, tx², tx³, ty¹, xtemp, tx_norm = predictor
    x⁰, x¹, x², x³ = vectors(tx³)
    y⁰, y¹ = vectors(ty¹)
    m = predictor.winding_number

    @inbounds for i = 1:length(xtemp)
        predictor.prev_tx¹[i, 1] = predictor.tx¹[i, 1]
        predictor.prev_tx¹[i, 2] = predictor.tx¹[i, 2]
    end

    predictor.prev_t, predictor.t = predictor.t, t

    if m > 1
        predictor.prev_s = predictor.s
        predictor.s = t_to_s_plane(t, m)
    end

    if isnothing(x̂)
        predictor.local_error = NaN
    else
        Δs = fast_abs(t - predictor.prev_t)
        predictor.local_error = norm(x̂, x) / Δs^predictor.order
    end

    x⁰ .= x
    tx_norm[1] = norm(x)
    if m > 1
        y⁰ .= x
    end

    # Compute ẋ
    trust, err = taylor!(u, Val(1), predictor, H, x, t)
    u .= .-u
    LA.ldiv!(xtemp, J, u)
    # Check error made in the linear algebra
    δ = fixed_precision_iterative_refinement!(xtemp, workspace(J), u, norm)
    predictor.cond_H_ẋ = δ / eps()
    tol_δ₁ = max(1e-10, err)
    δ > tol_δ₁ && iterative_refinement!(xtemp, J, u, InfNorm(); tol = tol_δ₁, max_iters = 5)
    prev_tx_norm2 = tx_norm[2]
    tx_norm[2] = norm(xtemp)
    x¹ .= xtemp

    if m > 1
        μ = m * predictor.s^(m - 1)
        y¹ .= μ .* xtemp
        @goto use_hermite
    end

    trust, err = taylor!(u, Val(2), predictor, H, tx¹, t)

    if !trust
        # Use cubic hermite
        @label use_hermite
        if isfinite(predictor.prev_t) && predictor.use_hermite
            predictor.method = PredictionMethod.Hermite
            predictor.order = 3
            predictor.trust_region = tx_norm[1] / tx_norm[2]
            if isnan(predictor.local_error)
                predictor.local_error = (tx_norm[2] / tx_norm[1])^3
            end
            return predictor
        else
            # don't use the previous local error estimate if we downgrade to Euler
            if predictor.method != PredictionMethod.Euler || isnan(predictor.local_error)
                predictor.local_error = (tx_norm[2] / tx_norm[1])^2
            end
            predictor.method = PredictionMethod.Euler
            predictor.order = 2
            predictor.trust_region = tx_norm[1] / tx_norm[2]
            return predictor
        end
    end
    u .= .-u
    LA.ldiv!(xtemp, J, u)
    tol_δ₂ = max(1e-8, err)
    δ > tol_δ₂ && iterative_refinement!(xtemp, J, u, norm; tol = tol_δ₂, max_iters = 4)
    tx_norm[3] = norm(xtemp)
    x² .= xtemp

    trust, err = taylor!(u, Val(3), predictor, H, tx², t)
    if !trust
        @label use_pade11
        # only trust tx_norm[2] / tx_norm[3] if tx_norm[1] / tx_norm[2] is of the same order
        if 0.1 ≤ (tx_norm[1] / tx_norm[2]) / (tx_norm[2] / tx_norm[3]) ≤ 10
            # Use Padé (1,1) predictor
            predictor.method = PredictionMethod.Pade11
            predictor.order = 3
            predictor.trust_region = tx_norm[2] / tx_norm[3]
            if isnan(predictor.local_error)
                predictor.local_error = (tx_norm[3] / tx_norm[2])^3
            end
            return predictor
        else
            @goto use_hermite
        end
    end

    u .= .-u
    LA.ldiv!(xtemp, J, u)
    tol_δ₃ = max(1e-6, err)
    δ > tol_δ₃ && iterative_refinement!(xtemp, J, u, norm; tol = tol_δ₃, max_iters = 3)
    tx_norm[4] = norm(xtemp)
    x³ .= xtemp

    if 0.1 ≤ (tx_norm[2] / tx_norm[3]) / (tx_norm[3] / tx_norm[4]) ≤ 10
        # Use Padé (2,1) predictor
        predictor.method = PredictionMethod.Pade21
        predictor.order = 4
        predictor.trust_region = tx_norm[3] / tx_norm[4]
        if isnan(predictor.local_error)
            predictor.local_error = (tx_norm[4] / tx_norm[3])^4
        end
    else
        @goto use_pade11
    end

    predictor
end

predict!(x̂, pred::Predictor, H, x, t, Δt) = predict!(x̂, pred, pred.method, H, x, t, Δt)
function predict!(x̂, pred::Predictor, method::PredictionMethod.methods, H, x, t, Δt)
    m = pred.winding_number

    if method == PredictionMethod.Pade21
        λ = pred.trust_region
        λ² = λ^2
        λ³ = λ^3
        tol = 1e-8
        for (i, (x, x¹, x², x³)) in enumerate(pred.tx³)
            c, c¹, c², c³ = fast_abs.((x, x¹, x², x³))
            # check if only taylor series is used
            τ = tol * √(c^2 + (c¹ * λ)^2 + (c² * λ²)^2 + (c³ * λ³)^2)
            if c³ * λ³ ≤ τ || c² * λ² ≤ τ
                x̂[i] = x + Δt * (x¹ + Δt * x²)
            else
                δᵢ = 1 - Δt * x³ / x²
                x̂[i] = x + Δt * (x¹ + Δt * x² / δᵢ)
            end
        end
    elseif method == PredictionMethod.Pade11
        λ = pred.trust_region
        λ² = λ^2
        tol = 1e-8
        for (i, (x, x¹, x²)) in enumerate(pred.tx²)
            c, c¹, c² = fast_abs.((x, x¹, x²))
            # check if only taylor series is used
            τ = tol * √(c^2 + (c¹ * λ)^2 + (c² * λ²)^2)
            if c¹ * λ ≤ τ || c² * λ² ≤ τ
                x̂[i] = x + Δt * x¹
            else
                δᵢ = 1 - Δt * x² / x¹
                x̂[i] = x + Δt * x¹ / δᵢ
            end
        end
    elseif method == PredictionMethod.Hermite
        if m > 1
            prev_s = t_to_s_plane(pred.prev_t, m)
            s = t_to_s_plane(t, m)
            s′ = t_to_s_plane(t + Δt, m)
            prev_sm = m * prev_s^(m - 1)
            for i = 1:length(pred.prev_tx¹)
                pred.prev_ty¹[i, 1] = pred.prev_tx¹[i, 1]
                pred.prev_ty¹[i, 2] = prev_sm * pred.prev_tx¹[i, 2]
            end
            sm = m * s^(m - 1)
            for i = 1:length(pred.prev_tx¹)
                pred.ty¹[i, 1] = pred.tx¹[i, 1]
                pred.ty¹[i, 2] = sm * pred.tx¹[i, 2]
            end
            cubic_hermite!(x̂, pred.prev_ty¹, prev_s, pred.ty¹, s, s′)
        else
            cubic_hermite!(x̂, pred.prev_tx¹, pred.prev_t, pred.tx¹, t, t + Δt)

            # ŷ = map(1:length(pred.tx¹)) do i
            #     t₀ = pred.prev_t# / t
            #     t₁ = t #1
            #     x₀, _ = pred.prev_tx¹[i]
            #     x₁, ẋ₁ = pred.tx¹[i]
            #     #ẋ₁ = t
            #     a₀, a₁, b₁ = Complex{BigFloat}[1 t₀ -x₀*t₀; 1 t₁ -x₁*t₁; 0 1  -x₁-ẋ₁*t₁] \ [x₀, x₁, ẋ₁]
            #     # a₀, a₁, b₁ = [1 t₀ -x₀ * t₀; 1 0 0; 0 1 -x₁] \ [x₀, x₁, ẋ₁]
            #     # @show a₀, a₁, b₁
            #
            #     # @show [1 t₀ -x₀*t₀; 1 t₁ -x₁*t₁; 0 1  -x₁-ẋ₁*t₁]
            #     (a₀ + ((t + Δt)) * a₁) / (1 + ((t + Δt)) * b₁)
            # end
            # @show InfNorm()(x̂, ŷ)
            # x̂ .= ŷ
        end
    elseif method == PredictionMethod.Euler
        if m > 1
            s = t_to_s_plane(t, m)
            s′ = t_to_s_plane(t + Δt, m)
            Δs = s′ - s
            sm = m * s^(m - 1)
            for (i, (x, x¹)) in enumerate(pred.tx¹)
                x̂[i] = x + Δs * sm * x¹
            end
        else
            for (i, (x, x¹)) in enumerate(pred.tx¹)
                x̂[i] = x + Δt * x¹
            end
        end
    end

    x̂
end


"""
    t_to_s_plane(t::Complex, m::Int; branch::Int = 0)

Convert ``t = r * \\exp(√(-1)θ)`` to ``t = r^{1/m} * \\exp(√(-1)θ/m)`` where for
``r^{1/m}`` the positive real root is choosen.
The `branch` keyword allows to choose a different branch cut.
"""
function t_to_s_plane(t::Complex, m::Int; branch::Int = 0)
    r = fast_abs(t)
    if isreal(t) && real(t) > 0
        if branch > 0
            complex(nthroot(r, m), 2branch * π / m)
        else
            complex(nthroot(r, m))
        end
    else
        # angle returns a value in [-π,π], but we want [0,π]
        θ = mod2pi(angle(t)) + branch * 2π
        return nthroot(r, m) * cis(θ / m)
    end
end


function cubic_hermite!(x̂, tx¹₀, t₀, tx¹₁, t₁, t)
    if isreal(t₀) && isreal(t₁) && isreal(t)
        _cubic_hermite!(x̂, tx¹₀, real(t₀), tx¹₁, real(t₁), real(t))
    else
        _cubic_hermite!(x̂, tx¹₀, t₀, tx¹₁, t₁, t)
    end
end
@inline function _cubic_hermite!(x̂, tx¹₀, t₀, tx¹₁, t₁, t)
    s = (t - t₀) / (t₁ - t₀)
    h₀₀ = (1 + 2s) * (1 - s)^2
    h₁₀ = (t - t₀) * (1 - s)^2
    h₀₁ = s^2 * (3 - 2s)
    h₁₁ = (t - t₀) * s * (s - 1)
    @inbounds for i in eachindex(x̂)
        x̂[i] = h₀₀ * tx¹₀[i, 1] + h₁₀ * tx¹₀[i, 2] + h₀₁ * tx¹₁[i, 1] + h₁₁ * tx¹₁[i, 2]
    end
    x̂
end
