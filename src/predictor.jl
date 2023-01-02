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
    Pade31
    Pade21
    Hermite
end
end

using .PredictionMethod: PredictionMethod

"""
    Predict(H::AbstractHomotopy)

Construct a predictor. The predicor uses a Padé approximant of order (2,1) to produce
an initial guess for the Newton corrector.
For this, the predictor tries to compute the local derivatives up to order 3.
This is done using automatic differentiation for the derivatives up to order `3`.
Otherwise numerical differentiation is used.
"""
Base.@kwdef mutable struct Predictor
    method::PredictionMethod.methods = PredictionMethod.Pade31
    order::Int = 5
    trust_region::Float64 = Inf
    local_error::Float64 = Inf
    cond_H_ẋ::Float64 = Inf

    tx⁰::TaylorVector{1,ComplexF64}
    tx¹::TaylorVector{2,ComplexF64}
    tx²::TaylorVector{3,ComplexF64}
    tx³::TaylorVector{4,ComplexF64}
    tx⁴::TaylorVector{5,ComplexF64}
    t::ComplexF64 = complex(NaN)
    tx_norm::Vector{Float64} = zeros(5)

    # finite diff
    xtemp::Vector{ComplexF64}
    u::Vector{ComplexF64}

    # # endgame predictor
    # winding_number::Int = 1
    # s::ComplexF64 = complex(NaN)
    # prev_s::ComplexF64 = complex(NaN)
    # ty¹::TaylorVector{2,ComplexF64}
    # prev_ty¹::TaylorVector{2,ComplexF64}
end


function Predictor(H::AbstractHomotopy)
    m, n = size(H)
    tx⁴ = TaylorVector{5}(ComplexF64, n)
    Predictor(
        tx⁰ = TaylorVector{1}(tx⁴),
        tx¹ = TaylorVector{2}(tx⁴),
        tx² = TaylorVector{3}(tx⁴),
        tx³ = TaylorVector{4}(tx⁴),
        tx⁴ = tx⁴,
        xtemp = zeros(ComplexF64, n),
        u = zeros(ComplexF64, m),
        # u₁ = zeros(ComplexF64, m),
        # u₂ = zeros(ComplexF64, m),
        # prev_tx¹ = TaylorVector{2}(ComplexF64, n),
        # ty¹ = TaylorVector{2}(ComplexF64, n),
        # prev_ty¹ = TaylorVector{2}(ComplexF64, n),
    )
end

function init!(predictor::Predictor)
    predictor.cond_H_ẋ = 1.0
    # predictor.winding_number = 1

    predictor.t = NaN
    # predictor.s = predictor.prev_s = NaN
    predictor.trust_region = NaN
    predictor.local_error = NaN
    # predictor.use_hermite = true
    predictor
end

order(predictor::Predictor) = predictor.order
local_error(predictor::Predictor) = predictor.local_error
trust_region(predictor::Predictor) = predictor.trust_region

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
    x̂ = nothing;
    winding_number::Int = 1,
)
    # The general strategy is as follows:
    #
    # If m = winding_number > 1 then we only use a hermite predictor.
    #
    # Otherwise we use a Padé predictor of order (2,1). For this we need the first three
    # derivativatives of the path x(t).
    # We obtain accurate derivatives by using automatic differentiation (AD) with the taylor!
    # api.

    @unpack u, tx⁰, tx¹, tx², tx³, tx⁴, xtemp, tx_norm = predictor

    x⁰, x¹, x², x³, x⁴ = vectors(tx⁴)

    x⁰ .= x
    tx_norm[1] = norm(x)

    # Compute first derivative
    taylor!(u, Val(1), H, x, t)
    u .= .-u
    LA.ldiv!(xtemp, J, u, norm)

    # establish linear algebra error using fixed precision error
    # The error of the fixed precision refinement is `2n * cond(J, ẋ) * eps()`
    δ, fixed_refinement_iters = iterative_refinement!(
        xtemp,
        J,
        u,
        norm;
        mixed_precision = false,
        tol = 1e-6,
        max_iters = 3,
    )
    cond_H_ẋ = δ / (2 * length(x) * eps())
    predictor.cond_H_ẋ = cond_H_ẋ

    # target 
    refine =
        () -> begin
            (fixed_refinement_iters == 1) && return
            iterative_refinement!(
                xtemp,
                J,
                u,
                norm;
                mixed_precision = cond_H_ẋ > 1e6,
                tol = 1e-6,
                max_iters = 3,
            )
        end
    if cond_H_ẋ > 1e6
        refine()
    end
    tx_norm[3] = norm(xtemp)
    x¹ .= xtemp


    taylor!(u, Val(2), H, tx¹, t, true)
    u .= .-u
    LA.ldiv!(xtemp, J, u, norm)
    refine()
    tx_norm[3] = norm(xtemp)
    x² .= xtemp


    taylor!(u, Val(3), H, tx², t, true)
    u .= .-u
    LA.ldiv!(xtemp, J, u, norm)
    refine()
    tx_norm[4] = norm(xtemp)
    x³ .= xtemp


    taylor!(u, Val(4), H, tx³, t, true)
    u .= .-u
    LA.ldiv!(xtemp, J, u, norm)
    refine()
    tx_norm[5] = norm(xtemp)
    x⁴ .= xtemp


    # Now detect whether we can actually use a taylor predictor of order 4

    c, c¹, c², c³, c⁴ = tx_norm
    λ = max(c, 1)
    c¹ /= λ
    c² /= λ^2
    c³ /= λ^3
    c⁴ /= λ^4
    tol = 1e-14 * c

    # check if all derivatives are zero effectively
    if c⁴ > tol
        predictor.order = 5
        τ = c³ / (λ * c⁴)
    elseif c³ > tol
        predictor.order = 4
        τ = c² / (λ * c³)
    elseif c² > tol
        predictor.order = 3
        τ = c¹ / (λ * c²)
    elseif c¹ > tol
        predictor.order = 2
        τ = c / (λ * c¹)
    else
        predictor.order = 1
        τ = 1.0
    end
    predictor.trust_region = τ

    predictor
end

function predict!(x̂, pred::Predictor, H, x, t, Δt; winding_number::Int = 1)
    @unpack u, tx⁰, tx¹, tx², tx³, tx⁴, xtemp, tx_norm = pred

    τ = pred.trust_region
    if (pred.order <= 2)
        for (i, (xᵢ, xᵢ¹)) in enumerate(tx²)
            x̂[i] = xᵢ + Δt * xᵢ¹
        end
    elseif (pred.order == 3)
        λ = τ
        for (i, (xᵢ, xᵢ¹, xᵢ²)) in enumerate(tx²)
            δᵢ = 1 - Δt * Base.FastMath.div_fast(λ * xᵢ², λ * xᵢ¹)
            x̂[i] = @fastmath xᵢ + Δt * xᵢ¹ / δᵢ
        end
        τ = norm(x¹) / norm(x²)
        return τ
    elseif (pred.order == 4)
        λ = τ^2
        for (i, (xᵢ, xᵢ¹, xᵢ², xᵢ³)) in enumerate(tx³)
            δᵢ = 1 - Δt * Base.FastMath.div_fast(λ * xᵢ³, λ * xᵢ²)
            x̂[i] = @fastmath xᵢ + Δt * (xᵢ¹ + Δt * xᵢ² / δᵢ)
        end
    elseif (pred.order == 5)
        λ = τ^3
        for (i, (xᵢ, xᵢ¹, xᵢ², xᵢ³, xᵢ⁴)) in enumerate(tx⁴)
            δᵢ = 1 - Δt * Base.FastMath.div_fast(λ * xᵢ⁴, λ * xᵢ³)
            x̂[i] = @fastmath xᵢ + Δt * (xᵢ¹ + Δt * (xᵢ² + Δt * xᵢ³ / δᵢ))
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


