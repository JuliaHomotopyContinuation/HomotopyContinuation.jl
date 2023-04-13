Base.@kwdef struct Predictor
    max_order::Int = 5
    order::RefValue{Int} = Ref(5)
    trust_region::RefValue{Float64} = Ref(Inf)
    t::RefValue{ComplexF64} = Ref(complex(NaN))

    ws::LinearSolveWorkspace
    norm::WeightedNorm{InfNorm}
    tx⁰::TaylorVector{1,ComplexF64}
    tx¹::TaylorVector{2,ComplexF64}
    tx²::TaylorVector{3,ComplexF64}
    tx³::TaylorVector{4,ComplexF64}
    tx⁴::TaylorVector{5,ComplexF64}
    scaling_factor::RefValue{Float64} = Ref(1.0)
    tx_norm::Vector{Float64} = zeros(5)

    x̂::Vector{ComplexF64}
    work_m::Vector{ComplexF64}
    work_n::Vector{ComplexF64}
end


function Predictor(
    H::AbstractHomotopy,
    ws::LinearSolveWorkspace,
    norm::WeightedNorm{InfNorm};
    max_order = 5,
)
    m, n = size(H)
    tx⁴ = TaylorVector{5}(ComplexF64, n)
    Predictor(
        norm = norm,
        max_order = max_order,
        tx⁰ = TaylorVector{1}(tx⁴),
        tx¹ = TaylorVector{2}(tx⁴),
        tx² = TaylorVector{3}(tx⁴),
        tx³ = TaylorVector{4}(tx⁴),
        tx⁴ = tx⁴,
        ws = ws,
        x̂ = zeros(ComplexF64, n),
        work_m = zeros(ComplexF64, m),
        work_n = zeros(ComplexF64, n),
    )
end

Base.show(io::IO, p::Predictor) =
    print(io, "Predictor(order=$(p.order[]), trust_region=$(p.trust_region[]))")

function init!(predictor::Predictor)
    predictor.t[] = NaN
    predictor.trust_region[] = NaN
    predictor.scaling_factor[] = 1.0
    predictor
end

order(predictor::Predictor) = 5
trust_region(predictor::Predictor) = predictor.trust_region[]

function assign_ẋ!(y, predictor::Predictor)
    _, ẋ = vectors(predictor.tx¹)
    α = predictor.scaling_factor[]
    y .= ẋ ./ α
end


function localize!(predictor::Predictor, H::AbstractHomotopy, x, t;)
    @unpack tx⁰, tx¹, tx², tx³, tx⁴, work_m, work_n, tx_norm, ws, max_order, norm =
        predictor
    x⁰, x¹, x², x³, x⁴ = vectors(tx⁴)



    x⁰ .= x
    tx_norm[1] = norm(x)

    predictor.scaling_factor[] = predictor.trust_region[]
    if isnan(predictor.scaling_factor[]) || isinf(predictor.scaling_factor[])
        predictor.scaling_factor[] = 1.0
    end
    α = predictor.scaling_factor[]

    # Compute first derivative
    taylor!(work_m, Val(1), H, x, (t, α), true)
    set_b!(ws, work_m)

    solve!(ws)

    # We do fixed refinement with just a single step to establish the forward error
    # Since we scale our linear system, we usually don't see that fixed precision refinement
    # works (the linear system is basically solved as good as it gets).
    # So our only option then is to do extended_precision refinement.
    # Here we assume quadratic convergence, i.e., we require target_forward_error for the refinment
    # of  √(target_forward_error / 100) where the factor 100 is for some safety margin
    # (basically avoid an expensive(ish) extended_prec evaluation just to establish an error bound)

    target_forward_error = 1e-8

    fixed_prec_res =
        refine!(ws; max_iters = 1, tol = target_forward_error, extended_precision = false)
    refine_needed = fixed_prec_res.initial_ferr > target_forward_error
    extended_prec = fixed_prec_res.ferr > target_forward_error
    if extended_prec
        refine_tol = sqrt(target_forward_error) / 10
    else
        refine_tol = target_forward_error
    end

    refine_needed &&
        extended_prec &&
        refine!(ws; max_iters = 2, tol = refine_tol, extended_precision = extended_prec)


    tx_norm[2] = norm(solution(ws))
    x¹ .= .-solution(ws)

    if max_order >= 3
        taylor!(work_m, Val(2), H, tx¹, (t, α), true)

        set_b!(ws, work_m)
        solve!(ws)
        refine_needed &&
            refine!(ws; max_iters = 2, tol = refine_tol, extended_precision = extended_prec)

        tx_norm[3] = norm(solution(ws))
        x² .= .-solution(ws)
    end

    if max_order >= 4
        taylor!(work_m, Val(3), H, tx², (t, α), true)
        set_b!(ws, work_m)
        solve!(ws)
        refine_needed &&
            refine!(ws; max_iters = 2, tol = refine_tol, extended_precision = extended_prec)

        tx_norm[4] = norm(solution(ws))
        x³ .= .-solution(ws)
    end


    if max_order >= 5
        taylor!(work_m, Val(4), H, tx³, (t, α), true)
        set_b!(ws, work_m)
        solve!(ws)
        refine_needed &&
            refine!(ws; max_iters = 2, tol = refine_tol, extended_precision = extended_prec)

        tx_norm[5] = norm(solution(ws))
        x⁴ .= .-solution(ws)
    end


    # Now detect whether we can actually use a taylor predictor of order 4

    τ = Inf
    for (i, (x, x¹, x², x³, x⁴)) in enumerate(tx⁴)
        _c, c¹, c², c³, c⁴ = fast_abs.((x, x¹, x², x³, x⁴))
        tol = 1e-14 * max(c¹, c², c³, c⁴)
        if !((c¹ ≤ tol && c² ≤ tol && c³ ≤ tol && c⁴ ≤ tol) || c³ ≤ tol)
            if (max_order == 5)
                τᵢ = (c³ / c⁴)
            elseif (max_order == 4)
                τᵢ = (c² / c³)
            elseif (max_order == 3)
                τᵢ = (c¹ / c²)
            end

            if τᵢ < τ
                τ = τᵢ
            end
        end
    end
    predictor.order[] = max_order
    if !isfinite(τ) && max_order == 5
        τ = tx_norm[4] / tx_norm[5]
    end
    if !isfinite(τ) && max_order == 4
        predictor.order[] = 4
        τ = tx_norm[3] / tx_norm[4]
    end
    if !isfinite(τ) && max_order == 3
        predictor.order[] = 3
        τ = tx_norm[2] / tx_norm[3]
    end
    if !isfinite(τ)
        predictor.order[] = 2
        τ = tx_norm[1] / tx_norm[2]
    end
    if !isfinite(τ)
        predictor.order[] = 1
    end

    predictor.trust_region[] = τ * α

    predictor
end

function predict!(pred::Predictor, Δs)
    @unpack x̂, tx⁰, tx¹, tx², tx³, tx⁴ = pred

    α = pred.scaling_factor[]
    Δt = Δs / α
    order = pred.order[]
    if (order <= 2)
        @inbounds for (i, (xᵢ, xᵢ¹)) in enumerate(tx²)
            x̂[i] = xᵢ + Δt * xᵢ¹
        end
    elseif (order == 3)
        @inbounds for (i, (xᵢ, xᵢ¹, xᵢ²)) in enumerate(tx²)
            δᵢ = 1 - Δt * xᵢ² / xᵢ¹
            x̂[i] = xᵢ + Δt * xᵢ¹ / ifelse(isnan(δᵢ), one(δᵢ))
        end
    elseif (order == 4)
        @inbounds for (i, (xᵢ, xᵢ¹, xᵢ², xᵢ³)) in enumerate(tx³)
            δᵢ = 1 - Δt * xᵢ³ / xᵢ²
            x̂[i] = xᵢ + Δt * (xᵢ¹ + Δt * xᵢ² / ifelse(isnan(δᵢ), one(δᵢ), δᵢ))
        end
    elseif (order == 5)
        @inbounds for (i, (xᵢ, xᵢ¹, xᵢ², xᵢ³, xᵢ⁴)) in enumerate(tx⁴)
            δᵢ = 1 - Δt * xᵢ⁴ / xᵢ³
            x̂[i] = xᵢ + Δt * (xᵢ¹ + Δt * (xᵢ² + Δt * xᵢ³ / ifelse(isnan(δᵢ), one(δᵢ), δᵢ)))
        end
    end
end

prediction(predictor::Predictor) = predictor.x̂
