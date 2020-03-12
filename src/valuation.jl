"""
    Valuation(x::PVector; affine = true)
    Valuation(x::AbstractVector)

A data structure for approximating the valuation of a solution path ``x(s)``.
It does this by computing an approximation function [`ν`](@ref) and approximating
the first and second derivatives ``ν̇`` and ``ν̈``.
If `affine = true` then the computed valuation is the valuation pulled back
to the affine space.
The valuation is estimated continously along a solution path. For this it is assumed that
the path is tracked in a **logarithmic time scale**.
"""
mutable struct Valuation
    val_x::Vector{Float64}
    val_tx¹::Vector{Float64}
    Δval_x::Vector{Float64}
    Δval_tx¹::Vector{Float64}
    δ_x::Vector{Float64}
end

function Valuation(n::Integer)
    Valuation((zeros(n) for i = 1:5)...)
end
Valuation(x::AbstractVector) = Valuation(length(x))

function Base.show(io::IO, val::Valuation)
    print(io, "Valuation :")
    for field in [:val_x, :val_tx¹, :Δval_x, :Δval_tx¹, :δ_x]
        vs = [Printf.@sprintf("%#.4g", v) for v in getfield(val, field)]
        print(io, "\n • ", field, " → ", "[", join(vs, ", "), "]")
    end
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", v::Valuation) = v

"""
    init!(val::Valuation)

Initialize the valuation estimator to start from scratch.
"""
function init!(val::Valuation)
    val.val_x .= 0
    val.val_tx¹ .= 0
    val.Δval_x .= 0
    val.Δval_tx¹ .= 0
    val.δ_x .= 0
    val
end

function func_ν(u::Real, u¹::Real, t)
    x² = u^2
    μ = u * u¹
    t * μ / x²
end

function ν(x::Number, ẋ::Number, t)
    u, v = reim(x)
    u¹, v¹ = reim(ẋ)
    x² = abs2(x)
    μ = u * u¹ + v * v¹
    t * μ / x²
end

function ν_ν¹(x, x¹, x², t)
    u, v = reim(x)
    u¹, v¹ = reim(x¹)
    u², v² = reim(x²)

    xx = abs2(x)

    μ = u * u¹ + v * v¹
    l = μ / xx

    μ¹ = u * u² + u¹^2 + v * v² + v¹^2
    l¹ = μ¹ / xx - 2 * l^2

    t * l, t * l¹ + l
end

function ν_ν¹_ν²(x, x¹, x², x³, t)
    u, v = reim(x)
    u¹, v¹ = reim(x¹)
    u², v² = reim(x²)
    u³, v³ = reim(x³)

    xx = abs2(x)

    μ = u * u¹ + v * v¹
    l = μ / xx

    μ¹ = u * u² + u¹^2 + v * v² + v¹^2
    l¹ = μ¹ / xx - 2 * l^2

    μ² = u * u³ + 3u¹ * u² + v * v³ + 3v¹ * v²
    l² = μ² / xx - 2 * μ * μ¹ / xx^2 - 4l * l¹

    t * l, t * l¹ + l, t * l² + 2l¹
end


function update!(
    val::Valuation,
    x::AbstractVector,
    x¹::AbstractVector,
    x²::AbstractVector,
    x³::AbstractVector,
    t::Real,
)
    for i in eachindex(x)
        xᵢ, x¹ᵢ, x²ᵢ, x³ᵢ = x[i], x¹[i], 2 * x²[i], 6 * x³[i]
        ν, ν¹, ν² = ν_ν¹_ν²(xᵢ, x¹ᵢ, x²ᵢ, x³ᵢ, t)
        # @show ν, ν¹, ν²
        ν_t = t * ν¹
        δ_x = t * (t * ν² + ν¹) / ν_t
        Δx = ν_t / δ_x
        val_x = ν
        if !isnan(Δx) && δ_x > 0
            val_x -= Δx
        end

        val.val_x[i] = val_x
        val.Δval_x[i] = Δx
        val.δ_x[i] = δ_x

        # valuation of tx¹
        ν, ν¹ = ν_ν¹(x¹ᵢ, x²ᵢ, x³ᵢ, t)
        Δtx¹ = t * ν¹ / δ_x
        val_tx¹ = ν + 1
        if !isnan(Δtx¹) && δ_x > 0
            val_tx¹ -= Δtx¹
        end

        val.val_tx¹[i] = val_tx¹
        val.Δval_tx¹[i] = Δtx¹
    end

    val
end


function mean_dev(a::Vararg{T,N}) where {T,N}
    m = sum(a) / N
    σ = sqrt(sum((a .- m) .^ 2)) / abs(m)
    m, σ
end

"""
    analyze(val::Valuation; tol, tol_at_infinity)::ValuationVerdict

Judge the current valuation to determine whether a path is diverging. Returns
a [`ValuationVerdict`](@ref).
A path is diverging if one entry of the valuation is negative. To assure that
the estimate is trust worthy we require that the first and second derivative
of [`ν``](@ref) is smaller than `tol_at_infinity`.
The `tol` is used to estimate whether a finite valuation is trustworthy.
"""
function analyze(
    val::Valuation;
    finite_tol::Float64 = throw(UndefKeywordError(:finite_tol)),
    at_infinity_tol::Float64 = throw(UndefKeywordError(:at_infinity_tol)),
    strict_at_infinity_tol::Float64 = at_infinity_tol^2,
    singular_tol::Float64 = throw(UndefKeywordError(:singular_tol)),
    max_winding_number::Int = throw(UndefKeywordError(:max_winding_number)),
    zero_is_finite::Bool = throw(UndefKeywordError(:zero_is_finite)),
)
    at_zero = false
    at_infinity = false
    singular = false
    indecisive = false

    @unpack val_x, val_tx¹, Δval_x, Δval_tx¹, δ_x = val

    ε = 1 / max_winding_number

    at_infinity = false
    strict_at_infinity = false
    at_zero = false
    strict_at_zero = false
    finite = true
    for (i, val_xᵢ) in enumerate(val_x)
        m = 0.5 * (val_xᵢ + val_tx¹[i])
        σ = max(abs(1 - val_xᵢ / m), abs(1 - val_tx¹[i] / m))
        Δ = max(abs(Δval_tx¹[i]), σ)

        # ∞ check
        if m + 2 * Δ < -finite_tol && δ_x[i] > 0
            if σ < at_infinity_tol
                at_infinity = true
            end
            if σ < strict_at_infinity_tol
                strict_at_infinity = true
            end
        end

        # 0 check
        if m - 2 * Δ > finite_tol && δ_x[i] > 0
            if σ < at_infinity_tol
                at_zero = true
            end
            if σ < strict_at_infinity_tol
                strict_at_zero = true
            end
        end

        # finite
        # Case a: val(x) = 0
        if abs(val_xᵢ) < finite_tol
            if δ_x[i] ≤ 0 || val_tx¹[i] < 0
                finite = false
            end
            # Case b: val(x) = val(tẋ) > 0
        elseif zero_is_finite && val_xᵢ > finite_tol
            # Δ = max(abs(Δval_x[i]), abs(Δval_tx¹[i]), abs(1 - val_xᵢ / m), abs(1 - val_tx¹[i] / m))
            if σ > finite_tol
                finite = false
            end
        else
            finite = false
        end
    end

    # Check singular only if finite
    winding_number_candidate = 1
    if finite
        # Compute best possible integer m
        d = Inf
        #TODO: Is there a way we can exclude some integers?
        for m = 1:max_winding_number
            d_m = 0.0
            for v in val_tx¹
                d_m = max(d_m, abs(round(m * v) - m * v))
            end
            if d_m < d
                winding_number_candidate, d = m, d_m
            end
        end

        singular =
            d < winding_number_candidate * singular_tol && winding_number_candidate > 1
    else
        singular = false
    end

    return (
        finite = finite,
        singular = singular,
        winding_number_candidate = winding_number_candidate,
        at_infinity = at_infinity,
        strict_at_infinity = strict_at_infinity,
        at_zero = at_zero,
        strict_at_zero = strict_at_zero,
    )
end
