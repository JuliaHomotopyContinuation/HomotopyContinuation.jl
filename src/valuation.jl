"""
    Valuation

A data structure for approximating the valuation of a solution path ``x(t)``.
"""
mutable struct Valuation
    val_x::Vector{Float64}
    val_tẋ::Vector{Float64}
    Δval_x::Vector{Float64}
    Δval_tẋ::Vector{Float64}

    # Data for computing the val_x and val_tẋ and ν̈ by finite differences
    val_x_data::NTuple{2,Vector{Float64}}
    val_ẋ_data::NTuple{2,Vector{Float64}}
    logx_data::NTuple{2,Vector{Float64}}
    logẋ_data::NTuple{2,Vector{Float64}}
    logt_data::NTuple{2,Float64}
end

Valuation(x::AbstractVector) = Valuation(length(x))
Valuation(n::Integer) =
    Valuation((zeros(n) for i = 1:4)..., ((zeros(n), zeros(n)) for i = 1:4)..., (NaN, NaN))

function Base.show(io::IO, val::Valuation)
    print(io, "Valuation :")
    for field in [:val_x, :val_tẋ, :Δval_x, :Δval_tẋ]
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
    val.val_x .= 0.0
    val.val_tẋ .= 0.0
    val.Δval_x .= 0.0
    val.Δval_tẋ .= 0.0
    val.val_x_data[1] .= val.val_x_data[2] .= 0.0
    val.val_ẋ_data[1] .= val.val_ẋ_data[2] .= 0.0
    val.logx_data[1] .= val.logx_data[2] .= 0.0
    val.logẋ_data[1] .= val.logẋ_data[2] .= 0.0
    val.logt_data = (NaN, NaN)
    val
end


function ν(x, ẋ, t)
    u, v = reim(x)
    u¹, v¹ = reim(ẋ)
    xx = abs2(x)

    μ = u * u¹ + v * v¹
    l = μ / xx

    t * l
end


function ν_ν¹(x, ẋ, x², t)
    u, v = reim(x)
    u¹, v¹ = reim(ẋ)
    u², v² = reim(x²)

    xx = abs2(x)

    μ = u * u¹ + v * v¹
    l = μ / xx

    μ¹ = u * u² + u¹^2 + v * v² + v¹^2
    l¹ = μ¹ / xx - 2 * l^2

    t * l, t * l¹ + l
end


logabs(x) = log(fast_abs(x))

function update!(val::Valuation, pred::Predictor{AD{N}}, t::Real) where {N}
    # use analytic expressions if derivatives are computed with automatic differentiation
    if N ≥ 3
        update_analytic!(val, pred.tx³, t)
        return val
    end

    # otherwise use finite difference scheme
    logt = log(t)
    ν₂, ν₁ = val.val_x_data
    logx₂, logx₁, = val.logx_data
    logẋ₂, logẋ₁, = val.logẋ_data
    logt₂, logt₁ = val.logt_data
    val_ẋ₂, val_ẋ₁ = val.val_ẋ_data
    for (i, (xᵢ, ẋᵢ)) in enumerate(pred.tx¹)
        logxᵢ = logabs(xᵢ)
        νᵢ = ν(xᵢ, ẋᵢ, t)
        Δνᵢ = finite_diff(νᵢ, logt, ν₂[i], logt₂, ν₁[i], logt₁)

        logẋᵢ = logabs(ẋᵢ)
        val_ẋᵢ = finite_diff(logẋᵢ, logt, logẋ₂[i], logt₂, logẋ₁[i], logt₁)
        Δval_ẋᵢ = finite_diff(val_ẋᵢ, logt, val_ẋ₂[i], logt₂, val_ẋ₁[i], logt₁)

        val.val_x[i] = νᵢ
        val.Δval_x[i] = Δνᵢ
        val.val_tẋ[i] = val_ẋᵢ + 1
        val.Δval_tẋ[i] = Δval_ẋᵢ

        ν₂[i], ν₁[i] = νᵢ, ν₂[i]
        logx₂[i], logx₁[i] = logxᵢ, logx₂[i]
        logẋ₂[i], logẋ₁[i] = logẋᵢ, logẋ₂[i]
        val_ẋ₂[i], val_ẋ₁[i] = val_ẋᵢ, val_ẋ₂[i]


    end
    val.logt_data = (logt, logt₂)

    val
end

function update_analytic!(val::Valuation, tx::TaylorVector, t::Real)
    for i in eachindex(tx)
        xᵢ, ẋᵢ, x²ᵢ, x³ᵢ = tx[i]

        ν, ν¹ = ν_ν¹(xᵢ, ẋᵢ, 2x²ᵢ, t)
        val.val_x[i] = ν
        val.Δval_x[i] = t * ν¹

        ν, ν¹ = ν_ν¹(ẋᵢ, 2x²ᵢ, 6x³ᵢ, t)
        val.val_tẋ[i] = ν + 1
        val.Δval_tẋ[i] = t * ν¹
    end

    val
end


"""
    finite_diff(ν₃, s₃, ν₂, s₂, ν₁, s₁)

Compute the derivative `ν̇` of a function ``ν(s)`` at `s₃` by using a finite
difference scheme. This uses the values `νᵢ = ν(sᵢ)` for ``i=1,…,3``.
Since we have a non-uniform grid, we need a more elaborate difference scheme.
The implementation follows the formulas derived in [^BS05].

[^BS05]: Bowen, M. K., and Ronald Smith. "Derivative formulae and errors for non-uniformly
  spaced points." Proceedings of the Royal Society A: Mathematical, Physical and Engineering
  Sciences 461.2059 (2005): 1975-1997.
"""
function finite_diff(ν, s, ν₂, s₂, ν₁, s₁)
    Δ₁, Δ₂, Δ₁₂ = s - s₁, s - s₂, s₁ - s₂
    ν̇ = (Δ₂ * ν₁) / (Δ₁₂ * Δ₁) - ((Δ₁₂ + Δ₂) * ν₂) / (Δ₁₂ * Δ₂) - (Δ₁₂ * ν) / (Δ₁ * Δ₂)
    # ν̈ = -2ν₁ / (Δ₁₂ * Δ₁) + 2ν₂ / (Δ₁₂ * Δ₂) + 2ν / (Δ₁ * Δ₂)
    ν̇
end



function analyze(
    val::Valuation;
    finite_tol::Float64 = throw(UndefKeywordError(:finite_tol)),
    at_infinity_tol::Float64 = throw(UndefKeywordError(:at_infinity_tol)),
    zero_is_finite::Bool = throw(UndefKeywordError(:zero_is_finite)),
)
    @unpack val_x, val_tẋ, Δval_x, Δval_tẋ = val

    at_infinity = false
    at_zero = false
    finite = true
    for (i, val_xᵢ) in enumerate(val_x)
        abs_Δ = max(abs(Δval_x[i]), abs(Δval_tẋ[i]))
        rel_Δ = abs_Δ / abs(val_xᵢ)
        m = max(abs(1.0 - val_tẋ[i] / val_xᵢ), rel_Δ)
        # ∞ check
        if val_xᵢ + m < -finite_tol && m < at_infinity_tol
            at_infinity = true
        end

        # 0 check
        if val_xᵢ - m > finite_tol && m < at_infinity_tol
            at_zero = true
        end

        # finite
        # Case a: val(x) = 0
        if abs(val_xᵢ) < finite_tol
            if abs_Δ > finite_tol || val_tẋ[i] ≤ 0
                finite = false
            end
            # Case b: val(x) = val(tẋ) > 0
        elseif zero_is_finite && val_xᵢ > finite_tol
            if m > finite_tol
                finite = false
            end
        else
            finite = false
        end
    end

    return (finite = finite, at_infinity = at_infinity, at_zero = at_zero)
end
