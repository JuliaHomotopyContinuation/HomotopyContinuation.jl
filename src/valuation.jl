"""
    Valuation

A data structure for approximating the valuation of a solution path ``x(t)``.
"""
mutable struct Valuation
    val_x::Vector{Float64}
    val_tx¹::Vector{Float64}
    Δval_x::Vector{Float64}
    Δval_tx¹::Vector{Float64}
end

Valuation(n::Integer) = Valuation((zeros(n) for i = 1:4)...)
Valuation(x::AbstractVector) = Valuation(length(x))

function Base.show(io::IO, val::Valuation)
    print(io, "Valuation :")
    for field in [:val_x, :val_tx¹, :Δval_x, :Δval_tx¹]
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
    val
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


function update!(
    val::Valuation,
    tx::TaylorVector,
    t::Real,
)
    for i in eachindex(tx)
        xᵢ, x¹ᵢ, x²ᵢ, x³ᵢ = tx[i]

        ν, ν¹ = ν_ν¹(xᵢ, x¹ᵢ, 2x²ᵢ, t)
        val.val_x[i] = ν
        val.Δval_x[i] = t * ν¹

        ν, ν¹ = ν_ν¹(x¹ᵢ, 2x²ᵢ, 6x³ᵢ, t)
        val.val_tx¹[i] = ν + 1
        val.Δval_tx¹[i] = t * ν¹
    end

    val
end

function analyze(
    val::Valuation;
    finite_tol::Float64 = throw(UndefKeywordError(:finite_tol)),
    at_infinity_tol::Float64 = throw(UndefKeywordError(:at_infinity_tol)),
    zero_is_finite::Bool = throw(UndefKeywordError(:zero_is_finite)),
)
    @unpack val_x, val_tx¹, Δval_x, Δval_tx¹ = val

    at_infinity = false
    at_zero = false
    finite = true
    for (i, val_xᵢ) in enumerate(val_x)
        abs_Δ = max(abs(Δval_x[i]), abs(Δval_tx¹[i]))
        rel_Δ = abs_Δ / abs(val_xᵢ)
        m = max(abs(1.0 - val_tx¹[i] / val_xᵢ), rel_Δ)
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
            if abs_Δ > finite_tol || val_tx¹[i] ≤ 0
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
