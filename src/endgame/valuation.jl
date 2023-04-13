"""
    Valuation

A data structure for approximating the valuation of a solution path ``x(t)``.
"""
struct Valuation
    val_x::Vector{Float64}
    val_ẋ::Vector{Float64}
    Δval_x::Vector{Float64}
    Δval_ẋ::Vector{Float64}
end

Valuation(x::AbstractVector) = Valuation(length(x))
Valuation(n::Integer) = Valuation((zeros(n) for i = 1:4)...)

function Base.show(io::IO, val::Valuation)
    print(io, "Valuation :")
    for field in [:val_x, :val_ẋ, :Δval_x, :Δval_ẋ]
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
    val.val_ẋ .= 0.0
    val.Δval_x .= 0.0
    val.Δval_ẋ .= 0.0
    val
end


function ν_ν¹(x, ẋ, x², t, scaling_factor::Float64)
    u, v = reim(x)
    u¹, v¹ = reim(ẋ)
    u², v² = reim(x²)

    xx = abs2(x)

    μ = u * u¹ + v * v¹
    l = μ / xx

    μ¹ = u * u² + u¹^2 + v * v² + v¹^2
    l¹ = μ¹ / xx - 2 * l^2

    l /= scaling_factor
    l¹ /= scaling_factor^2

    t * l, t * l¹ + l
end


logabs(x) = log(fast_abs(x))

function update!(val::Valuation, pred::Predictor, t::Real)
    α = pred.scaling_factor[]

    for (i, (xᵢ, ẋᵢ, x²ᵢ, x³ᵢ)) in enumerate(pred.tx³)
        νᵢ, ν¹ = ν_ν¹(xᵢ, ẋᵢ, 2x²ᵢ, t, α)
        val.val_x[i] = νᵢ
        val.Δval_x[i] = t * ν¹

        val_ẋᵢ, ν¹ = ν_ν¹(ẋᵢ, 2x²ᵢ, 6x³ᵢ, t, α)
        val.val_ẋ[i] = val_ẋᵢ
        val.Δval_ẋ[i] = t * ν¹
    end

    val
end

module ValuationVerdict

@enum code begin
    Diverging
    Finite
    FiniteZero
    Unknown
end

end


Base.@kwdef struct ValuationAnalyzeResult
    verdict::ValuationVerdict.code
    coordinate_verdicts::Vector{ValuationVerdict.code}
end

function analyze_valuation(
    val::Valuation;
    tol_Δ::Float64 = 0.01,
    zero_is_finite::Bool = true,
)
    coordinate_verdicts = [ValuationVerdict.Diverging for _ = 1:length(val.val_x)]
    for i in eachindex(val.val_x)
        # Check if negative valuation
        if val.val_x[i] < -tol_Δ &&
           val.val_ẋ[i] < -tol_Δ &&
           abs(val.Δval_x[i]) < tol_Δ &&
           abs(val.Δval_ẋ[i]) < tol_Δ &&
           abs(val.val_x[i] - (val.val_ẋ[i] + 1)) < tol_Δ
            coordinate_verdicts[i] = ValuationVerdict.Diverging

            # Check if positive valuation
        elseif abs(val.val_x[i]) > tol_Δ &&
               abs(val.Δval_x[i]) < tol_Δ &&
               abs(val.Δval_ẋ[i]) < tol_Δ &&
               abs(1 - (val.val_ẋ[i] + 1) / val.val_x[i]) < tol_Δ

            coordinate_verdicts[i] =
                zero_is_finite ? ValuationVerdict.FiniteZero : ValuationVerdict.Diverging

            # Check if valuation 0
        elseif abs(val.val_x[i]) < tol_Δ && abs(val.Δval_x[i]) < tol_Δ
            coordinate_verdicts[i] = ValuationVerdict.Finite
        else
            coordinate_verdicts[i] = ValuationVerdict.Unknown
        end
    end

    verdict = if all(t -> t == ValuationVerdict.Finite, coordinate_verdicts)
        ValuationVerdict.Finite
    elseif all(
        t -> t == ValuationVerdict.Finite || t == ValuationVerdict.FiniteZero,
        coordinate_verdicts,
    )
        ValuationVerdict.FiniteZero
    elseif any(t -> t == ValuationVerdict.Diverging, coordinate_verdicts)
        ValuationVerdict.Diverging
    else
        ValuationVerdict.Unknown
    end

    ValuationAnalyzeResult(verdict, coordinate_verdicts)
end

is_finite(r::ValuationAnalyzeResult) =
    r.verdict == ValuationVerdict.Finite || r.verdict == ValuationVerdict.FiniteZero


function estimate_winding_number(
    val::Valuation;
    max_winding_number::Int = 10,
    tol_Δ::Float64 = 0.01,
)
    m = 0
    min_err = Inf
    for k = 1:max_winding_number
        err = check_winding_number(val, k; tol_Δ = tol_Δ)
        if err < min_err
            m = k
            min_err = err
        end
    end
    m, min_err
end

function check_winding_number(val::Valuation, m::Int; tol_Δ::Float64 = 0.01)
    err = 0.0
    for vᵢ in val.val_ẋ
        mvᵢ = m * vᵢ
        errᵢ = abs(round(mvᵢ) - mvᵢ)
        err = max(errᵢ, err)
    end
    err
end


function check_valuation_accuracy(val::Valuation, m::Int; tol_Δ::Float64 = 0.01)
    for i in eachindex(val.val_x)
        vᵢ = val.val_ẋ[i] + 1
        Δvᵢ = val.Δval_ẋ[i]

        # Valuation of ẋ has to be at least (1 / m) - 1
        if vᵢ < 1 / m - tol_Δ
            return false
        end

        mvᵢ = m * vᵢ
        Δmvᵢ = abs(m * Δvᵢ)

        if Δmvᵢ > tol_Δ * mvᵢ
            return false
        end
    end
    return true
end

function check_valuation_regular_finite(val::Valuation; tol_Δ::Float64 = 0.01)
    for i in eachindex(val.val_x)
        regular_finite =
            val.val_x[i] > -tol_Δ && # valuation non-negative
            abs(val.val_x[i] % 1) < tol_Δ && # valuation integer
            abs(val.val_ẋ[i] - (val.val_x[i] + 1)) < tol_Δ &&
            abs(val.Δval_x[i]) < tol_Δ &&
            abs(val.Δval_ẋ[i]) < tol_Δ
        if !regular_finite
            return false
        end
    end
    return true
end
