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
    Δval_x::Vector{Float64}
    val_x¹::Vector{Float64}
    Δval_x¹::Vector{Float64}
    val_x²::Vector{Float64}
    Δval_x²::Vector{Float64}
    val_x³::Vector{Float64}
    t::Float64
end

function Valuation(n::Integer)
    Valuation((zeros(n) for i = 1:7)..., NaN)
end
Valuation(x::AbstractVector) = Valuation(length(x))

function Base.show(io::IO, val::Valuation)
    print(io, "Valuation :")
    for name in [:val_x, :val_x¹, :val_x², :val_x³, :Δval_x, :Δval_x¹, :Δval_x²]
        print(
            io,
            "\n • ",
            name == :val_x ? "val_x " : name,
            " → ",
            "[",
            join(
                [Printf.@sprintf("%#.5g", v) for v in getfield(val, name)],
                ", ",
            ),
            "]",
        )
    end
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", v::Valuation) = v

"""
    init!(val::Valuation)

Initialize the valuation estimator to start from scratch.
"""
function init!(val::Valuation)
    val.val_x .= 0
    val.Δval_x .= 0
    val.val_x¹ .= 0
    val.Δval_x¹ .= 0
    val.val_x² .= 0
    val.Δval_x² .= 0
    val.t = NaN
    val
end

function ν(x::Number, ẋ::Number, t)
    u, v = reim(x)
    u¹, v¹ = reim(ẋ)
    x² = abs2(x)
    μ = u * u¹ + v * v¹
    t * μ / x²
end

function ν_Δν(x::Number, ẋ::Number, ẍ::Number, t)
    u, v = reim(x)
    u¹, v¹ = reim(ẋ)
    u², v² = reim(ẍ)

    x² = abs2(x)
    μ = u * u¹ + v * v¹
    μ̇ = u * u² + u¹^2 + v * v² + v¹^2
    ν = μ / x²
    # ν̇ = μ̇ / x² - 2(μ/x²)² and  μ/x² = ν
    ν̇ = μ̇ / x² - 2 * ν^2

    ν * t, (ν̇ * t + ν) * t
end

function update!(
    val::Valuation,
    x::AbstractVector,
    x¹::AbstractVector,
    x²::AbstractVector,
    x³::AbstractVector,
    x⁴::AbstractVector,
    t::Real,
)
    for i in eachindex(x)
        val_xᵢ, Δval_xᵢ = ν_Δν(x[i], x¹[i], 2 * x²[i], t)
        val.val_x[i] = val_xᵢ
        val.Δval_x[i] = abs(Δval_xᵢ)

        val_x¹ᵢ, Δval_x¹ᵢ = ν_Δν(x¹[i], 2 * x²[i], 6 * x³[i], t)
        val.val_x¹[i] = val_x¹ᵢ + 1
        val.Δval_x¹[i] = abs(Δval_x¹ᵢ)

        val_x²ᵢ, Δval_x²ᵢ = ν_Δν(2 * x²[i], 6 * x³[i], 24 * x⁴[i], t)
        val.val_x²[i] = val_x²ᵢ + 2
        val.Δval_x²[i] = abs(Δval_x²ᵢ)

        val.val_x³[i] = ν(6 * x³[i], 24 * x⁴[i], t) + 3
    end
    val.t = t

    val
end

"""
    enum ValuationVerdict

The possible states the [`judge`](@ref) function returns:

* `VAL_INDECISIVE`: The estimates are not trustworthy enough.
* `VAL_FINITE`: The valuation indicates that the path is finite.
* `VAL_AT_INFINITY`: The valuation indicates that the path is diverging.
"""
@enum ValuationVerdict begin
    VAL_INDECISIVE
    VAL_FINITE
    VAL_FINITE_SINGULAR
    VAL_AT_INFINITY
end

"""
    judge(val::Valuation; tol, tol_at_infinity)::ValuationVerdict

Judge the current valuation to determine whether a path is diverging. Returns
a [`ValuationVerdict`](@ref).
A path is diverging if one entry of the valuation is negative. To assure that
the estimate is trust worthy we require that the first and second derivative
of [`ν``](@ref) is smaller than `tol_at_infinity`.
The `tol` is used to estimate whether a finite valuation is trustworthy.
"""
function judge(
    val::Valuation,
    tol::Float64;
    at_infinity_tol::Float64 = throw(UndefKeywordError(:at_infinity_tol)),
    max_winding_number::Int = throw(UndefKeywordError(:max_winding_number)),
)
    singular = false
    at_infinity = false
    indecisive = false

    @unpack val_x, Δval_x, val_x¹, Δval_x¹, val_x², Δval_x², val_x³ = val

    ε = 1 / max_winding_number

    for i = 1:length(val_x)
        # if one coordinate seems to diverge ignore all others
        # in particular don't care about indecisive results on other
        # coordinates

        # At infinity check
        # Idea: If val ≠ 0 then then all val_x values coincide. We consider
        #       to have correct
        #
        val_mean = (val_x[i] + val_x¹[i] + val_x²[i] + val_x³[i]) / 4
        if val_mean < -ε && abs(val_x[i] - val_mean) < at_infinity_tol &&
           abs(val_x¹[i] - val_mean) < at_infinity_tol &&
           abs(val_x²[i] - val_mean) < at_infinity_tol &&
           abs(val_x³[i] - val_mean) < at_infinity_tol &&
           min(Δval_x[i], Δval_x¹[i], Δval_x²[i]) < at_infinity_tol
            return VAL_AT_INFINITY
        end

        # If we have a finite value then val_x[i] ≈ 0 and
        # the others have to conincide
        if max(Δval_x[i], Δval_x¹[i], Δval_x²[i]) > tol
            indecisive = true
        else
            # check solution in torus
            if val_x[i] > -ε + tol
                if val_x³[i] > ε + tol &&
                   # check whether val is fraction
                   abs(round(Int, val_x³[i]) - val_x³[i]) > ε
                    singular = true
                end
            else
                # probably diverging but not yet sure
                indecisive = true
            end
        end
    end
    indecisive && return VAL_INDECISIVE
    singular && return VAL_FINITE_SINGULAR

    VAL_FINITE
end
