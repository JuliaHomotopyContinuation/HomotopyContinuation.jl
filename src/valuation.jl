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
    t::Float64
end

function Valuation(n::Integer)
    Valuation(zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), NaN)
end
Valuation(x::AbstractVector) = Valuation(length(x))

function Base.show(io::IO, val::Valuation)
    println(io, "Valuation :")
    for name in [:val_x, :Δval_x, :val_x¹, :Δval_x¹, :val_x², :Δval_x²]
        println(
            io,
            " • ",
            name,
            " → ",
            "[",
            join(
                [Printf.@sprintf("%.4g", v) for v in getfield(val, name)],
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

"""
    ν(x, ẋ)

Computes the function ``ν(x) = (uu¹ + vv¹) / |x|²`` where `u,v = reim(x)`.
This is equivalent to computing ``d/dt log|x(t)|``.
"""
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
    tol::Float64,
    # at_infinity_tol::Float64 = throw(UndefKeywordError(:tol)),
)

    w, Δω, σ, Δσ = val.val_x, val.Δval_x, val.val_x¹, val.Δval_x¹

    singular = false
    at_infinity = false
    indecisive = false

    for (wᵢ, Δωᵢ, σᵢ, Δσᵢ) in zip(w, Δω, σ, Δσ)
        # if one coordinate seems to diverge ignore all others
        # in particular don't care about indecisive results on other
        # coordinates
        if wᵢ + Δωᵢ < -tol && Δωᵢ < tol && Δσᵢ < tol && abs(wᵢ - σᵢ) < tol
            return VAL_AT_INFINITY
        end

        if Δωᵢ > tol
            indecisive = true
        else
            # check singular
            # assume m < 20, otherwise we probably cannot get a
            # solution anyway due to accuracy limits
            if wᵢ ≥ 0.05 ||
               (-0.05 < wᵢ < 0.05 && abs(round(Int, σᵢ) - σᵢ) - Δσᵢ > 0.05)
                # test if σᵢ is fractional
                singular = true
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
