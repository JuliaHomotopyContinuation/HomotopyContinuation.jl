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
    w::Vector{Float64}
    Δω::Vector{Float64}
    σ::Vector{Float64}
    Δσ::Vector{Float64}
    t::Float64
end

Valuation(n::Integer) = Valuation(zeros(n), zeros(n), zeros(n), zeros(n), NaN)
Valuation(x::AbstractVector) = Valuation(length(x))

function Base.show(io::IO, val::Valuation)
    println(io, "Valuation :")
    for name in [:w, :Δω, :σ, :Δσ]
        println(
            io,
            " • ",
            name,
            " → ",
            round.(getfield(val, name); sigdigits = 8),
        )
    end
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", v::Valuation) = v

"""
    init!(val::Valuation)

Initialize the valuation estimator to start from scratch.
"""
function init!(val::Valuation)
    val.w .= 0
    val.Δω .= 0
    val.σ .= 0
    val.Δσ .= 0
    val.t = NaN
    val
end

"""
    ν(x, ẋ)

Computes the function ``ν(x) = (uu¹ + vv¹) / |x|²`` where `u,v = reim(x)`.
This is equivalent to computing ``d/dt log|x(t)|``.
"""
function ν_ν̇(x::Number, ẋ::Number, ẍ::Number)
    u, v = reim(x)
    u¹, v¹ = reim(ẋ)
    u², v² = reim(ẍ)

    x² = abs2(x)
    μ = u * u¹ + v * v¹
    μ̇ = u * u² + u¹^2 + v * v² + v¹^2
    ν = μ / x²
    # ν̇ = μ̇ / x² - 2(μ/x²)² and  μ/x² = ν
    ν̇ = μ̇ / x² - 2 * ν^2
    ν, ν̇
end

function update!(
    val::Valuation,
    x::AbstractVector,
    ẋ::AbstractVector,
    ẍ::AbstractVector,
    x3::AbstractVector,
    t::Real,
)
    for i in eachindex(x)
        νᵢ, ν̇ᵢ = ν_ν̇(x[i], ẋ[i], 2 * ẍ[i])
        wᵢ = νᵢ * t
        val.w[i] = wᵢ

        ν_σᵢ, ν̇_σᵢ = ν_ν̇(ẋ[i], 2 * ẍ[i], 6 * x3[i])
        val.σ[i] = ν_σᵢ * t + 1

        wᵢ, Δωᵢ = val.w[i], val.Δω[i]
        val.Δω[i] = abs((ν̇ᵢ * t + νᵢ) * t) / abs(val.σ[i])
        val.Δσ[i] = abs((ν̇_σᵢ * t + ν_σᵢ) * t) / abs(val.σ[i])
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

    @unpack w, Δω, σ, Δσ = val

    singular = false
    at_infinity = false
    indecisive = false

    for (wᵢ, Δωᵢ, σᵢ, Δσᵢ) in zip(w, Δω, σ, Δσ)
        # if one coordinate seems to diverge ignore all others
        # in particular don't care about indecisive results on other
        # coordinates
        if wᵢ + Δωᵢ < -tol && Δωᵢ < tol && abs(wᵢ - σᵢ) < tol
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
