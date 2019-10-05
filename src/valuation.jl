###############
## VALUATION ##
###############

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
struct Valuation
    # Estimate
    ν::Vector{Float64}
    ν̇::Vector{Float64}
    ν̈::Vector{Float64}
    s::Base.RefValue{Float64}
    affine::Bool
    # Data for computing the ν̇ and ν̈ by finite differences
    ν_data::NTuple{2,Vector{Float64}}
    s_data::Base.RefValue{NTuple{2,Float64}}
end

function Valuation(n::Integer, affine::Bool = true)
    v = zeros(n)
    ν̇ = zeros(n)
    ν̈ = zeros(n)
    ν_data = (zeros(n), zeros(n))
    s_data = (NaN, NaN)
    Valuation(v, ν̇, ν̈, Ref(NaN), affine, ν_data, Ref(s_data))
end
function Valuation(x::ProjectiveVectors.PVector; affine::Bool = true)
    if affine
        Valuation(length(x) - length(ProjectiveVectors.dims(x)), affine)
    else
        Valuation(length(x), affine)
    end
end
Valuation(x::Vector; affine::Bool = true) = Valuation(length(x))

function Base.show(io::IO, val::Valuation)
    println(io, typeof(val), ":")
    println(io, " • s → ", round.(val.s[]; sigdigits = 4))
    for name in [:ν, :ν̇, :ν̈]
        println(io, " • ", name, " → ", round.(getfield(val, name); sigdigits = 4))
    end
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", v::Valuation) = v

"""
    init!(val::Valuation)

Initialize the valuation estimator to start from scratch.
"""
function init!(val::Valuation)
    val.ν .= NaN
    val.ν̇ .= NaN
    val.ν̈ .= NaN
    val.s[] = NaN
    val.s_data[] = (NaN, NaN)
    val.ν_data[1] .= NaN
    val.ν_data[2] .= NaN
    val
end

#############################
## Valuation approximation ##
#############################

# First to the coordinate wise functions

"""
    ν(z, ż)

Computes the function ``ν(z) = -(xẋ + yẏ) / |z|²`` where `x,y = reim(z)`.
This is equivalent to computing ``d/ds log|x(s)|``.
"""
function ν(z::Complex, ż::Complex)
    x, y = reim(z)
    ẋ, ẏ = reim(ż)
    -(x * ẋ + y * ẏ) / abs2(z)
end

"""
    ν_ν̇_ν̈(z::Complex, ż::Complex, z̈::Complex, z³::Complex)

Computes the function [`ν`](@ref) and the derivaties ``ν̇``, ``ν̈`` by using the analytic
derivatives.
"""
function ν_ν̇_ν̈(z::Complex, ż::Complex, z̈::Complex, z³::Complex)
    x, y = reim(z)
    ẋ, ẏ = reim(ż)
    ẍ, ÿ = reim(z̈)
    x³, y³ = reim(z³)

    μ = x * ẋ + y * ẏ
    # The following is just applying product rule properly
    μ̇ = x * ẍ + ẋ^2 + y * ÿ + ẏ^2
    μ̈ = x * x³ + 3 * ẋ * ẍ + y * y³ + 3 * ẏ * ÿ

    z² = abs2(z)
    ν = -μ / z²
    # ν̇ = -μ̇ / z² + 2μ²/(z²)² and  μ²/(z²)² = ν^2
    ν̇ = 2 * ν^2 - μ̇ / z²
    ν̈ = 4 * ν * ν̇ - μ̈ / z² + 2 * μ̇ * μ / z²^2
    ν, ν̇, ν̈
end

"""
    ν̇_ν̈(z::Complex, ż::Complex, s::Float64, ν₂, s₂, ν₁, s₁)

Computes the function [`ν`](@ref) and the derivaties ``ν̇``, ``ν̈`` by using a finite
difference scheme. This uses the value `ν₂` of ``ν`` at `s₂` and `ν₁` of ``ν`` at `s₁` with
``s > s₂ > s₁``.
Since we have a non-uniform grid, we need a more elaborate difference scheme.
The implementation follows the formulas derived in [^BS05].

[^BS05]: Bowen, M. K., and Ronald Smith. "Derivative formulae and errors for non-uniformly
  spaced points." Proceedings of the Royal Society A: Mathematical, Physical and Engineering
  Sciences 461.2059 (2005): 1975-1997.
"""
function ν̇_ν̈(ν::Float64, s::Float64, ν₂::Float64, s₂::Float64, ν₁::Float64, s₁::Float64)
    Δ₁, Δ₂, Δ₁₂ = s - s₁, s - s₂, s₁ - s₂
    ν̇ = (Δ₂ * ν₁) / (Δ₁₂ * Δ₁) - ((Δ₁₂ + Δ₂) * ν₂) / (Δ₁₂ * Δ₂) - (Δ₁₂ * ν) / (Δ₁ * Δ₂)
    ν̈ = -2ν₁ / (Δ₁₂ * Δ₁) + 2ν₂ / (Δ₁₂ * Δ₂) + 2ν / (Δ₁ * Δ₂)
    ν̇, ν̈
end

"""
    update!(val::Valuation, x, ẋ, s,
            predictor::AbstractPredictorCache = NullPredictorCache())

Compute new approximations of the valuation of the path ``x(s)`` from `x` and `ẋ` at `s`.
This queries the given predictor cache `predictor` for the second and third derivative
of the path ``x(s)``. If this information is available it computes the derivatives
of ``ν`` analytically otherwise a finite differene scheme is used.
"""
function update!(
    val::Valuation,
    z::AbstractVector,
    ż,
    s,
    predictor::AbstractPredictorCache = NullPredictorCache(),
)
    @unpack ν̇, ν̈, ν_data, s_data = val

    ν₂, ν₁ = ν_data
    s₂, s₁ = s_data[]

    if val.affine && isa(z, PVector)
        k = 1
        for (rⱼ, j) in ProjectiveVectors.dimension_indices_homvars(z)
            νⱼ = ν(z[j], ż[j])
            for i in rⱼ
                ν_k = ν(z[i], ż[i]) - νⱼ
                if !isnan(s₁)
                    ν̇_k, ν̈_k = ν̇_ν̈(ν_k, s, ν₂[k], s₂, ν₁[k], s₁)
                    ν̇[k], ν̈[k] = ν̇_k, ν̈_k
                end
                # update ν and shift the previous ones
                val.ν[k], ν₂[k], ν₁[k] = ν_k, val.ν[k], ν₂[k]
                k += 1
            end
        end
    else
        for i = 1:length(val.ν)
            # try to compute the values using higher derivatives from the predictor
            z̈ = second_derivative(predictor)
            z³ = third_derivative(predictor)
            if !isnothing(z̈) && !isnothing(z³)
                νᵢ, ν̇ᵢ, ν̈ᵢ = ν_ν̇_ν̈(z[i], ż[i], z̈[i], z³[i])
                if !isnan(s₁)
                    # The estimates can be sensitive to numerical errors, so we still compare
                    # against the numerical derivatives since this is cheap anyway
                    ν̇ᵢ′, ν̈ᵢ′ = ν̇_ν̈(νᵢ, s, ν₂[i], s₂, ν₁[i], s₁)
                    ν̇[i] = abs(ν̇ᵢ) < abs(ν̇ᵢ′) ? ν̇ᵢ : ν̇ᵢ′
                    ν̈[i] = abs(ν̈ᵢ) < abs(ν̈ᵢ′) ? ν̈ᵢ : ν̈ᵢ′
                else
                    ν̇[i], ν̈[i] = ν̇ᵢ, ν̈ᵢ
                end
            else
                νᵢ = ν(z[i], ż[i])
                if !isnan(s₁)
                    ν̇ᵢ, ν̈ᵢ = ν̇_ν̈(νᵢ, s, ν₂[i], s₂, ν₁[i], s₁)
                    ν̇[i], ν̈[i] = ν̇ᵢ, ν̈ᵢ
                end
            end
            # update ν and shift the previous ones
            val.ν[i], ν₂[i], ν₁[i] = νᵢ, val.ν[i], ν₂[i]
        end
    end

    # Shift the s
    s_data[] = (val.s[], s₂)
    val.s[] = s
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
    VAL_AT_INFINITY
end

"""
    judge(val::Valuation; tol = 1e-3, tol_at_infinity = 1e-4)::ValuationVerdict

Judge the current valuation to determine whether a path is diverging. Returns
a [`ValuationVerdict`](@ref).
A path is diverging if one entry of the valuation is negative. To assure that
the estimate is trust worthy we require that the first and second derivative
of [`ν``](@ref) is smaller than `tol_at_infinity`.
The `tol` is used to estimate whether a finite valuation is trustworthy.
"""
function judge(val::Valuation; tol::Float64 = 1e-3, tol_at_infinity::Float64 = 1e-4)
    @unpack ν, ν̇, ν̈ = val
    finite = true
    indecisive = false
    for (i, νᵢ) in enumerate(ν)
        if νᵢ < -max(0.01, tol)
            finite = false
        end

        if νᵢ < -0.05 && abs(ν̇[i]) < tol_at_infinity && abs(ν̈[i]) < tol_at_infinity
            return VAL_AT_INFINITY
        end

        if abs(ν̇[i]) > tol || abs(ν̈[i]) > tol
            indecisive = true
        end
    end
    finite && !indecisive && return VAL_FINITE

    VAL_INDECISIVE
end
