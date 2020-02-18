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
    w::Vector{Float64}
    ẇ::Vector{Float64}
    δ::Vector{Float64}
    δ̇::Vector{Float64}
end

Valuation(n::Integer) = Valuation(zeros(n), zeros(n), zeros(n), zeros(n))
Valuation(x::AbstractVector) = Valuation(length(x))

function Base.show(io::IO, val::Valuation)
    println(io, "Valuation :")
    for name in [:w, :ẇ, :δ, :δ̇]
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
    val.w .= NaN
    val.ẇ .= NaN
    val
end

"""
    ν(x, ẋ)

Computes the function ``ν(x) = (uu̇ + vv̇) / |x|²`` where `u,v = reim(x)`.
This is equivalent to computing ``d/dt log|x(t)|``.
"""
function ν(x::Complex, ẋ::Complex)
    u, v = reim(x)
    u̇, v̇ = reim(ẋ)
    (u * u̇ + v * v̇) / abs2(x)
end

function ν_ν̇(x::Complex, ẋ::Complex, ẍ::Complex)
    u, v = reim(x)
    u̇, v̇ = reim(ẋ)
    ü, v̈ = reim(ẋ)

    x² = abs2(x)
    μ = u * u̇ + v * v̇
    μ̇ = u * ü + u̇^2 + v * v̈ + v̇^2
    ν = μ / x²
    # ν̇ = μ̇ / x² - 2μ²/(x²)² and  μ²/(x²)² = ν^2
    ν̇ = μ̇ / x² - 2 * ν^2
    ν, ν̇
end

function update!(
    val::Valuation,
    x::AbstractVector,
    ẋ::AbstractVector,
    ẍ::AbstractVector,
    x3::AbstractVector,
    t::Real;
    logarithmic_scale::Bool = false,
)
    for i in eachindex(x)
        νᵢ, ν̇ᵢ = ν_ν̇(x[i], ẋ[i], 2ẍ[i])
        if logarithmic_scale
            val.w[i] = νᵢ
            val.ẇ[i] = ν̇ᵢ
        else
            val.w[i] = νᵢ * t
            val.ẇ[i] = ν̇ᵢ * t^2 / (2 - 2log(t))
        end

        if logarithmic_scale
            val.δ[i] = νᵢ
            val.δ̇[i] = ν̇ᵢ
        else
            ν_δᵢ, ν̇_δᵢ = ν_ν̇(ẋ[i], 2ẍ[i], 6x3[i])
            val.δ[i] = ν_δᵢ * t + 1
            val.δ̇[i] = ν̇_δᵢ * t^2 / (2 - 2log(t))
        end
    end
    val
end
