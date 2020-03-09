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
    val_x¹::Vector{Float64}
    val_x²::Vector{Float64}
    Δval_x::Vector{Float64}
    Δval_x¹::Vector{Float64}
    Δval_x²::Vector{Float64}
    δ_x::Vector{Float64}
    δ_x¹::Vector{Float64}
end

function Valuation(n::Integer)
    Valuation((zeros(n) for i = 1:8)...)
end
Valuation(x::AbstractVector) = Valuation(length(x))

function Base.show(io::IO, val::Valuation)
    print(io, "Valuation :")
    fs = [:val_x, :val_x¹, :val_x², :Δval_x, :Δval_x¹, :Δval_x², :δ_x, :δ_x¹]
    for field in fs
        print(
            io,
            "\n • ",
            field == :val_x ? "val_x " : field,
            " → ",
            "[",
            join(
                [Printf.@sprintf("%#.4g", v) for v in getfield(val, field)],
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
    val.val_x¹ .= 0
    val.val_x² .= 0
    val.Δval_x .= 0
    val.Δval_x¹ .= 0
    val.Δval_x² .= 0
    val.δ_x .= 0
    val.δ_x¹ .= 0
    val
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
    x⁴::AbstractVector,
    t::Real,
)
    for i in eachindex(x)
        xᵢ, x¹ᵢ, x²ᵢ, x³ᵢ, x⁴ᵢ = x[i], x¹[i], 2 * x²[i], 6 * x³[i], 24 * x⁴[i]
        ν, ν¹, ν² = ν_ν¹_ν²(xᵢ, x¹ᵢ, x²ᵢ, x³ᵢ, t)
        ν_t = t * ν¹
        ν_t¹ = t * ν² + ν¹
        δ_x = t * ν_t¹ / ν_t
        Δx = ν_t / δ_x
        val_x = ν
        if !isnan(Δx) && δ_x > 0
            val_x -= Δx
        end

        val.val_x[i] = val_x
        val.Δval_x[i] = Δx
        val.δ_x[i] = δ_x

        # valuation of t x¹
        ν, ν¹, ν² = ν_ν¹_ν²(x¹ᵢ, x²ᵢ, x³ᵢ, x⁴ᵢ, t)
        ν_t = t * ν¹
        ν_t¹ = t * ν² + ν¹

        δ_x¹ = t * ν_t¹ / ν_t
        Δx¹ = ν_t / δ_x¹
        val_x¹ = ν + 1
        if !isnan(Δx¹) && δ_x¹ > 0
            val_x¹ -= Δx¹
        end

        val.val_x¹[i] = val_x¹
        val.Δval_x¹[i] = Δx¹
        val.δ_x¹[i] = δ_x¹
        # valuation of t^2 x²
        ν, ν¹ = ν_ν¹(x²ᵢ, x³ᵢ, x⁴ᵢ, t)
        ν_t = t * ν¹
        Δx² = ν_t / δ_x¹
        # reuse δ_x¹ from tx¹
        val_x² = ν + 2
        if !isnan(Δx²) && δ_x¹ > 0
            val_x² -= Δx²
        end
        val.val_x²[i] = val_x²
        val.Δval_x²[i] = Δx²
    end

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

function mean_dev(a::Vararg{T,N}) where {T,N}
    m = sum(a) / N
    σ = sqrt(sum((a .- m) .^ 2)) / abs(m)
    m, σ
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
    val::Valuation;
    in_torus_tol::Float64 = throw(UndefKeywordError(:in_torus_tol)),
    at_infinity_tol::Float64 = throw(UndefKeywordError(:at_infinity_tol)),
    strict_at_infinity_tol::Float64 = throw(UndefKeywordError(
        :strict_at_infinity_tol,
    )),
    singular_tol::Float64 = throw(UndefKeywordError(:singular_tol)),
    max_winding_number::Int = throw(UndefKeywordError(:max_winding_number)),
)
    in_torus = true
    at_zero = false
    at_infinity = false
    strict_at_zero = false
    strict_at_infinity = false
    singular = false
    indecisive = false

    @unpack val_x, val_x¹, val_x², Δval_x, Δval_x¹, Δval_x², δ_x, δ_x¹ = val

    ε = 1 / max_winding_number

    for (i, val_xᵢ) in enumerate(val_x)
        if abs(val_xᵢ) < in_torus_tol &&
           (val_x¹[i] < -in_torus_tol || val_x²[i] < -in_torus_tol)
            indecisive = true
        end

        # at infinity check (val(x) < 0)
        #
        # -->  val_xᵢ, val_x¹[i], val_x²[i] coincide

        # check that they coincide by computing relative std. derivation
        m, σ = mean_dev(val_xᵢ, val_x¹[i], val_x²[i])
        if m < -in_torus_tol && σ < at_infinity_tol && δ_x[i] > 0 && δ_x¹[i] > 0
            at_infinity = true
            in_torus = false
        end
        if m < -in_torus_tol &&
           σ < strict_at_infinity_tol && δ_x[i] > 0 && δ_x¹[i] > 0
            strict_at_infinity = true
            in_torus = false
        end

        # zero check (val(x) > 0)

        # possible cases
        # 1) 0 < val(x) && val(x) ∉ ℕ₊
        #   -->  val_xᵢ, val_x¹[i], val_x²[i] coincide
        # 2) val(x) ∈ ℕ₊
        #   -->
        #   a) val_xᵢ = 1 : val_x¹[i]  -> 1, val_x²[i] -> 2
        #   b) val_xᵢ ≥ 2 : val_xᵢ, val_x¹[i], val_x²[i] coincide

        # Case 1) and 2)b
        if m > in_torus_tol && (
            σ < at_infinity_tol ||
            # Case 2b)
            at_infinity_tol >
            √((1 - val_xᵢ)^2 + (1 - val_x¹[i])^2 + (1 - 0.5 * val_x²[i])^2)
        ) && δ_x[i] > 0 && δ_x¹[i] > 0
            at_zero = true
            in_torus = false
            # torus check

            # val(x) = 0
        elseif abs(val_xᵢ) > in_torus_tol
            in_torus = false
            indecisive = true
        end

        if m > in_torus_tol && (
            σ < strict_at_infinity_tol || (
                # Case 2b)
                strict_at_infinity_tol >
                √((1 - val_xᵢ)^2 + (1 - val_x¹[i])^2 + (1 - 0.5 * val_x²[i])^2)
            )
        ) && δ_x[i] > 0 && δ_x¹[i] > 0
            strict_at_zero = true
            in_torus = false
        end


        # torus singular check (val(x) = 0, m > 1)
        #  the case if val_x¹[i], val_x²[i] coincide and val_x¹[i] ∉ ℕ₊
        if in_torus
            m₃, σ₃ = mean_dev(val_x¹[i], val_x²[i])
            is_fractional =
                !isnan(m₃) && abs(round(Int, m₃) - m₃) ≥ ε - singular_tol
            if m₃ > in_torus_tol &&
               is_fractional && σ₃ < singular_tol && δ_x[i] > 0 && δ_x¹[i] > 0
                singular = true
            end
        end

        # zero singular check (val(x) > 0, m > 1)
        #  the case if
        #    val_xᵢ > 0, val_xᵢ ∉ ℕ₊, and val_xᵢ, val_x¹[i], val_x²[i] coincide
        is_fractional = !isnan(m) && abs(round(Int, m) - m) ≥ ε - singular_tol
        if m > in_torus_tol &&
           is_fractional && σ < singular_tol && δ_x[i] > 0 && δ_x¹[i] > 0
            singular = true
        end
    end

    return (
        in_torus = in_torus,
        at_zero = at_zero,
        strict_at_zero = strict_at_zero,
        finite = !(at_infinity || indecisive),
        at_infinity = at_infinity,
        strict_at_infinity = strict_at_infinity,
        singular = singular,
    )
end
