module DoubleDouble

import Printf

export DoubleF64, ComplexDF64

import Base: +, -, *, /, ^, <, ==, <=

# Basic building blocks
"""
    quick_two_sum(a, b)

Computes `s = fl(a+b)` and `e = err(a+b)`. Assumes `|a| ≥ |b|`.
"""
@inline function quick_two_sum(a::Float64, b::Float64)
    s = a + b
    e = b - (s - a)
    s, e
end

"""
    two_sum(a, b)

Computes `s = fl(a+b)` and `e = err(a+b)`.
"""
@inline function two_sum(a::Float64, b::Float64)
    s = a + b
    v = s - a
    e = (a - (s - v)) + (b - v)

    s, e
end

"""
    quick_two_diff(a, b)

Computes `s = fl(a-b)` and `e = err(a-b)`.  Assumes `|a| ≥ |b|`.
"""
@inline function quick_two_diff(a::Float64, b::Float64)
    s = a - b
    e = (a - s) - b
    s, e
end

"""
    two_diff(a, b)

Computes `s = fl(a-b)` and `e = err(a-b)`.
"""
@inline function two_diff(a::Float64, b::Float64)
    s = a - b
    v = s - a
    e = (a - (s - v)) - (b + v)

    s, e
end

"""
    two_prod(a, b)

Computes `s = fl(a*b)` and `e = err(a*b)`.
"""
@inline function two_prod(a::Float64, b::Float64)
    p = a * b
    e = fma(a, b, -p)
    p, e
end

"""
    two_square(a)

Computes `s = fl(a*a)` and `e = err(a*a)`. Faster than [`two_prod(a, a)`](@ref).
"""
@inline function two_square(a::Float64)
    p = a * a
    e = fma(a, a, -p)
    p, e
end


"""
    DoubleF64(x [, mode::ComputeMode]) <: AbstractFloat
"""
struct DoubleF64 <: AbstractFloat
    hi::Float64
    lo::Float64
end

DoubleF64(x::DoubleF64) = DoubleF64(x.hi, x.lo)
function DoubleF64(x::Float64)
    DoubleF64(x, isinf(x) ? Inf : 0.0)
end

function DoubleF64(x::Float32)
    DoubleF64(convert(Float64, x), isinf(x) ? Inf : 0.0)
end
function DoubleF64(x::Float16)
    DoubleF64(convert(Float64, x), isinf(x) ? Inf : 0.0)
end
function DoubleF64(x::Integer)
    DoubleF64(convert(Float64, x), isinf(x) ? Inf : 0.0)
end
function DoubleF64(x::BigFloat)
    z = convert(Float64, x)
    DoubleF64(z, convert(Float64, x - z))
end
function DoubleF64(x::Irrational)
    DoubleF64(big(x))
end
function DoubleF64(x::Rational)
    DoubleF64(numerator(x)) / DoubleF64(denominator(x))
end

const ComplexDF64 = Complex{DoubleF64}

hi(x::DoubleF64) = x.hi
lo(x::DoubleF64) = x.lo

Base.isbits(::DoubleF64) = true

Base.zero(::DoubleF64) = DoubleF64(0.0, 0.0)
Base.zero(::Type{DoubleF64}) = DoubleF64(0.0, 0.0)

Base.one(::DoubleF64) = DoubleF64(1.0, 0.0)
Base.one(::Type{DoubleF64}) = DoubleF64(1.0, 0.0)

Base.convert(::Type{T}, a::DoubleF64) where {T<:AbstractFloat} = convert(T, a.hi)
Base.convert(::Type{BigFloat}, a::DoubleF64) = big(a.hi) + big(a.lo)
Base.convert(::Type{T}, a::DoubleF64) where {T<:Integer} = convert(T, a.hi)
Base.convert(::Type{Integer}, a::DoubleF64) = convert(Int64, a.hi)
Base.convert(::Type{BigInt}, a::DoubleF64) =
    convert(BigInt, big(a.hi) + big(a.lo))

Base.convert(::Type{DoubleF64}, x::DoubleF64) = x
Base.convert(::Type{DoubleF64}, x::AbstractFloat) = DoubleF64(x)
Base.convert(::Type{DoubleF64}, x::Irrational) = DoubleF64(x)
Base.convert(::Type{DoubleF64}, x::Integer) = DoubleF64(x)

Base.promote_rule(::Type{DoubleF64}, ::Type{<:Integer}) = DoubleF64
Base.promote_rule(::Type{DoubleF64}, ::Type{BigInt}) = BigFloat
Base.promote_rule(::Type{DoubleF64}, ::Type{BigFloat}) = BigFloat
Base.promote_rule(::Type{DoubleF64}, ::Type{Float64}) = DoubleF64
Base.promote_rule(::Type{DoubleF64}, ::Type{Float32}) = DoubleF64
Base.promote_rule(::Type{DoubleF64}, ::Type{Float16}) = DoubleF64

Base.big(x::DoubleF64) = big(x.hi) + big(x.lo)

const double_π = DoubleF64(pi)
const double_2pi = DoubleF64(2 * BigFloat(Base.pi))
const double_pi = DoubleF64(BigFloat(Base.pi))
const double_pi2 = DoubleF64(BigFloat(Base.pi) * 0.5)
const double_pi4 = DoubleF64(BigFloat(Base.pi) * 0.25)
const double_pi16 = DoubleF64(BigFloat(Base.pi) * (1 / 16))
const double_3pi4 = DoubleF64(BigFloat(Base.pi) * 0.75)
const double_nan = DoubleF64(NaN, NaN)
const double_inf = DoubleF64(Inf)
const double_e = convert(DoubleF64, Base.MathConstants.e)

const double_log2 = convert(DoubleF64, log(BigFloat(2.0)))
const double_log10 = convert(DoubleF64, log(BigFloat(10.0)))

const double_eps = 4.93038065763132e-32 # 2^-104

#
# ADDITION
#
"""
    wide_add(a::Float64, b::Float64)

Add two `Float64`s with `DoubleF64` precision.
"""
@inline function wide_add(a::Float64, b::Float64)
    hi, lo = two_sum(a, b)

    DoubleF64(hi, lo)
end

@inline function +(a::DoubleF64, b::Float64)
    hi, lo = two_sum(a.hi, b)
    lo += a.lo
    hi, lo = quick_two_sum(hi, lo)

    DoubleF64(hi, lo)
end
+(a::Float64, b::DoubleF64) = b + a
+(a::Integer, b::DoubleF64) = b + float(a)
+(a::DoubleF64, b::Integer) = a + float(b)

@inline function +(a::DoubleF64, b::DoubleF64)
    hi, lo = two_sum(a.hi, b.hi)
    lo += (a.lo + b.lo)
    hi, lo = quick_two_sum(hi, lo)

    DoubleF64(hi, lo)
end

#
# SUBSTRACTION
#
"""
    wide_sub(a::Float64, b::Float64)

Subtract two `Float64`s with `DoubleF64` precision.
"""
@inline function wide_sub(a::Float64, b::Float64)
    hi, lo = two_diff(a, b)

    DoubleF64(hi, lo)
end

@inline function -(a::DoubleF64, b::Float64)
    hi, lo = two_diff(a.hi, b)
    lo += a.lo
    hi, lo = quick_two_sum(hi, lo)

    DoubleF64(hi, lo)
end

@inline function -(a::Float64, b::DoubleF64)
    hi, lo = two_diff(a, b.hi)
    lo -= b.lo
    hi, lo = quick_two_sum(hi, lo)

    DoubleF64(hi, lo)
end

@inline function -(a::DoubleF64, b::DoubleF64)
    hi, lo = two_diff(a.hi, b.hi)
    lo += a.lo
    lo -= b.lo
    hi, lo = quick_two_sum(hi, lo)

    DoubleF64(hi, lo)
end
-(a::DoubleF64) = DoubleF64(-a.hi, -a.lo)


#
# MULTIPLICATION
#
"""
    wide_mutiply(a::Float64, b::Float64)

Multiply two `Float64`s with `DoubleF64` precision.
"""
@inline function wide_mul(a::Float64, b::Float64)
    hi, lo = two_prod(a, b)
    DoubleF64(hi, lo)
end

@inline function *(a::DoubleF64, b::Float64)
    p1, p2 = two_prod(a.hi, b)
    p2 += a.lo * b
    p1, p2 = quick_two_sum(p1, p2)

    DoubleF64(p1, p2)
end
*(a::Float64, b::DoubleF64) = b * a
*(a::Integer, b::DoubleF64) = b * float(a)
*(a::Bool, b::DoubleF64) = b * convert(Float64, a)
*(a::DoubleF64, b::Integer) = a * float(b)


@inline function *(a::DoubleF64, b::DoubleF64)
    p1, p2 = two_prod(a.hi, b.hi)
    p2 += a.hi * b.lo + a.lo * b.hi
    p1, p2 = quick_two_sum(p1, p2)

    DoubleF64(p1, p2)
end


#
# DIVISION
#
"""
    wide_div(a::Float64, b::Float64)

Divide two `Float64`s with `DoubleF64` precision.
"""
@inline function wide_div(a::Float64, b::Float64)
    q1 = a / b
        # Compute a - q1 * b
    p1, p2 = two_prod(q1, b)
    s, e = two_diff(a, p1)
    e -= p2

        # get next approximation
    q2 = (s + e) / b

    s, e = quick_two_sum(q1, q2)

    DoubleF64(s, e)
end

@inline function /(a::DoubleF64, b::Float64)
    q1 = a.hi / b

        # Compute  this - q1 * d
    p1, p2 = two_prod(q1, b)
    s, e = two_diff(a.hi, p1)
    e += a.lo
    e -= p2

        # get next approximation.
    q2 = (s + e) / b

        # renormalize
    hi, lo = quick_two_sum(q1, q2)

    DoubleF64(hi, lo)
end


@inline function /(a::DoubleF64, b::DoubleF64)
    q1 = a.hi / b.hi  # approximate quotient

        # compute  this - q1 * dd
    r = b * q1
    s1, s2 = two_diff(a.hi, r.hi)
    s2 -= r.lo
    s2 += a.lo

        # get next approximation
    q2 = (s1 + s2) / b.hi

        # renormalize
    hi, lo = quick_two_sum(q1, q2)

    DoubleF64(hi, lo)
end

/(a::Float64, b::DoubleF64) = DoubleF64(a) / b
/(a::Integer, b::DoubleF64) = float(a) / b
/(a::DoubleF64, b::Integer) = a / float(b)


Base.inv(a::DoubleF64) = 1.0 / a

Base.rem(a::DoubleF64, b::DoubleF64) = a - round(a / b) * b
@inline function Base.divrem(a::DoubleF64, b::DoubleF64)
    n = round(a / b)
    n, a - n * b
end

function Base.mod(x::DoubleF64, y::DoubleF64)
    n = round(a / b)
    return (a - b * n)
end

#
# POWERS
#

"""
    square(x::DoubleF64)

Compute `x * x` in a more efficient way.
"""
@inline function square(a::DoubleF64)
    p1, p2 = two_square(a.hi)
    p2 += 2.0 * a.hi * a.lo
    p2 += a.lo * a.lo
    hi, lo = quick_two_sum(p1, p2)
    return DoubleF64(hi, lo)
end
square(a::Float64) = wide_square(a)
"""
    wide_square(x::Float64)

Convert `x` to a DoubleF64 and then compute `x*x`.
"""
@inline function wide_square(a::Float64)
    hi, lo = two_square(a)
    return DoubleF64(hi, lo)
end
    # Implementation adapted from Base
@inline function power_by_squaring(x::DoubleF64, p::Integer)
    if p == 1
        return copy(x)
    elseif p == 0
        return one(x)
    elseif p == 2
        return square(x)
    end
    P = abs(p)
    t = trailing_zeros(P) + 1
    P >>= t
    while (t -= 1) > 0
        x = square(x)
    end
    y = x
    while P > 0
        t = trailing_zeros(P) + 1
        P >>= t
        while (t -= 1) >= 0
            x = square(x)
        end
        y *= x
    end
    if (p < 0)
        return 1.0 / y
    end
    y
end

^(a::DoubleF64, p::Integer) = power_by_squaring(a, p)
^(a::DoubleF64, b::Real) = iszero(a) ? one(a) : exp(b * log(a))

@inline function Base.sqrt(a::DoubleF64)
    # Strategy:  Use Karp's trick:  if x is an approximation
    # to sqrt(a), then
    #
    #  sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)
    #
    # The approximation is accurate to twice the accuracy of x.
    # Also, the multiplication (a*x) and [-]*x can be done with
    # only half the precision.
    if a.hi < 0
        throw(DomainError("sqrt will only return a complex result if called with a complex argument."))
    end
    if iszero(a)
        return zero(a)
    end

    x = inv(sqrt(a.hi))
    ax = a.hi * x

    wide_add(ax, (a - square(ax)).hi * (x * 0.5))
end

"""
    wide_sqrt(x::Float64)

Convert `x` to a DoubleF64 and then compute `sqrt(x)`.
"""
wide_sqrt(a::Float64) = sqrt(DoubleF64(a))


    # COMPARISON & EQUALITY
<(a::DoubleF64, b::DoubleF64) = a.hi + a.lo < b.hi + b.lo
<(a::DoubleF64, b::Float64) = a.hi < b || (a.hi == b) && a.lo < 0.0
<(a::Float64, b::DoubleF64) = a < b.hi || (a == b.hi) && b.lo > 0.0

Base.isless(a::DoubleF64, b::DoubleF64) = isless(a.hi + a.lo, b.hi + b.lo)

<=(a::DoubleF64, b::DoubleF64) = !(b < a)
<=(a::DoubleF64, b::Float64) = !(b < a)
<=(a::Float64, b::DoubleF64) = !(b < a)


==(a::DoubleF64, b::Float64) = a.hi == b && a.lo == 0.0
==(a::Float64, b::DoubleF64) = b == a
Base.iszero(a::DoubleF64) = a.hi == 0.0

Base.isone(a::DoubleF64) = a.hi == 1.0 && a.lo == 0.0
Base.abs(a::DoubleF64) = a.hi < 0.0 ? -a : a

Base.eps(::DoubleF64) = 4.93038065763132e-32 # 2^-104
Base.eps(::Type{DoubleF64}) = 4.93038065763132e-32 # 2^-104

    # Base.realmin(::Type{DoubleF64}) = 2.0041683600089728e-292 # = 2^(-1022 + 53)
    # Base.realmin(::Type{DoubleF64}) = 2.0041683600089728e-292 # = 2^(-1022 + 53)
    # Base.realmax(::Type{DoubleF64}) = DoubleF64(1.79769313486231570815e+308, 9.97920154767359795037e+291);
    # Base.realmax(::Type{DoubleF64}) = DoubleF64(1.79769313486231570815e+308, 9.97920154767359795037e+291);

Base.isnan(a::DoubleF64) = isnan(a.hi) || isnan(a.lo)
Base.isinf(a::DoubleF64) = isinf(a.hi)
Base.isfinite(a::DoubleF64) = isfinite(a.hi)

    #
    # ROUNDING
    #
@inline function Base.round(a::DoubleF64, r::RoundingMode = RoundNearest())
    hi = round(a.hi, r)
    lo = 0.0

    if hi == a.hi
            # High word is an integer already.  Round the low word.
        lo = round(a.lo, r)

            # Renormalize. This is needed if hi = some integer, lo = 1/2.
        hi, lo = quick_two_sum(hi, lo)
    else
            # High word is not an integer.
        if abs(hi - a.hi) == 0.5 && a.lo < 0.0
                # There is a tie in the high word, consult the low word to break the tie.
            hi -= 1.0
        end
    end

    DoubleF64(hi, lo)
end

@inline function Base.floor(a::DoubleF64)
    hi = floor(a.hi)
    lo = 0.0

    if hi == a.hi
        lo = floor(a.lo)
        hi, lo = quick_two_sum(hi, lo)
    end

    DoubleF64(hi, lo)
end

@inline function Base.floor(::Type{I}, a::DoubleF64)  where {I<:Integer}
    hi = floor(I, a.hi)
    lo = zero(I)

    if hi == a.hi
        lo = floor(I, a.lo)
    end

    hi + lo
end


@inline function Base.ceil(a::DoubleF64)
    hi = ceil(a.hi)
    lo = 0.0

    if hi == a.hi
        lo = ceil(a.lo)
        hi, lo = quick_two_sum(hi, lo)
    end

    DoubleF64(hi, lo)
end

@inline function Base.ceil(::Type{I}, a::DoubleF64)  where {I<:Integer}
    hi = ceil(I, a.hi)
    lo = zero(I)

    if hi == a.hi
        lo = ceil(I, a.lo)
    end

    hi + lo
end

Base.trunc(a::DoubleF64) = a.hi ≥ 0.0 ? floor(a) : ceil(a)
Base.trunc(::Type{I}, a::DoubleF64) where {I<:Integer} =
    a.hi ≥ 0.0 ? floor(I, a) : ceil(I, a)
Base.isinteger(x::DoubleF64) = iszero(x - trunc(x))

    # import Random: AbstractRNG, GLOBAL_RNG
    # function Base.rand(rng::AbstractRNG, S::Type{DoubleF64})
    #     u = rand(rng, UInt64)
    #     f = Float64(u)
    #     uf = UInt64(f)
    #     ur = uf > u ? uf - u : u - uf
    #     DoubleF64(5.421010862427522e-20 * f, 5.421010862427522e-20 * Float64(ur))
    # end
    # Base.rand(::Type{DoubleF64}) = rand(GLOBAL_RNG, DoubleF64)
    #
    # function Base.rand(rng::AbstractRNG, T::Type{<:DoubleF64}, dims::Vararg{Int, N}) where N
    #     rands = Array{T}(dims)
    #     for l in eachindex(rands)
    #         rands[l] = rand(rng, T)
    #     end
    #     rands
    # end
    # Base.rand(T::Type{<:DoubleF64}, dims::Vararg{Int, N}) where N = rand(GLOBAL_RNG, T, dims)
    #
    # Base.rand(::Type{Complex{DoubleF64}}) = rand(Complex{DoubleF64})
    # Base.rand(::Type{Complex{DoubleF64}}, dims::Vararg{Int, N}) where N = rand(Complex{DoubleF64}, dims)
    # Base.rand(rng::AbstractRNG, ::Type{Complex{DoubleF64}}) = rand(rng, Complex{DoubleF64})
    # Base.rand(rng::AbstractRNG, ::Type{Complex{DoubleF64}}, dims::Vararg{Int, N}) where N = rand(rng, Complex{DoubleF64}, dims)

function Base.decompose(a::DoubleF64)::Tuple{Int128,Int,Int}
    hi, lo = a.hi, a.lo
    num1, pow1, den1 = Base.decompose(hi)
    num2, pow2, den2 = Base.decompose(lo)

    num = Int128(num1)

    pdiff = pow1 - pow2
    shift = min(pdiff, 52)
    signed_num = den1 * (Int128(num) << shift) # den1 is +1/-1
    signed_num += den2 * (num2 >> (pdiff - shift)) # den2 is +1/-1

    num = abs(signed_num)
    den = signed_num ≥ 0 ? 1 : -1
    pow = pow1 - shift

    num, pow, den
end

# Precomputed valuee
const inv_fact = [DoubleF64(1.0 / BigFloat(factorial(Int64(k)))) for k = 3:17]
const ninv_fact = length(inv_fact)


Base.ldexp(a::DoubleF64, exp::Int) = DoubleF64(ldexp(a.hi, exp), ldexp(a.lo, exp))
"""
	mul_pwr2(a::DoubleF64, b::Float64)

`a` * `b`,  where `b` is a power of 2.
"""
mul_pwr2(a::DoubleF64, b::Float64) = DoubleF64(a.hi * b, a.lo * b)

function Base.exp(a::DoubleF64)
    # Strategy:  We first reduce the size of x by noting that
    #
    #       exp(kr + m * log(2)) = 2^m * exp(r)^k
    #
    # where m and k are integers.  By choosing m appropriately
    # we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
    # evaluated using the familiar Taylor series.  Reducing the
    # argument substantially speeds up the convergence.

    k = 512.0
    inv_k = 1.0 / k

    if a.hi ≤ -709.0
        return zero(a)
    end

    if a.hi ≥ 709.0
        convert(DoubleF64, Inf)
    end

    if iszero(a)
        return one(a)
    end

    if isone(a)
        return double_e
    end

    m = floor(a.hi / double_log2.hi + 0.5)
    r = mul_pwr2(a - double_log2 * m, inv_k)

    #
    p = square(r)
    s = r + mul_pwr2(p, 0.5)
    p *= r
    t = p * inv_fact[1]


    i = 1
    while true
        s += t
        p *= r
        i += 1
        t = p * inv_fact[i]

        if abs(convert(Float64, t)) < inv_k * double_eps || i > 6
            break
        end
    end

    s += t

    s = mul_pwr2(s, 2.0) + square(s)
    s = mul_pwr2(s, 2.0) + square(s)
    s = mul_pwr2(s, 2.0) + square(s)
    s = mul_pwr2(s, 2.0) + square(s)
    s = mul_pwr2(s, 2.0) + square(s)
    s = mul_pwr2(s, 2.0) + square(s)
    s = mul_pwr2(s, 2.0) + square(s)
    s = mul_pwr2(s, 2.0) + square(s)
    s = mul_pwr2(s, 2.0) + square(s)
    s += 1.0

    return ldexp(s, convert(Int, m))
end

Base.log10(a::DoubleF64) = log(a) / double_log10
Base.log2(a::DoubleF64) = log(a) / double_log2
function Base.log(a::DoubleF64)
     # Strategy.  The Taylor series for log converges much more
     # slowly than that of exp, due to the lack of the factorial
     # term in the denominator.  Hence this routine instead tries
     # to determine the root of the function
         #
     #     f(x) = exp(x) - a
         #
     # using Newton iteration.  The iteration is given by
         #
     #     x' = x - f(x)/f'(x)
     #        = x - (1 - a * exp(-x))
     #        = x + a * exp(-x) - 1.
         #
     # Only one iteration is needed, since Newton's iteration
     # approximately doubles the number of digits per iteration.

    if isone(a)
        return zero(DoubleF64)
    end

    if a.hi ≤ 0.0
        throw(DomainError("Passed non-positive argument to log."))
    end

    x = DoubleF64(log(a.hi)) # Initial approximation

    x + a * exp(-x) - 1.0
end


const sin_table = [DoubleF64(sin(k * big(π) * 0.0625)) for k = 1:4]
const cos_table = [DoubleF64(cos(k * big(π) * 0.0625)) for k = 1:4]


"""
	sin_taylor(a::DoubleF64)

Computes sin(a) using Taylor series. Assumes |a| <= π/32.
"""
function sin_taylor(a::DoubleF64)
    thresh = 0.5 * abs(convert(Float64, a)) * double_eps


    if iszero(a)
        return zero(a)
    end

    i = 1
    x = -square(a)
    s = a
    r = a
    while true
        r *= x
        t = r * inv_fact[i]
        s += t
        i += 2
        if i > ninv_fact || abs(convert(Float64, t)) < thresh
            break
        end
    end

    s
end

"""
	cos_taylor(a::DoubleF64)

Computes cos(a) using Taylor series. Assumes |a| <= π/32.
"""
function cos_taylor(a::DoubleF64)
    thresh = 0.5 * double_eps

    if iszero(a)
        return one(a)
    end

    x = -square(a)
    r = x
    s = 1.0 + mul_pwr2(r, 0.5)
    i = 2
    while true
        r *= x
        t = r * inv_fact[i]
        s += t
        i += 2
        if i > ninv_fact || abs(convert(Float64, t)) < thresh
            break
        end
    end

    s
end

function sincos_taylor(a::DoubleF64)
    if iszero(a)
        zero(a), one(a)
    end

    sin_a = sin_taylor(a)
    cos_a = sqrt(1.0 - square(sin_a))

    sin_a, cos_a
end

function Base.sin(a::DoubleF64)
        # Strategy.  To compute sin(x), we choose integers a, b so that
        #
        #    x = s + a * (pi/2) + b * (pi/16)
        #
        #  and |s| <= pi/32.  Using the fact that
        #
        #    sin(pi/16) = 0.5 * sqrt(2 - sqrt(2 + sqrt(2)))
        #
        #  we can compute sin(x) from sin(s), cos(s).  This greatly
        #  increases the convergence of the sine Taylor series.

    if iszero(a)
        return zero(a)
    end

    # approximately reduce modulo 2*pi
    z = round(a / double_2pi)
    r = a - double_2pi * z

    # approximately reduce modulo pi/2 and then modulo pi/16.

    q = floor(r.hi / double_pi2.hi + 0.5)
    t = r - double_pi2 * q
    j = convert(Int, q)
    q = floor(t.hi / double_pi16.hi + 0.5)
    t -= double_pi16 * q
    k = convert(Int, q)
    abs_k = abs(k)

    if j < -2 || j > 2
    # Cannot reduce modulo pi/2.
        return double_nan
    end

    if (abs_k > 4)
    # Cannot reduce modulo pi/16.
        return double_nan
    end

    if k == 0
        if j == 0
            return sin_taylor(t)
        elseif j == 1
            return cos_taylor(t)
        elseif j == -1
            return -cos_taylor(t)
        else
            return -sin_taylor(t)
        end
    end

    u = cos_table[abs_k]
    v = sin_table[abs_k]
    sin_t, cos_t = sincos_taylor(t)

    if j == 0
        r = k > 0 ? u * sin_t + v * cos_t : u * sin_t - v * cos_t
    elseif j == 1
        r = k > 0 ? u * cos_t - v * sin_t : u * cos_t + v * sin_t
    elseif j == -1
        r = k > 0 ? v * sin_t - u * cos_t : -u * cos_t - v * sin_t
    else
        r = k > 0 ? -u * sin_t - v * cos_t : v * cos_t - u * sin_t
    end

    r
end

function Base.cos(a::DoubleF64)
    if iszero(a)
        return one(a)
    end

    # approximately reduce modulo 2*pi
    z = round(a / double_2pi)
    r = a - double_2pi * z

    # approximately reduce modulo pi/2 and then modulo pi/16.

    q = floor(r.hi / double_pi2.hi + 0.5)
    t = r - double_pi2 * q
    j = convert(Int, q)
    q = floor(t.hi / double_pi16.hi + 0.5)
    t -= double_pi16 * q
    k = convert(Int, q)
    abs_k = abs(k)

    if j < -2 || j > 2
    # Cannot reduce modulo pi/2.
        return double_nan
    end

    if (abs_k > 4)
    # Cannot reduce modulo pi/16.
        return double_nan
    end

    if k == 0
        if j == 0
            return cos_taylor(t)
        elseif j == 1
            return -sin_taylor(t)
        elseif j == -1
            return sin_taylor(t)
        else
            return -cos_taylor(t)
        end
    end

    u = cos_table[abs_k]
    v = sin_table[abs_k]
    sin_t, cos_t = sincos_taylor(t)

    if j == 0
        r = k > 0 ? u * cos_t - v * sin_t : u * cos_t + v * sin_t
    elseif j == 1
        r = k > 0 ? -u * sin_t - v * cos_t : v * cos_t - u * sin_t
    elseif j == -1
        r = k > 0 ? u * sin_t + v * cos_t : u * sin_t - v * cos_t
    else
        r = k > 0 ? v * sin_t - u * cos_t : -u * cos_t - v * sin_t
    end

    r
end


function Base.sincos(a::DoubleF64)
    if iszero(a)
        return zero(a), one(a)
    end

    # approximately reduce modulo 2*pi
    z = round(a / double_2pi)
    r = a - double_2pi * z

    # approximately reduce modulo pi/2 and then modulo pi/16.
    q = floor(r.hi / double_pi2.hi + 0.5)
    t = r - double_pi2 * q
    j = convert(Int, q)
    abs_j = abs(j)
    q = floor(t.hi / double_pi16.hi + 0.5)
    t -= double_pi16 * q
    k = convert(Int, q)
    abs_k = abs(k)

    if abs_j > 2
    # Cannot reduce modulo pi/2.
        return double_nan, double_nan
    end

    if abs_k > 4
    # Cannot reduce modulo pi/16.
        return double_nan, double_nan
    end

    sin_t, cos_t = sincos_taylor(t)

    if abs_k == 0
        s = sin_t
        c = cos_t
    else
        u = cos_table[abs_k]
        v = sin_table[abs_k]

        s = k > 0 ? u * sin_t + v * cos_t : u * sin_t - v * cos_t
        c = k > 0 ? u * cos_t - v * sin_t : u * cos_t + v * sin_t
    end

    if j == 0
        s, c
    elseif j == 1
        c, -s
    elseif j == -1
        -c, s
    else
        -s, -c
    end
end

Base.atan(a::DoubleF64) = atan2(a, one(a))

function Base.atan(y::DoubleF64, x::DoubleF64)
     # Strategy: Instead of using Taylor series to compute
         # arctan, we instead use Newton's iteration to solve
         # the equation
         #
         #    sin(z) = y/r    or    cos(z) = x/r
         #
         # where r = sqrt(x^2 + y^2).
         # The iteration is given by
         #
         #    z' = z + (y - sin(z)) / cos(z)          (for equation 1)
         #    z' = z - (x - cos(z)) / sin(z)          (for equation 2)
         #
         # Here, x and y are normalized so that x^2 + y^2 = 1.
         # If |x| > |y|, then first iteration is used since the
         # denominator is larger.  Otherwise, the second is used.
    if iszero(x)
        if iszero(y)
            double_nan
        end

        return y.hi > 0.0 ? double_pi2 : -double_pi2
    elseif iszero(y)
        return x.hi > 0.0 ? zero(x) : double_pi
    end

    if x == y
        return y.hi > 0.0 ? double_pi4 : -double_3pi4
    end

    if x == -y
        return y.hi > 0.0 ? double_3pi4 : -double_pi4
    end

    r = sqrt(square(x) + square(y))
    xx = x / r
    yy = y / r

    # Compute double precision approximation to atan.
    z = DoubleF64(atan2(convert(Float64, y), convert(Float64, x)))

    if abs(xx.hi) > abs(yy.hi)
    # Use Newton iteration 1.  z' = z + (y - sin(z)) / cos(z)
        sin_z, cos_z = sincos(z)
        z += (yy - sin_z) / cos_z
    else
    # Use Newton iteration 2.  z' = z - (x - cos(z)) / sin(z)
        sin_z, cos_z = sincos(z)
        z -= (xx - cos_z) / sin_z
    end

    z
end

function Base.tan(a::DoubleF64)
    s, c = sincos(a)
    s / c
end

function Base.asin(a::DoubleF64)
    abs_a = abs(a)

    if abs_a > 1.0
        throw(DomainError())
    end

    if isone(abs_a)
        return a.hi > 0.0 ? double_pi2 : -double_pi2
    end

    atan2(a, sqrt(1.0 - square(a)))
end

function Base.acos(a::DoubleF64)
    abs_a = abs(a)

    if abs_a > 1.0
        throw(DomainError())
    end

    if isone(abs_a)
        return a.hi > 0.0 ? zero(a) : -double_pi
    end

    atan2(sqrt(1.0 - square(a)), a)
end


    # function Base.sinh(a::DoubleF64)
    # 	if iszero(a)
    # 		return zero(a)
    # 	end
    # 	if a > 0.05 || a < -0.05
    # 		ea = exp(a)
    # 		return mul_pwr2(ea - inv(ea), 0.5)
    # 	end
    #
    # 	# The computation below didn't yield a satisfying accuracy
    # 	# Thus we just fallback to the Float64 formula which yields at least
    # 	# a
    # 	return DoubleF64(sinh(convert(Float64, a)))
    #
    # 	#
    # 	#
    # 	# #=  since a is small, using the above formula gives
    #     # a lot of cancellation.  So use Taylor series. =#
    # 	# s = a
    # 	# t = a
    # 	# r = square(t)
    # 	# m = 1.0
    # 	# thresh = abs(convert(Float64, a)) * double_eps
    # 	#
    # 	# while true
    # 	# 	m += 2.0
    # 	# 	t *= r
    # 	# 	t /= (m - 1) * m
    # 	#
    # 	# 	s += t
    # 	#
    # 	# 	if abs(t) < thresh
    # 	# 		break
    # 	# 	end
    # 	# end
    # 	#
    # 	# s
    # end

function Base.cosh(a::DoubleF64)
    if iszero(a)
        return one(a)
    end
    ea = exp(a)
    mul_pwr2(ea + inv(ea), 0.5)
end
    #
    # """
    # 	sincosh(x)
    #
    # Compute `(sinh(x), cosh(x))`. This is faster than computing the values
    # separetly.
    # """
    # function sincosh(a::DoubleF64)
    # 	if abs(convert(Float64, a)) ≤ 0.05
    # 		s = sinh(a)
    # 		c = sqrt(1.0 + square(s))
    # 	else
    # 		ea = exp(a)
    # 		inv_ea = inv(ea)
    # 		s = mul_pwr2(ea - inv_ea, 0.5)
    # 		c = mul_pwr2(ea + inv_ea, 0.5)
    # 	end
    # 	s, c
    # end

function Base.tanh(a::DoubleF64)
    if iszero(a)
        return zero(a)
    end
    if abs(convert(Float64, a)) > 0.05
        ea = exp(a)
        inv_ea = inv(ea)
        return (ea - inv_ea) / (ea + inv_ea)
    else
        s = sinh(a)
        c = sqrt(1.0 + square(s))
        return s / c
    end
end

function Base.show(io::IO, x::DoubleF64)
    Printf.@printf io "%.32g" big(x)  # crude approximation to valid number of digits
end
end # module
