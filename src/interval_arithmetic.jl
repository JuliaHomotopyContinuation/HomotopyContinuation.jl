module IntervalArithmetic

export Interval,
    IComplex, IComplexF64, mid, diam, rad, mig, mag, hull, isinterior, isdisjoint

import Base: *, /, +, -, ^

@static if VERSION ≥ v"1.5.0-"
    import Base: isdisjoint
end

import Printf

# Implementation of interval arithmetic following "Interval Analysis" - Mayer.

struct Interval{T<:Real}
    lo::T
    hi::T
end

Interval{T}(a::Real) where {T<:Real} = Interval(convert(T, a))
Interval{T}(a::Interval{S}) where {T<:Real,S<:Real} =
    Interval(convert(T, a.lo), convert(T, a.hi))
Interval(a::T) where {T<:Real} = Interval(a, a)
Interval(a::T, b::S) where {T<:Real,S<:Real} = Interval(promote(a, b)...)
Interval(a::T, b::T) where {T<:Integer} = Interval(float(a), float(b))

Base.convert(::Type{Interval{T}}, a::Interval{S}) where {T,S} =
    Interval(convert(T, a.lo), convert(T, a.hi))
Base.convert(::Type{Interval{T}}, a::Real) where {T} = Interval{T}(a)

function interval(a::Real, b::Real)
    is_valid_interval(a, b) || invalid_interval_error(a, b)
    Interval(a, b)
end
is_valid_interval(a::Real, b::Real) = isfinite(a) && isfinite(b) && a ≤ b
@noinline function invalid_interval_error(a, b)
    throw(ArgumentError("`[$a, $b]` is not a valid interval. Need `a ≤ b` to construct `interval(a, b)`."))
end

Base.hash(x::Interval, h::UInt) = hash(x.hi, hash(x.lo, u))
Base.:(==)(a::Interval, b::Interval) = a.lo == b.lo && a.hi == b.hi
Base.eltype(x::Interval{T}) where {T} = T

Base.promote_rule(::Type{Interval{T}}, ::Type{Interval{S}}) where {T<:Real,S<:Real} =
    Interval{promote_type(T, S)}
Base.promote_rule(::Type{Interval{T}}, ::Type{S}) where {T<:Real,S<:Real} =
    Interval{promote_type(T, S)}

round_up(a::AbstractFloat) = nextfloat(a)
round_down(a::AbstractFloat) = prevfloat(a)
macro round(a, b)
    :(Interval(round_down($(esc(a))), round_up($(esc(b)))))
end

mid(a::Interval) = (a.lo + a.hi) / 2
diam(a::Interval) = round_up(a.hi - a.lo)
function rad(a::Interval)
    m = mid(a)
    return round_up(max(m - a.lo, a.hi - m))
end
function mag(a::Interval{T}) where {T}
    max(abs(a.lo), abs(a.hi))
end
function mig(a::Interval{T}) where {T}
    zero(T) ∈ a && return zero(T)
    min(abs(a.lo), abs(a.hi))
end
hull(a::Interval, b::Interval) = Interval(min(a.lo, b.lo), max(a.hi, b.hi))

Base.issubset(a::Interval, b::Interval) = (a.lo ≥ b.lo) && (a.hi ≤ b.hi)
isinterior(a::Interval, b::Interval) = (a.lo > b.lo) && (a.hi < b.hi)
isdisjoint(a::Interval, b::Interval) = (b.hi < a.lo) || (a.hi < b.lo)
function Base.intersect(a::Interval, b::Interval)
    c = Interval(max(a.lo, b.lo), min(a.hi, b.hi))
    c.lo ≤ c.hi || return Interval(convert(eltype(c), NaN))
    c
end
Base.isempty(a::Interval) = isnan(a.lo) || isnan(a.hi)
Base.in(x::Number, a::Interval) = a.lo ≤ x ≤ a.hi

function Base.show(io::IO, a::Interval)
    print(io, mid(a), " ± ")
    Printf.@printf(io, "%.5g", rad(a))
end

# arithmetic

# zero, one, typemin, typemax
Base.zero(a::Interval{T}) where {T} = Interval(zero(T))
Base.zero(::Type{Interval{T}}) where {T} = Interval(zero(T))
Base.one(a::Interval{T}) where {T} = Interval(one(T))
Base.one(::Type{Interval{T}}) where {T} = Interval(one(T))

# The arithmetic of Interval is based on the IntervalArithmetic.jl package
# https://github.com/JuliaIntervals/IntervalArithmetic.jl/
# IntervalArithmetic.jl is licensed under the MIT "Expat" License

## Addition and subtraction
+(a::Interval) = a
+(a::Interval, b) = @round(a.lo + b, a.hi + b)
+(b, a::Interval) = a + b
function +(a::Interval, b::Interval)
    @round(a.lo + b.lo, a.hi + b.hi)
end

-(a::Interval) = Interval(-a.hi, -a.lo)
-(a::Interval, b) = @round(a.lo - b, a.hi - b)
-(b, a::Interval) = @round(b - a.hi, b - a.lo)
-(a::Interval, b::Interval) = @round(a.lo - b.hi, a.hi - b.lo)


## Multiplication
function *(x, a::Interval)
    (iszero(a) || iszero(x)) && return zero(a.lo * x)
    if x ≥ 0.0
        return @round(a.lo * x, a.hi * x)
    else
        return @round(a.hi * x, a.lo * x)
    end
end
*(a::Interval, x) = x * a
function *(a::Interval, b::Interval)
    if b.lo >= zero(b.lo)
        a.lo >= zero(a.lo) && return @round(*(a.lo, b.lo), *(a.hi, b.hi))
        a.hi <= zero(a.hi) && return @round(*(a.lo, b.hi), *(a.hi, b.lo))
        return @round(a.lo * b.hi, a.hi * b.hi)   # zero(T) ∈ a
    elseif b.hi <= zero(b.hi)
        a.lo >= zero(a.lo) && return @round(*(a.hi, b.lo), *(a.lo, b.hi))
        a.hi <= zero(a.hi) && return @round(*(a.hi, b.hi), *(a.lo, b.lo))
        return @round(a.hi * b.lo, a.lo * b.lo)   # zero(T) ∈ a
    else
        a.lo > zero(a.lo) && return @round(*(a.hi, b.lo), *(a.hi, b.hi))
        a.hi < zero(a.hi) && return @round(*(a.lo, b.hi), *(a.lo, b.lo))
        return @round(min(*(a.lo, b.hi), *(a.hi, b.lo)), max(*(a.lo, b.lo), *(a.hi, b.hi)))
    end
end

function Base.muladd(a::Interval, b::Interval, c::Interval)
    lo = let
        lo1 = muladd(a.lo, b.lo, c.lo)
        lo2 = muladd(a.lo, b.hi, c.lo)
        lo3 = muladd(a.hi, b.lo, c.lo)
        lo4 = muladd(a.hi, b.hi, c.lo)
        round_down(Base.FastMath.min_fast(lo1, lo2, lo3, lo4))
    end

    hi = let
        hi1 = muladd(a.lo, b.lo, c.hi)
        hi2 = muladd(a.lo, b.hi, c.hi)
        hi3 = muladd(a.hi, b.lo, c.hi)
        hi4 = muladd(a.hi, b.hi, c.hi)
        round_up(Base.FastMath.max_fast(hi1, hi2, hi3, hi4))
    end

    Interval(lo, hi)
end
function Base.muladd(a::Interval{T}, b::T, c::Interval) where {T}
    lo = let
        lo1 = muladd(a.lo, b, c.lo)
        lo2 = muladd(a.hi, b, c.lo)
        round_down(Base.FastMath.min_fast(lo1, lo2))
    end

    hi = let
        hi1 = muladd(a.lo, b, c.hi)
        hi2 = muladd(a.hi, b, c.hi)
        round_up(Base.FastMath.max_fast(hi1, hi2))
    end

    Interval(lo, hi)
end
function Base.muladd(a::T, b::Interval{T}, c::Interval) where {T}
    lo = let
        lo1 = muladd(a, b.lo, c.lo)
        lo2 = muladd(a, b.hi, c.lo)
        round_down(Base.FastMath.min_fast(lo1, lo2))
    end

    hi = let
        hi1 = muladd(a, b.lo, c.hi)
        hi2 = muladd(a, b.hi, c.hi)
        round_up(Base.FastMath.max_fast(hi1, hi2))
    end

    Interval(lo, hi)
end

function inv(a::Interval{T}) where {T<:Real}
    if zero(T) ∈ a
        return Interval(convert(T, NaN))
    end
    @round(inv(a.hi), inv(a.lo))
end

function /(a::Interval{T}, x::T) where {T<:Real}
    iszero(a) && return zero(Interval{T})

    if x ≥ 0.0
        return @round(a.lo / x, a.hi / x)
    else
        return @round(a.hi / x, a.lo / x)
    end
end
function /(a::Interval{T}, b::Interval{T}) where {T<:Real}
    if iszero(a)
        S = typeof(a.lo / b.lo)
        return zero(Interval{S})
    end

    # b contains zero
    if b.lo > zero(T) # b strictly positive
        a.lo >= zero(T) && return @round(a.lo / b.hi, a.hi / b.lo)
        a.hi <= zero(T) && return @round(a.lo / b.lo, a.hi / b.hi)
        return @round(a.lo / b.lo, a.hi / b.lo)  # zero(T) ∈ a

    elseif b.hi < zero(T) # b strictly negative

        a.lo >= zero(T) && return @round(a.hi / b.hi, a.lo / b.lo)
        a.hi <= zero(T) && return @round(a.hi / b.lo, a.lo / b.hi)
        return @round(a.hi / b.hi, a.lo / b.hi)  # zero(T) ∈ a
    else   # b contains zero, but is not zero(b)
        S = typeof(a.lo / b.lo)
        return Interval(convert(S, NaN))
    end
end

function sqr(a::Interval{T}) where {T<:Real}
    if a.lo ≥ zero(T)
        return @round(a.lo^2, a.hi^2)

    elseif a.hi ≤ zero(T)
        return @round(a.hi^2, a.lo^2)
    end

    return @round(mig(a)^2, mag(a)^2)
end
Base.literal_pow(::typeof(^), a::Interval, ::Val{2}) = sqr(a)

function ^(x::Interval, n::Integer)  # fast integer power
    if n < 0
        return inv(pow(x, -n))
    end
    isempty(x) && return x
    if iseven(n) && 0 ∈ x
        return hull(
            zero(x),
            hull(
                Base.power_by_squaring(Interval(mig(x)), n),
                Base.power_by_squaring(Interval(mag(x)), n),
            ),
        )

    else
        return hull(
            Base.power_by_squaring(Interval(x.lo), n),
            Base.power_by_squaring(Interval(x.hi), n),
        )
    end
end


struct IComplex{T} <: Number
    re::Interval{T}
    im::Interval{T}
end
IComplex(x::Real, y::Interval) = IComplex(promote(x, y)...)
IComplex(x::Interval, y::Real) = IComplex(promote(x, y)...)
function IComplex(x::Real, y::Real)
    ix, iy = promote(x, y)
    IComplex(Interval(ix), Interval(iy))
end
IComplex(x::Union{Interval,Real}) = IComplex(x, zero(x))
IComplex(z::Complex) = IComplex(real(z), imag(z))
IComplex(z::IComplex) = z
Base.complex(x::Interval, y::Real) = IComplex(x, y)
Base.complex(x::Real, y::Interval) = IComplex(x, y)
Base.complex(x::Interval, y::Interval) = IComplex(x, y)
Base.complex(x::IComplex) = x

const IComplexF64 = IComplex{Float64}
IComplex{T}(x::Complex{S}) where {T<:Real,S<:Real} =
    IComplex(Interval{T}(real(x)), Interval{T}(imag(x)))
IComplex{T}(x::IComplex{S}) where {T<:Real,S<:Real} =
    IComplex(Interval{T}(real(x)), Interval{T}(imag(x)))
IComplex{T}(x) where {T<:Real} = IComplex(Interval{T}(x))
# (::Type{IComplex{T}})(z::IComplex{T}) where {T<:Real} = z

Base.zero(a::IComplex{T}) where {T} = IComplex(zero(Interval{T}))
Base.zero(::Type{IComplex{T}}) where {T} = IComplex(zero(Interval{T}))
Base.one(a::IComplex{T}) where {T} = IComplex(one(Interval{T}))
Base.one(::Type{IComplex{T}}) where {T} = IComplex(one(Interval{T}))

Base.show(io::IO, c::IComplex) = print(io, "(", c.re, ") + (", c.im, ")im")
Base.Base.broadcastable(z::IComplex) = z
Base.promote_rule(::Type{IComplex{T}}, ::Type{S}) where {T,S<:Real} =
    IComplex{promote_type(T, S)}
Base.promote_rule(::Type{IComplex{T}}, ::Type{Complex{S}}) where {T,S} =
    IComplex{promote_type(T, S)}
Base.promote_rule(::Type{IComplex{T}}, ::Type{Interval{S}}) where {T,S<:Real} =
    IComplex{promote_type(T, S)}
Base.promote_rule(::Type{IComplex{T}}, ::Type{IComplex{S}}) where {T,S} =
    IComplex{promote_type(T, S)}

Base.widen(::Type{IComplex{T}}) where {T} = IComplex{widen(T)}

Base.real(z::IComplex) = z.re
Base.real(C::Type{IComplex{T}}) where {T} = Interval{T}

Base.imag(z::IComplex) = z.im
Base.imag(C::Type{IComplex{T}}) where {T} = Interval{T}
Base.reim(z::IComplex) = (real(z), imag(z))

Base.conj(z::IComplex) = IComplex(real(z), -imag(z))

+(z::IComplex) = z
+(x::Union{Interval,Real}, z::IComplex) = IComplex(x + real(z), imag(z))
+(z::IComplex, x::Union{Interval,Real}) = IComplex(x + real(z), imag(z))

+(z::IComplex, w::IComplex) = IComplex(real(z) + real(w), imag(z) + imag(w))
-(z::IComplex) = IComplex(-real(z), -imag(z))
-(z::IComplex, w::IComplex) = IComplex(real(z) - real(w), imag(z) - imag(w))
-(x::Union{Interval,Real}, z::IComplex) = IComplex(x - real(z), -imag(z))
-(z::IComplex, x::Union{Interval,Real}) = IComplex(real(z) - x, imag(z))

*(z::IComplex, w::IComplex) =
    IComplex(real(z) * real(w) - imag(z) * imag(w), real(z) * imag(w) + imag(z) * real(w))
*(x::Union{Interval,Real}, z::IComplex) = IComplex(x * real(z), x * imag(z))
*(z::IComplex, x::Union{Interval,Real}) = IComplex(x * real(z), x * imag(z))

Base.muladd(z::Union{Complex,IComplex}, w::Union{Complex,IComplex}, x::IComplex) = IComplex(
    muladd(real(z), real(w), real(x)) - imag(z) * imag(w),
    muladd(real(z), imag(w), muladd(imag(z), real(w), imag(x))),
)

/(a::R, z::S) where {R<:Real,S<:IComplex} = (T = promote_type(R, S); a * inv(T(z)))
/(z::IComplex, x::Union{Interval,Real}) = IComplex(real(z) / x, imag(z) / x)
function /(a::IComplex{T}, b::IComplex{T}) where {T}
    are, aim = reim(a)
    bre, bim = reim(b)
    denom = bre^2 + bim^2
    IComplex((are * bre + aim * bim) / denom, (aim * bre - are * bim) / denom)
end

mid(z::IComplex) = Complex(mid(real(z)), mid(imag(z)))
diam(z::IComplex) = max(diam(real(z)), diam(imag(z)))
rad(z::IComplex) = max(rad(real(z)), rad(imag(z)))
isinterior(a::IComplex, b::IComplex) =
    isinterior(real(a), real(b)) && isinterior(imag(a), imag(b))
Base.issubset(a::IComplex, b::IComplex) = real(a) ⊆ real(b) && imag(a) ⊆ imag(b)
isdisjoint(a::IComplex, b::IComplex) =
    isdisjoint(real(a), real(b)) || isdisjoint(imag(a), imag(b))
Base.intersect(a::IComplex, b::IComplex) =
    IComplex(intersect(real(a), real(b)), intersect(imag(a), imag(b)))
Base.in(x::Number, a::IComplex) = in(real(x), real(a)) && in(imag(x), imag(a))

end #module
