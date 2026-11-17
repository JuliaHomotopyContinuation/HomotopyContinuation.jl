# Accumulate a short sum of double-double products before normalizing. 
# Complex perations otherwise normalize each real product and every subsequent sum.
@inline function _dd_product_sum(
    a::NTuple{N,DoubleF64},
    b::NTuple{N,DoubleF64},
    c::Union{Nothing,DoubleF64} = nothing,
) where {N}
    s, e = two_prod(a[1].hi, b[1].hi)
    e += muladd(a[1].hi, b[1].lo, a[1].lo * b[1].hi)
    @inbounds for k = 2:N
        p, pe = two_prod(a[k].hi, b[k].hi)
        s, se = two_sum(s, p)
        e += se + pe + muladd(a[k].hi, b[k].lo, a[k].lo * b[k].hi)
    end
    if c !== nothing
        s, se = two_sum(s, c.hi)
        e += se + c.lo
    end
    # Cancellation can make |e| > |s|, so quick_two_sum is not sufficient here.
    hi, lo = two_sum(s, e)
    DoubleF64(hi, lo)
end

# Fused complex double-double operations. Real and imaginary part are each computed
# as one sum of products with `_dd_product_sum`:
# * `op_muladd(x, y, z) = x*y + z`
# * `op_mulsub(x, y, z) = x*y - z`
# * `op_submul(x, y, z) = z - x*y`
# * `op_mulmuladd(x, y, z, w) = x*y + z*w`
# * `op_mulmulsub(x, y, z, w) = x*y - z*w`
@inline function op_muladd(x::ComplexDF64, y::ComplexDF64, z::ComplexDF64)
    a, b = reim(x)
    c, d = reim(y)
    Complex(
        _dd_product_sum((a, -b), (c, d), real(z)),
        _dd_product_sum((a, b), (d, c), imag(z)),
    )
end

@inline op_mulsub(x::ComplexDF64, y::ComplexDF64, z::ComplexDF64) = op_muladd(x, y, -z)
@inline op_submul(x::ComplexDF64, y::ComplexDF64, z::ComplexDF64) = op_muladd(-x, y, z)

@inline function op_mulmuladd(
    x::ComplexDF64,
    y::ComplexDF64,
    z::ComplexDF64,
    w::ComplexDF64,
)
    a, b = reim(x)
    c, d = reim(y)
    e, f = reim(z)
    g, h = reim(w)
    Complex(
        _dd_product_sum((a, -b, e, -f), (c, d, g, h)),
        _dd_product_sum((a, b, e, f), (d, c, h, g)),
    )
end

@inline op_mulmulsub(x::ComplexDF64, y::ComplexDF64, z::ComplexDF64, w::ComplexDF64) =
    op_mulmuladd(x, y, -z, w)
