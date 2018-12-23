import DoubleFloats

"""
    multiply_add(z, w, x)

Compute `z * w + x`. The precision with which is computed is determined by the highest of all inputs.
"""
multiply_add(z::ComplexF64, w::ComplexF64, x::ComplexF64) = muladd(z, w, x)
function multiply_add(z::ComplexF64, w::ComplexF64, x::Complex{Double64})
    re = (DoubleFloats.mul2(real(z), real(w)) ⊕ real(x)) ⊖ DoubleFloats.mul2(imag(z), imag(w))
    im = DoubleFloats.mul2(real(z), imag(w)) ⊕ DoubleFloats.mul2(imag(z), real(w)) ⊕ imag(x)
    Complex{Double64}(re, im)
end

# These are a workaround for the bad performance in DoubleFloats due to the NaN handling
# See also https://github.com/JuliaMath/DoubleFloats.jl/issues/34
⊕(x, y) = add(x, y)
add(x, y) = x + y
add(x::Float64, y::Double64) = DoubleFloats.add_fpdb_db_nonfinite(x, y)
add(x::Double64, y::Float64) = DoubleFloats.add_dbfp_db_nonfinite(x, y)
add(x::Double64, y::Double64) = DoubleFloats.add_dbdb_db_nonfinite(x, y)
add(x::Complex, y::Complex) = Complex(add(real(x),real(y)), add(imag(x),imag(y)))

sub(x, y) = x - y
sub(x::Float64, y::Double64) = DoubleFloats.sub_fpdb_db_nonfinite(x, y)
sub(x::Double64, y::Float64) = DoubleFloats.sub_dbfp_db_nonfinite(x, y)
sub(x::Double64, y::Double64) = DoubleFloats.sub_dbdb_db_nonfinite(x, y)
sub(x::Complex, y::Complex) = Complex(sub(real(x),real(y)), sub(imag(x),imag(y)))

mul(x, y) = x * y
mul(x::Float64, y::Double64) = DoubleFloats.mul_fpdb_db_nonfinite(x, y)
mul(x::Double64, y::Float64) = DoubleFloats.mul_dbfp_db_nonfinite(x, y)
mul(x::Double64, y::Double64) = DoubleFloats.mul_dbdb_db_nonfinite(x, y)
