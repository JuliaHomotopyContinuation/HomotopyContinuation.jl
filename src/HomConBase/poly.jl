

@inline evaluate(f::MP.AbstractPolynomial, x) = f(MP.variables(f)=>x)
@inline gradient(f::MP.AbstractPolynomial) = MP.differentiate(f, MP.variables(f))
is_homogenous(f::MP.AbstractPolynomial) = MP.mindeg(f) == MP.maxdeg(f)

function homogenize(f::MP.AbstractPolynomial, variable::MP.AbstractVariable)
    deg = MP.maxdeg(f)
    terms = MP.terms(f)
    promoted_type = promote_type(typeof(variable), eltype(terms))
    newterms = map(term -> begin
        homvar_deg = deg - MP.deg(term)
        if homvar_deg == 0
            convert(promoted_type, term)
       else
            convert(promoted_type, term * variable^homvar_deg)
        end
    end, terms)
    MP.polynomial(newterms)
end

@inline deg(f::MP.AbstractPolynomial) = MP.maxdeg(f)
@inline nvars(f::MP.AbstractPolynomial) = MP.nvariables(f)

exponents(f::MP.AbstractPolynomial) = MP.exponents.(MP.terms(f))
coefficients(f::MP.AbstractPolynomial) = MP.coefficient.(MP.terms(f))

"""
    weyl_dot(f [, g])

computes the Bombieri-Weyl dot product between `FixedPoly`s `f` and `g`. If `g` is omitted, ``<f,f>`` is calculated.
Assumes that `f` and `g` are homogenous. See [here](https://en.wikipedia.org/wiki/Bombieri_norm)
for more details.
"""
function weyl_dot(f::MP.AbstractPolynomial{T}, g::MP.AbstractPolynomial{T}) where {T<:Complex}
    # if we have the identical argument we can simplify the computation
    # Maybe even just use `==`?
    if (f === g)
        return sum(x -> abs2(x[1]) / _multinomial(x[2]), zip(coefficients(f), exponents(f)))
    end
    result = zero(T)
    for (c_f, exp_f) in zip(coefficients(f), exponents(f))
        normalizer = _multinomial(exp_f)
        for (c_g, exp_g) in zip(coefficients(g), exponents(g))
            if exp_f == exp_g
                result += (c_f * conj(c_g)) / normalizer
                break
            end
        end
    end
    result
end

"""
    multinomial(k::Vector{Int})

Computes the multinomial coefficient (|k| \\over k)
"""
function _multinomial(iter)
    s = 0
    result = 1
    @inbounds for i in iter
        s += i
        result *= binomial(s, i)
    end
    result
end
# createpoly(exponents, coeffs, variables) = MPoly.FixedPoly(exponents, coeffs, variables)