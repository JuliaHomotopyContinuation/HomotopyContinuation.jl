"""
```
FixedPoly
```

A structure for fast multivariate FixedPolynomial evaluation.

### Fields
* `exponents::Matrix{Int}`: Each column represents the exponent of a term. The columns are sorted by total degree.
* `coeffs::Vector{T}`: List of the coefficients.
* `vars::Vector{Symbol}`: List of variables.

### Example
```
FixedPoly: 3XYZ^2 - 2X^3Y
exponents:
    [ 3 1
      1 1
      0 2 ]
coeffs: [-2.0, 3.0]
vars: [:x, :y, :z]
```
"""
struct FixedPoly{T<:Number}
    exponents::Matrix{Int}
    coeffs::Vector{T}
    vars::Vector{Symbol}

    function FixedPoly{T}(exponents::Matrix{Int}, coeffs::Vector{T}, vars::Vector{Symbol}) where {T<:Number}
        sorted_cols = sort!([1:size(exponents,2);], lt=((i, j) -> lt_total_degree(exponents[:,i], exponents[:,j])), rev=true)
        if length(vars) != size(exponents, 1)
             error("Number of vars and rows of exponents doesn't match")
         end
        new(exponents[:, sorted_cols], coeffs[sorted_cols], vars)
    end
end
FixedPoly{T}(exponents::Matrix{Int}, coeffs::Vector{T}, vars::Vector{Symbol})=FixedPoly{T}(exponents, coeffs, vars)

function FixedPoly{T,N}(terms::Dict{NTuple{N,Int}, T}, vars::NTuple{N,Symbol})
    term_keys = terms |> keys |> collect
    exponents = hcat(map(collect, term_keys)...)
    coeffs = map(key -> terms[key], term_keys)

    FixedPoly(exponents, coeffs, collect(vars))
end
FixedPoly(p::Poly) = FixedPoly(terms(p), vars(p))

system(p::Vector{Poly{T,N}}) where {T,N} = FixedPoly.(p)
system(p...) = system(collect(p))

*(a, p::FixedPoly) =  FixedPoly(p.exponents, a .* p.coeffs, p.vars)
*(p::FixedPoly, a) =  FixedPoly(p.exponents, a .* p.coeffs, p.vars)
# .*{T<:Number}(a::Number, P::AbstractArray{FixedPoly{T}}) =  map(p -> a * p, P)
# .*{T<:Number}(P::AbstractArray{FixedPoly{T}}, a::Number) =  map(p -> p * a, P)


# helpers
function lt_total_degree(a::Vector{T}, b::Vector{T}) where {T<:Real}
    sum_a = sum(a)
    sum_b = sum(b)
    if sum_a < sum_b
        return true
    elseif sum_a > sum_b
        return false
    else
        for i in eachindex(a)
            if a[i] < b[i]
                return true
            elseif a[i] > b[i]
                return false
            end
        end
    end
    false
end

Base.eltype(p::FixedPoly{T}) where {T} = T

"""
    exponents(p::FixedPoly)

Returns the exponents matrix
"""
@inline exponents(p::FixedPoly) = p.exponents

"""
    coeffs(p::FixedPoly)

Returns the coefficient vector
"""
@inline coeffs(p::FixedPoly) = p.coeffs

"""
    vars(p::FixedPoly)

Returns the variables vector
"""
@inline vars(p::FixedPoly) = p.vars


"""
    nterms(p::FixedPoly)

Returns the number of terms of p
"""
@inline nterms(p::FixedPoly) = size(exponents(p), 2)

"""
    nvars(p::FixedPoly)

Returns the number of vars of p
"""
@inline nvars(p::FixedPoly) = length(vars(p))

"""
    deg(p::FixedPoly)

Returns the degree of p
"""
@inline deg(p::FixedPoly) = sum(exponents(p)[:,1])


# ITERATOR
start(p::FixedPoly) = (1, nterms(p))
function next(p::FixedPoly, state::Tuple{Int,Int})
    (i, limit) = state
    newstate = (i + 1, limit)
    val = (coeffs(p)[i], exponents(p)[:,i])

    (val, newstate)
end
done(p::FixedPoly, state) = state[1] > state[2]
length(p::FixedPoly) = nterms(p)

function getindex(f::FixedPoly, I)
    A = f.exponents
    res = 1:size(A,2)
    for i in eachindex(x)
        res = res[find(y -> y == x[i], A[i,res])]
        if isempty(res)
            return zero(eltype(f))
        end
    end

    f.coeffs[res[1]]
end

"""
    evaluate(p::FixedPoly{T}, x::Vector{T})

Evaluates `p` at `x`, i.e. p(x)
"""
function evaluate(p::FixedPoly{T}, x::Vector{T})::T where {T<:Number}
    cfs = coeffs(p)
    exps = exponents(p)
    nvars, nterms = size(exps)
    res = zero(T)
    for j = 1:nterms
        term = cfs[j]
        for i = 1:nvars
            k = exps[i, j]
            term *= x[i]^k
        end
        res += term
    end
    res
end
(p::FixedPoly{T})(x::Vector{T}) where {T<:Number} = evaluate(p, x)


"""
    gradient(p::FixedPoly)

Differentiates FixedPoly `p`. Returns the gradient vector.
"""
function gradient(p::FixedPoly)
    p_exps = exponents(p)
    p_coeffs = coeffs(p)
    p_vars = vars(p)
    n_vars, n_terms = size(p_exps)

    grad = Vector{FixedPoly{eltype(p)}}()
    for i_var=1:n_vars
        exps = copy(p_exps)
        cfs = copy(p_coeffs)
        for j=1:n_terms
            k = exps[i_var, j]
            if k > 0
                exps[i_var, j] = max(0, k - 1)
                cfs[j] *= k
            else
                exps[:,j] = zeros(Int, n_vars, 1)
                cfs[j] = zero(eltype(p))
            end
        end
        push!(grad, FixedPoly(exps, cfs, copy(p_vars)))
    end
    grad
end


"""
    is_homogenous(p::FixedPoly)

homogenizes `p`
"""
function is_homogenous(p::FixedPoly)
    monomials_degree = sum(exponents(p), 1)
    max_deg = monomials_degree[1]
    all(x -> x == max_deg, monomials_degree)
end

"""
    homogenize(p::FixedPoly)

homogenizes `p`
"""
function homogenize(p::FixedPoly, var::Symbol=:⬣)
    monomials_degree = sum(exponents(p), 1)
    max_deg = monomials_degree[1]

    FixedPoly([max_deg - monomials_degree; exponents(p)], coeffs(p), [var; vars(p)])
end

"""
    dehomogenize(p::FixedPoly)

dehomogenizes `p`
"""
dehomogenize(p::FixedPoly) = FixedPoly(exponents(p)[2:end,:], coeffs(p), vars(p)[2:end])

"""
    multinomial(k::Vector{Int})

Computes the multinomial coefficient (|k| \\over k)
"""
function multinomial(k::Vector{Int})
    s = 0
    result = 1
    @inbounds for i in k
        s += i
        result *= binomial(s, i)
    end
    result
end

"""
    weyl_dot(f [, g])

computes the Bombieri-Weyl dot product between `FixedPoly`s `f` and `g`. If `g` is omitted, ``<f,f>`` is calculated.
Assumes that `f` and `g` are homogenous. See [here](https://en.wikipedia.org/wiki/Bombieri_norm)
for more details.
"""
function weyl_dot{T<:Complex}(f::FixedPoly{T},g::FixedPoly{T})
    if (f === g)
        return sum(x -> abs2(x[1]) / multinomial(x[2]), f)
    end
    result = zero(T)
    for (c_f, exp_f) in f
        normalizer = multinomial(exp_f)
        for (c_g, exp_g) in g
            if exp_f == exp_g
                result += (c_f * conj(c_g)) / normalizer
                break
            end
        end
    end
    result
end

"""
    weyl_norm(f::FixedPoly)

computes the Bombieri-Weyl norm for `f`. Assumes that `f` is homogenous.
See [here](https://en.wikipedia.org/wiki/Bombieri_norm) for more details.
"""
weyl_norm{T<:Complex}(f::FixedPoly{T}) = √weyl_dot(f,f)
