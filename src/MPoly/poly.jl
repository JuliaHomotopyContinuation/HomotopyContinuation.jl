"""
```
Poly{T<:Number,N}
```

### Fields
* `terms::Dict{NTuple{N,Int}, T}`: Dict where keys are exponents and values are the coefficients.
* `vars:NTuple{N,Symbol}`: List of the variables.
"""
struct Poly{T<:Number,N}
    terms::Dict{NTuple{N,Int},T}
    vars::NTuple{N,Symbol}
end

Base.eltype(p::Poly{T,N}) where {T,N} = T

"""
    vars(p::Poly)

Returns a list of the variables
"""
vars(p::Poly) = p.vars


Terms{N,T} = Dict{NTuple{N,Int},T}
# Primitives
"""
    zeropoly(::Type{T<:Number}, vars::Vararg{Symbol,N})

Creates a poly representing 0 in the ring T[vars]
"""
zeropoly(::Type{T}, vars::NTuple{N,Symbol}) where {T<:Number,N} = Poly(Dict{NTuple{N,Int},T}(), vars)
zeropoly(::Type{T}, vars::Vararg{Symbol,N}) where {T<:Number,N} = Poly(Dict{NTuple{N,Int},T}(), vars)
zero(p::Poly{T,N}) where {T<:Number,N} = zeropoly(T, vars(p))
"""
    zeropolys(::Type{T<:Number}, vars::Vararg{Symbol,N}, dim::Int...)

Creates an Array of polys representing 0 in the ring T[vars]
"""
function zeropolys(::Type{T}, vars::NTuple{N,Symbol}, dims::Int...) where {T<:Number,N}
    reshape([zeropoly(T, vars) for _ in 1:prod(dims)], dims)
end
zeros(p::Poly{T,N}, dims::Int...) where {T<:Number,N}= zeropolys(T, vars(p), dims...)

"""
    onepoly(::Type{T<:Number}, vars::Vararg{Symbol,N})

Creates poly representing 1 in the ring T[vars]
"""
onepoly(::Type{T}, vars::NTuple{N,Symbol}) where {T<:Number,N} = Poly(Dict(ntuple(_ -> 0, N) => one(T)), vars)
onepoly(::Type{T}, vars::Vararg{Symbol,N}) where {T<:Number,N} = Poly(Dict(ntuple(_ -> 0, N) => one(T)), vars)

one(p::Poly{T,N}) where {T<:Number,N} = onepoly(T, vars(p))


"""

    terms(p::Poly)

Returns a dict of the exponents / coefficients of p
"""
terms(p::Poly) = p.terms

"""
    numvars(p::Poly)

Returns the number of variables
"""
numvars(p::Poly{T,N}) where {T,N} = N

"""
    deg(p::Poly)

Returns the degree of the polynomial. This is a O(number of terms) operation.
"""
deg(p::Poly) = maximum(map(sum, keys(p.terms)))


"""
    clean!(p::Poly)

Removes all "zero" coefficients c with norm(c) < 1e-14
"""
function clean!(p::Poly)
    for (key, val) in p.terms
        if norm(val) < 1e-14
            delete!(p.terms, key)
        end
    end
    return p
end

# ITERATOR
start(p::Poly) = start(p.terms)
function next(p::Poly, state::Int)
    ((power, coeff), state) = next(p.terms, state)

    ((coeff, power), state)
end
done(p::Poly, state) = done(p.terms, state)

# INDEXING
getindex(p::Poly{T,N}, I::Vararg{Int,N}) where {T<:Number,N} = try p.terms[I] catch zero(T) end
getindex(p::Poly{T,N}, I::NTuple{N,Int}) where {T<:Number,N} = try p.terms[I] catch zero(T) end
setindex!(p::Poly{T,N}, v::T, I::Vararg{Int,N}) where {T<:Number,N} = p.terms[I] = v
setindex!(p::Poly{T,N}, v::T, I::NTuple{N,Int}) where {T<:Number,N} = p.terms[I] = v
setindex!(p::Poly{T,N}, v::T, I::Vector{Int}) where {T<:Number,N} = setindex!(p, v, I...)

# # Copying
Base.deepcopy(p::Poly) = Poly(deepcopy(p.terms), vars(p))


"""
    generator(::Type{T<:Number}, var::Symbol)

Returns the generator of the ring T[var], i.e. X for T[X]
"""
function generator(::Type{T}, var::Symbol) where {T<:Number}
    p = zeropoly(T, var)
    p[1] = one(T)
    p
end

"""
    generators(::Type{T<:Number}, vars::Symbol...)

Returns the generators of the ring T[vars], i.e. (X, Y) for T[X,Y]
"""
function generators(::Type{T}, vars::Symbol...) where {T<:Number}
    N = length(vars)
    map(vars) do var
        p = zeropoly(T, vars)
        var_idx = findfirst(v -> v == var, p.vars)
        p[ntuple(i -> i == var_idx ? 1 : 0, N)] = one(T)
        p
    end
end


"""
    evaluate(p::Poly{T}, x::Vector{T})

evaluates p at x
"""
@inline function evaluate(p::Poly{T,N}, x::Vector{T}) where {T<:Number,N}
    value = zero(T)
    for (coeff, exponent) in p
        term = coeff
        for i in eachindex(exponent)
            term *= x[i]^exponent[i]
        end
        value += term
    end
    value
end
@inline function evaluate(p::Poly{T,N}, xs::Vararg{T,N}) where {T<:Number,N}
    value = zero(T)
    for (coeff, exponent) in p
        term = coeff
        for i in eachindex(exponent)
            term *= xs[i]^exponent[i]
        end
        value += term
    end
    value
end

(p::Poly{T,N})(x::Vararg{T,N}) where {T<:Number,N} = evaluate(p, collect(x))
(p::Poly{T,N})(x::Vector{T}) where {T<:Number,N} = evaluate(p, x)

@inline evaluate(ps::AbstractArray{Poly{T,N}}, x::Vector{T}) where {T,N} = map(p -> evaluate(p, x), ps)
@inline function evaluate(ps::AbstractArray{Poly{T,N}}, x::Vararg{T,N}) where {T,N}
    x_ = collect(x)
    map(p -> evaluate(p, x_), ps)
end

"""
    eval_var(p::Poly{T}, var::Symbol, x::T)

Evaluates the variable `var` with `x` for `p`.
"""
function eval_var(p::Poly{T,N}, var::Symbol, value::T) where {T,N}
    elim_var_idx = findfirst(v -> v == var, p.vars)
    if elim_var_idx == 0
        return p
    end
    new_vars = deleteat!(collect(p.vars), elim_var_idx)
    new_p = zeropoly(T, elim_vars...)

    # transfer coefficients and evalute var
    for (coeff, exponent) in p
        vec_exponent = collect(exponent)
        k = exponent[elim_var_idx]
        deleteat!(vec_exponent, elim_var_idx)
        new_exponent = tuple(vec_exponent...)
        new_p[new_exponent] += coeff * (value^k)
    end

    return new_p
end


"""
    embed(p::Poly{T,M}, vars::Vector{Symbol})

embeds a poly `p` into the ring T[vars].
"""
function embed(p::Poly{T,M}, newvars::Vector{Symbol}) where {T,M}
    pvars = vars(p)

    # create new poly of pvars
    poly = zeropoly(T, tuple(newvars...))
    polyvars = vars(poly)
    # store exponents
    exponent_storage = Dict(zip(polyvars, zeros(Int, length(newvars))))

    # convert the powers
    for (coeff, exponents) in p

        for (var, exponent) in zip(pvars, exponents)
            exponent_storage[var] = exponent
        end
        new_exponents = map(var -> exponent_storage[var], polyvars)

        poly[new_exponents] = coeff
    end
    poly
end

"""
    homogenize(p::Poly, var::Symbol)

Homogenizes p with the given variable, e.g. homogenize(X^2+X, :Y) == X^2+XY
"""
function homogenize(p::Poly{T,N}, newvar::Symbol) where {T,N}
    pvars = vars(p)

    # create new poly of pvars
    poly = zeropoly(T, tuple(newvar, pvars...))
    polyvars = vars(poly)
    # store exponents
    exponent_storage = Dict(zip(polyvars, zeros(Int, N + 1)))

    deg_p = deg(p)

    # convert the powers
    for (coeff, exponents) in p

        for (var, exponent) in zip(pvars, exponents)
            exponent_storage[var] = exponent
        end
        exponent_storage[newvar] = deg_p - sum(exponents)

        new_exponents = map(var -> exponent_storage[var], polyvars)
        poly[new_exponents] = coeff
    end
    poly
end

"""
    dehomogenize(p::Poly, var::Symbol)

Homogenizes p for the given variable, e.g. dehomogenize(X^2+XY, :Y) == X^2+X
"""
dehomogenize(p::Poly{T,N}, var::Symbol) where {T,N} = eval_var(p, var, one(T))

function diff_exp(a::NTuple{N,Int}, i::Int) where {N}
    k = a[i]

    if k > 0
        new_exp = ntuple(idx -> idx == i ? max(0, k - 1) : a[idx], N)

        return k, new_exp
    else
        return 0, ntuple(_ -> 0, N)
    end
end

"""
    gradient(p::Poly)

Returns the gradient of p.
"""
function gradient(p::Poly{T,N})::Vector{Poly{T,N}} where {T,N}
    polys = zeros(p, N)
    for (i, poly) in enumerate(polys)
        for (coeff, pow) in p
            k, new_pow = diff_exp(pow, i)
            new_coeff = coeff * k
            if new_coeff != zero(T)
                poly[new_pow] = new_coeff
            end
        end
    end
    polys
end

"""
    jacobi(p::Vector{Poly})

Returns the jacobi matrix of p.
"""
jacobi(F::Vector{Poly{T,N}}) where {T,N} = permutedims(hcat(map(diff, ps)...), [2, 1])
