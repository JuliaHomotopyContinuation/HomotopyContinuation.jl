# Addition
function embed_polys{T,N,M}(p::Poly{T,N}, q::Poly{T,M})
    merged_vars = unique([vars(p)..., vars(q)...])
    embed(p, merged_vars), embed(q, merged_vars)
end

# assume p and q have the same variables
function add!{T,N}(p::Poly{T,N}, q::Poly{T,N})
    for (c, pow) in q
        p[pow] += c
    end
end

@inline function add{T,N}(p::Poly{T,N}, q::Poly{T,N})::Poly{T,N}
    result = deepcopy(p)
    add!(result, q)
    result
end


@inline function +{T,N,M}(p::Poly{T,N}, q::Poly{T,M})
    p_, q_ = embed_polys(p, q)
    add(p_,q_)
end

@inline function +{T,N}(p::Poly{T,N}, q::Poly{T,N})
    if p.vars == q.vars
        add(p, q)
    else
        p_, q_ = embed_polys(p, q)
        add(p_,q_)
    end
end


function +{T,N}(polys::Poly{T,N}...)
    result = deepcopy(first(polys))
    for p in polys[2:end]
        if p.vars === result.vars
            add!(result, p)
        else
            result, promoted_p = embed_polys(result, p)
            add!(result, promoted_p)
        end
    end
    result
end

function +{T,N}(p::Poly{T,N}, a::Number)::Poly{T,N}
    result = deepcopy(p)
    result[ntuple(_ -> 0, N)] += convert(T, a)
    result
end

+{T,N}(a::Number, p::Poly{T,N})::Poly{T,N} = p + a



# SUBSTRACTION
# assume p and q have the same variables
function minus!{T,N}(p::Poly{T,N}, q::Poly{T,N})
    for (c, pow) in q
        p[pow] -= c
    end
end

@inline function minus{T,N}(p::Poly{T,N}, q::Poly{T,N})::Poly{T,N}
    result = deepcopy(p)
    minus!(result, q)
    result
end


@inline function -{T,N,M}(p::Poly{T,N}, q::Poly{T,M})
    p_, q_ = embed_polys(p, q)
    minus(p_,q_)
end

@inline function -{T,N}(p::Poly{T,N}, q::Poly{T,N})
    if p.vars == q.vars
        minus(p, q)
    else
        p_, q_ = embed_polys(p, q)
        minus(p_,q_)
    end
end


function -{T,N}(polys::Poly{T,N}...)
    result = deepcopy(first(polys))
    for p in polys[2:end]
        if p.vars == result.vars
            minus!(result, p)
        else
            result, promoted_p = embed_polys(result, p)
            minus!(result, promoted_p)
        end
    end
    result
end

function -(p::Poly)
    result = zero(p)
    for (c, pow) in p
        result[pow] = -c
    end
    result
end

function -{T,N}(p::Poly{T,N}, a::Number)::Poly{T,N}
    result = deepcopy(p)
    result[ntuple(_ -> 0, N)] -= convert(T, a)
    result
end

-{T,N}(a::Number, p::Poly{T,N})::Poly{T,N} = -p + a


# MULTIPLICATION

add{T<:Number, N}(a::NTuple{N, T}, b::NTuple{N, T}) = ntuple(i -> a[i] + b[i], N)

# assume p and q have the same variables
function mul{T,N}(p::Poly{T,N}, q::Poly{T,N})::Poly{T,N}
    result = zero(p)
    for (a, a_exp) in p
        for (b, b_exp) in q
            result[add(a_exp, b_exp)] += a * b
        end
    end
    result
end

# TODO: Implement Karatsuba
function *{T,N,M}(p::Poly{T,N}, q::Poly{T,M})
    if p.vars == q.vars
        mul(p, q)
    else
        a, b = embed_polys(p, q)
        mul(a, b)
    end
end

function *{T<:Number, N}(a::Number, p::Poly{T,N})
    result = zero(p)
    for (coeff, pow) in p
        result[pow] = convert(T,a) * coeff
    end
    result
end
*{T<:Number, N}(p::Poly{T,N}, a::Number) = a * p


#powers
function ^(p::Poly, pow::Int)
    if (pow == 0)
        return one(p)
    elseif (pow == 1)
        return p
    elseif (pow > 1)
        a, r = divrem(pow, 2)
        return p^(a+r) * p^a
    else
        throw(ArgumentError("Can't take a negative power of polynomial!"))
    end
end
