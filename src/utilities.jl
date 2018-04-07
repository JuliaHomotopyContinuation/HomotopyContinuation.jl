module Utilities

import MultivariatePolynomials
const MP = MultivariatePolynomials

export allvariables,
    nvariables,
    ishomogenous,
    uniquevar,
    homogenize,
    totaldegree,
    solve_with_lu_inplace!


"""
    allvariables(polys)

Returns a sorted list of all variables occuring in `polys`.
"""
function allvariables(polys::Vector{<:MP.AbstractPolynomialLike})
    sort!(union(Iterators.flatten(MP.variables.(polys))), rev=true)
end

"""
    nvariables(polys)

Returns the number of variables occuring in `polys`.
"""
nvariables(polys::Vector{<:MP.AbstractPolynomialLike}) = length(allvariables(polys))


"""
    ishomogenous(f::MP.AbstractPolynomialLike)

Checks whether `f` is homogenous.

    ishomogenous(polys::Vector{MP.AbstractPolynomialLike})

Checks whether each polynomial in `polys` is homogenous.
"""
ishomogenous(f::MP.AbstractPolynomialLike) = MP.mindegree(f) == MP.maxdegree(f)
ishomogenous(F::Vector{<:MP.AbstractPolynomialLike}) = all(ishomogenous, F)

"""
    uniquevar(f::MP.AbstractPolynomialLike, tag=:x0)
    uniquevar(F::Vector{<:MP.AbstractPolynomialLike}, tag=:x0)

Creates a unique variable.
"""
uniquevar(f::MP.AbstractPolynomialLike, tag=:x0) = MP.similarvariable(f, gensym(tag))
uniquevar(F::Vector{<:MP.AbstractPolynomialLike}, tag=:x0) = uniquevar(F[1], tag)

"""
    homogenize(f::MP.AbstractPolynomial, variable=uniquevar(f))

Homogenize the polynomial `f` by using the given variable `variable`.

    homogenize(F::Vector{<:MP.AbstractPolynomial}, variable=uniquevar(F))

Homogenize each polynomial in `F` by using the given variable `variable`.
"""
function homogenize(f::MP.AbstractPolynomialLike, var=uniquevar(f))
    d = MP.maxdegree(f)
    MP.polynomial(map(t -> var^(d - MP.degree(t)) * t, MP.terms(f)))
end
homogenize(F::Vector{<:MP.AbstractPolynomialLike}, var=uniquevar(F)) = homogenize.(F, var)


"""
    totaldegree(F::Vector{<:MP.AbstractPolynomialLike})

Construct a total degree start system for the system `F`.
This is the system
```math
\\begin{align*}
    z_1^{d_1} &- 1\\\\
    z_1^{d_2} &- 1\\\\
    &\\vdots \\\\
    z_n^{d_n} &- 1\\\\
\\end{align*}
```
where ``d_i`` is the degree of the ``i``-th polynomial of `F`.

## Example
```julia
julia> @polyvar x y;
julia> totaldegree([x^2 + y + 1, x^2*y^2 - 2])
[x^2 - 1, y^4 - 1]
```
"""
function totaldegree(F::Vector{<:MP.AbstractPolynomialLike})
    out = zeros(F)
    vars = allvariables(F)
    if ishomogenous(F)
        @assert length(out) + 1 ≥ length(vars)
        for k=2:length(vars)
            d = MP.maxdegree(F[k - 1])
            out[k - 1] = vars[k]^d - vars[1]^d
        end
    else
        @assert length(out) ≥ length(vars)
        for k=1:length(vars)
            out[k] = vars[k]^MP.maxdegree(F[k]) - 1
        end
    end

    out
end

"""
    solve_with_lu_inplace!(A, b)

Solves ``Ax =b`` inplace. The result is stored in `b`. This method also overrides
the contents of `A`.
"""
solve_with_lu_inplace!(A, b) = A_ldiv_B!(lufact!(A), b)

end






filter_kwargs(predicate, kwargs) = filter( x -> predicate(first(x)), kwargs)

# we define our own widen
_widen(T) = widen(T)
# Wait until DoubleFloat64 is released
# _widen(::Type{Float64}) = FastDouble

affine(xs::AbstractVector) = xs[2:end] ./ x[1]

"""
    embed_projective_if_necessary(x, H)

Embeds a vector into the projective space if necessary, i.e. if it's length is one less
than the number of variables of `H`. `H` is assumed to be homogenized. After the (eventual)
embedding the value is normalized.
"""
function embed_projective_if_necessary!(x::AbstractVector{T}, H::AbstractHomotopy{T}) where T
    N = Homotopies.nvariables(H)
    n = length(x)
    if N - 1 == n
        unshift!(x, one(T))
    elseif N != n
        error("A start value has length $n. Excepted length $N or $(N-1).")
    end
    x
end

"""
    projectivenorm2(a, b)

Calculate the squared infinity norm of |a-b|, but first brings `a` and `b` on the same patch.
Brings fist `a` and `b` on the same patch by finding diving `a` through it's maximum value
(w.r.t. to the absolute value) with index `i` and diving `b` through `b[i]`.
Then computes the norm of the differences.
"""
function projectivenorm2(a::AbstractVector{T}, b::AbstractVector{T}) where T
    maxind = 1
    maxval = abs2(first(a))
    for i = 2:length(a)
        val = abs2(a[i])
        if val > maxval
            maxind = i
            maxval = val
        end
    end
    out = real(zero(T))
    adiv = a[maxind]
    bdiv = b[maxind]

    for i=1:length(a)
        out = max(out, abs2(a[i] / adiv - b[i] / bdiv))
    end
    out
end


"""
    UnitRootsIterator(r, n)

Construct an infinite iterator which returns the `n`-th scaled unit roots, i.e.
the values ``r⋅exp(i2πk/n)`` for ``k=0,1,...``.
"""
struct UnitRootsIterator
    radius::Float64
    order::Float64
end
UnitRootsIterator(r::Real, order::Real) = UnitRootsIterator(float(r), float(order))

Base.start(::UnitRootsIterator) = 0
Base.next(loop::UnitRootsIterator, k::Int) = (loop.radius *  exp(im * 2π * k / loop.order), k + 1)
Base.done(::UnitRootsIterator, ::Int) = false
Base.iteratorsize(::UnitRootsIterator) = Base.IsInfinite()
Base.iteratoreltype(::UnitRootsIterator) = Base.HasEltype()
Base.eltype(::UnitRootsIterator) = Complex128
