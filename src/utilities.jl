module Utilities

import MultivariatePolynomials
const MP = MultivariatePolynomials

export allvariables,
    nvariables,
    ishomogenous,
    uniquevar,
    homogenize,
    totaldegree,
    ldiv_lu!, blas_ldiv_lu!,
    infinity_norm,
    unsafe_infinity_norm,
    logabs,
    fastlog,
    batches,
    randomish_gamma,
    filterkwargs


"""
    ldiv_lu!(A, b)

Solve ``Ax=b`` inplace using a pure Julia LU-factorization. This overwrites `A` and `b`
and stores the result in `b`. This is faster than the BLAS version [`blas_ldiv_lu!`](@ref) but seems to loose
around 2 digits of precision.
"""
function ldiv_lu!(A::StridedMatrix, b::StridedVector)
    A_ldiv_B!(LinAlg.generic_lufact!(A), b)
    b
end

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
    infinity_norm(z)

Compute the ∞-norm of `z`. If `z` is a complex vector this is more efficient
than `norm(z, Inf)`.

    infinity_norm(z₁, z₂)

Compute the ∞-norm of `z₁-z₂`.
"""
infinity_norm(z::AbstractVector{<:Complex}) = sqrt(maximum(abs2, z))
function infinity_norm(z₁::AbstractVector{<:Complex}, z₂::AbstractVector{<:Complex})
    m = abs2(z₁[1] - z₂[1])
    n₁, n₂ = length(z₁), length(z₂)
    if n₁ ≠ n₂
        return convert(typeof(m), Inf)
    end
    @inbounds for k=2:n₁
        m = max(m, abs2(z₁[k] - z₂[k]))
    end
    sqrt(m)
end
unsafe_infinity_norm(v, w) = infinity_norm(v, w)


"""
    logabs(z)

The log absolute map `log(abs(z))`.
"""
logabs(z::Complex) = 0.5 * fastlog(abs2(z))
logabs(x) = fastlog(abs(x))

fastlog(z) = Base.Math.JuliaLibm.log(z)

function randomish_gamma()
    # Usually values near 1, i, -i, -1 are not good randomization
    # Therefore we artificially constrain the choices
    theta = rand() * 0.30 + 0.075 + (rand(Bool) ? 0.0 : 0.5)
    cis(2π * theta)
end

"""
    filterkwargs(kwargs, allowed_kwargs)

Remove all keyword arguments out of `kwargs` where the keyword is not contained
in `allowed_kwargs`.
"""
function filterkwargs(kwargs, allowed_kwargs)
    [kwarg for kwarg in kwargs if any(kw -> kw == first(kwarg), allowed_kwargs)]
end

# Parallelization

get_num_BLAS_threads() = convert(Int, _get_num_BLAS_threads())
# This is into 0.7 but we need it for 0.6 as well
const _get_num_BLAS_threads = function() # anonymous so it will be serialized when called
    blas = Base.LinAlg.BLAS.vendor()
    # Wrap in a try to catch unsupported blas versions
    try
        if blas == :openblas
            return ccall((:openblas_get_num_threads, Base.libblas_name), Cint, ())
        elseif blas == :openblas64
            return ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
        elseif blas == :mkl
            return ccall((:MKL_Get_Max_Num_Threads, Base.libblas_name), Cint, ())
        end

        # OSX BLAS looks at an environment variable
        if Sys.isapple()
            return ENV["VECLIB_MAXIMUM_THREADS"]
        end
    end

    return nothing
end

end
