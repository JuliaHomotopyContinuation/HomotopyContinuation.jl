module Utilities

import LinearAlgebra
import MultivariatePolynomials
const MP = MultivariatePolynomials

export allvariables,
    nvariables,
    ishomogenous,
    uniquevar,
    homogenize,
    ldiv_lu!, blas_ldiv_lu!,
    infinity_norm,
    unsafe_infinity_norm,
    fubini_study,
    logabs,
    batches,
    randomish_gamma,
    filterkwargs,
    splitkwargs,
    solve!,
    set_num_BLAS_threads,
    get_num_BLAS_threads,
    check_zero_dimensional,
    randseed

"""
    check_zero_dimensional(F::Vector{<:MP.AbstractPolynomial}, homvar)

Check that the given polynomial system can have zero dimensional components.
"""
function check_zero_dimensional(F::Vector{<:MP.AbstractPolynomial}, homvar)
    N = homvar === nothing ? MP.nvariables(F) : MP.nvariables(F) - 1
    n = length(F)

    if n ≥ N ||
       n == N - 1 && ishomogenous(F)
        return
    end
    throw(AssertionError("The input system will not result in a finite number of solutions."))
end

"""
    randseed(range=1_000:1_000_000)

Return a random seed in the range `range`.
"""
randseed(range=1_000:1_000_000) = rand(range)

"""
    solve!(A, b)

Solve ``Ax=b`` inplace. This overwrites `A` and `b`
and stores the result in `b`.
"""
function solve!(A::StridedMatrix, b::StridedVecOrMat)
    m, n = size(A)
    if m == n
        LinearAlgebra.ldiv!(LinearAlgebra.generic_lufact!(A), b)
    else
        LinearAlgebra.ldiv!(LinearAlgebra.qr!(A), b)
    end
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
    ishomogenous(f::MP.AbstractPolynomialLike, v::Vector{<:MP.AbstractVariable})

Checks whether `f` is homogenous in the variables `v`.

    ishomogenous(polys::Vector{MP.AbstractPolynomialLike}, v::Vector{<:MP.AbstractVariable})

Checks whether each polynomial in `polys` is homogenous in the variables `v`.
"""
function ishomogenous(f::MP.AbstractPolynomialLike, variables::Vector{T}) where {T<:MP.AbstractVariable}
 var_indices = findall(in(variables), MP.variables(f))
 degrees = map(t -> sum(MultivariatePolynomials.exponents(t)[var_indices]), f)
 minimum(degrees) == maximum(degrees)
end
function ishomogenous(F::Vector{<:MP.AbstractPolynomialLike}, variables::Vector{T}) where {T<:MP.AbstractVariable}
 all(f -> ishomogenous(f, variables), F)
end

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
homogenize(F::Vector{<:MP.AbstractPolynomialLike}, var=uniquevar(F)) = homogenize.(F, Ref(var))

"""
    homogenize(f::MP.AbstractPolynomial, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(f))

Homogenize the variables `v` in the polynomial `f` by using the given variable `variable`.

    homogenize(F::Vector{<:MP.AbstractPolynomial}, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(F))

Homogenize the variables `v` in each polynomial in `F` by using the given variable `variable`.
"""
function homogenize(f::MP.AbstractPolynomialLike, variables::Vector{T}, var=uniquevar(f)) where {T<:MP.AbstractVariable}
    var_indices = findall(in(variables), MP.variables(f))
    degrees = map(t -> sum(MultivariatePolynomials.exponents(t)[var_indices]), f)
    d = maximum(degrees)
    MP.polynomial(map(i -> var^(d - degrees[i]) * f[i], 1:length(f)))
end
function homogenize(F::Vector{<:MP.AbstractPolynomialLike}, variables::Vector{T}, var=uniquevar(F)) where {T<:MP.AbstractVariable}
    map(f -> homogenize(f, variables, var), F)
end


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
    fubini_study(x, y)

Computes the Fubini-Study norm of `x` and `y`.
"""
fubini_study(x,y) = acos(min(1.0, abs(LinearAlgebra.dot(x,y))))

"""
    logabs(z)

The log absolute map `log(abs(z))`.
"""
logabs(z::Complex) = 0.5 * log(abs2(z))
logabs(x) = log(abs(x))


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

"""
    splitkwargs(kwargs, supported_keywords)

Split the vector of `kwargs` in two vectors, the first contains all `kwargs`
whose keywords appear in `supported_keywords` and the rest the other one.
"""
function splitkwargs(kwargs, supported_keywords)
    supported = []
    rest = []
    for kwarg in kwargs
        if any(kw -> kw == first(kwarg), supported_keywords)
            push!(supported, kwarg)
        else
            push!(rest, kwarg)
        end
    end
    supported, rest
end

# Parallelization

set_num_BLAS_threads(n) = LinearAlgebra.BLAS.set_num_threads(n)
get_num_BLAS_threads() = convert(Int, _get_num_BLAS_threads())
# This is into 0.7 but we need it for 0.6 as well
const _get_num_BLAS_threads = function() # anonymous so it will be serialized when called
    blas = LinearAlgebra.BLAS.vendor()
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
    catch
    end

    return nothing
end

end
