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
    randomish_gamma


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
    totaldegree_solutions(F::Vector{<:MP.AbstractPolynomialLike})

Returns an iterator of the solutions of the total degree startsystem of `F`.
See [`totaldegree`](@ref) for more details.
"""
function totaldegree_solutions(F::Vector{<:MP.AbstractPolynomialLike})
    TotalDegreeSolutionIterator(MP.maxdegree.(F), ishomogenous(F))
end


"""
    TotalDegreeSolutionIterator(degrees, homogenous::Bool)

Given the `Vector`s `degrees` `TotalDegreeSolutionIterator` enumerates all solutions
of the system
```math
\\begin{align*}
    z_1^{d_1} - 1 &= 0 \\\\
    z_1^{d_2} - 1 &= 0 \\\\
    &\\vdots \\\\
    z_n^{d_n} - 1 &= 0 \\\\
\\end{align*}
```
where ``d_i`` is `degrees[i]`. If `homogenous` is `true` then
the solutions of the system will be embedded in projective space by the map
`x → [1.0; x]`.
"""
struct TotalDegreeSolutionIterator{Iter}
    degrees::Vector{Int}
    homogenous::Bool
    iterator::Iter
end
function TotalDegreeSolutionIterator(degrees::Vector{Int}, homogenous::Bool)
    iterator = Base.Iterators.product(map(d -> 0:d-1, degrees)...)
    TotalDegreeSolutionIterator(degrees, homogenous, iterator)
end

Base.start(iter::TotalDegreeSolutionIterator) = start(iter.iterator)
function Base.next(iter::TotalDegreeSolutionIterator, state)
    indices, nextstate = next(iter.iterator, state)

    value = Complex{Float64}[]
    if iter.homogenous
        push!(value, complex(1.0, 0.0))
    end
    for (i, k) in zip(1:length(indices), indices)
        push!(value, cis(2π * k / iter.degrees[i]))
    end
    value, nextstate
end
Base.done(iter::TotalDegreeSolutionIterator, state) = done(iter.iterator, state)
Base.length(iter::TotalDegreeSolutionIterator) = length(iter.iterator)
Base.eltype(iter::TotalDegreeSolutionIterator) = Vector{Complex{Float64}}

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
