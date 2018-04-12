module Utilities

import MultivariatePolynomials
const MP = MultivariatePolynomials

export allvariables,
    nvariables,
    ishomogenous,
    uniquevar,
    homogenize,
    totaldegree,
    solve_with_lu_inplace!,
    ProjectiveVector,
    sine_distance,
    dₛᵢₙ,
    riemann_distance,
    distance,
    batches


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
    solve_with_lu_inplace!(A, b)

Solves ``Ax =b`` inplace. The result is stored in `b`. This method also overrides
the contents of `A`.
"""
solve_with_lu_inplace!(A, b) = A_ldiv_B!(lufact!(A), b)


"""
    ProjectiveVector(v::Vector)

Wrap a vector `v` to indicate that it is an projective vector.
The vector `v` is *not* copied. It will also be normalized.
"""
struct ProjectiveVector{T} <: AbstractVector{T}
    data::Vector{T}

    ProjectiveVector{T}(data::Vector{T}) where T = new(normalize!(data))
end
ProjectiveVector(data::Vector{T}) where T = ProjectiveVector{T}(data)

Base.size(v::ProjectiveVector) = size(v.data)
Base.getindex(v::ProjectiveVector, i::Int) = getindex(v.data, i)
# Base.setindex!(v::ProjectiveVector, x, i::Int) = setindex!(v.data, x, Int)
Base.IndexStyle(::Type{<:ProjectiveVector}) = Base.IndexLinear()
Base.length(v::ProjectiveVector) = length(v.data)

"""
    riemann_distance(v::ProjectiveVector, w::ProjectiveVector)

Compute the Riemann distance between `v` and `w`. This is defined as
``\\arccos|⟨v, w⟩|``.

For details see: P. 283, Prop. 14.12 in
Bürgisser, Peter, and Felipe Cucker. Condition: The geometry of numerical algorithms. Vol. 349. Springer Science & Business Media, 2013.
"""
riemann_distance(v::ProjectiveVector, w::ProjectiveVector) = acos(clamp(abs(v ⋅ w), 0.0, 1.0))

"""
    sine_distance(v::ProjectiveVector, w::ProjectiveVector)

Compute the sine distance between `v` and `w`. This is defined as
``\\sin \\arccos|⟨v, w⟩|``. This defines a metric.

For details see: P. 284 in
Bürgisser, Peter, and Felipe Cucker. Condition: The geometry of numerical algorithms. Vol. 349. Springer Science & Business Media, 2013.
"""
sine_distance(v::ProjectiveVector, w::ProjectiveVector) = sin(riemann_distance(v, w))

"""
    dₛᵢₙ(v, w)

See  [`sine_distance`](v, w).
"""
dₛᵢₙ(v::ProjectiveVector, w::ProjectiveVector) = sine_distance(v, w)


"""
    distance(v, w)

Compute the distance between `v` and `w`.
"""
distance(v::ProjectiveVector, w::ProjectiveVector) = dₛᵢₙ(v, w)


# Parallelization

# This is into 0.7 but we need it for 0.6 as well
const get_num_BLAS_threads = function() # anonymous so it will be serialized when called
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


"""
    batches(nthreads, n)

Create a list of unit ranges with exponential dropoff to distribute `n` tasks
over `nthreads` Threads.
"""
function batches(nthreads, n)
    packages = Vector{UnitRange{Int}}()
    i = 1
    k = 1
    while i < n
        package_length = max(ceil(Int, n * 2.0^(-k) / nthreads), 1)
        for tid=1:nthreads
            r = i:min(n, i + package_length)
            push!(packages, r)
            i += length(r)
            if i ≥ n
                break
            end
        end
        k += 1
    end
    packages
end

end
