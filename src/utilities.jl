module Utilities

import LinearAlgebra
import MultivariatePolynomials
const MP = MultivariatePolynomials

export @moduleenum,
    nthroot,
    print_fieldnames,
    infinity_norm,
    unsafe_infinity_norm,
    fubini_study,
    logabs,
    batches,
    randomish_gamma,
    filterkwargs,
    splitkwargs,
    set_num_BLAS_threads,
    get_num_BLAS_threads,
    randseed,
    check_kwargs_empty,
    start_solution_sample,
    ComplexSegment

include("utilities/unique_points.jl")
include("utilities/linear_algebra.jl")
include("utilities/polynomials.jl")


"""
    nthroot(x, n)

Compute the `n`-th root of `x`.
"""
function nthroot(x::Real, N::Integer)
    if N == 3
        cbrt(x)
    elseif N == 2
        √(x)
    elseif N == 4
        √(√(x))
    elseif N == 1
        x
    elseif N == 0
        one(x)
    else
        x^(1/N)
    end
end

"""
    @modulenum(name, block)

This is a modification of `@enum` and creates an intermediate enum.

## Example

The definition

```julia
@moduleenum Car begin
    audi
    volkswagen
    bmw
end
```

expands into

```julia
module Car
    @enum t begin
        audi
        volkswagen
        bmw
    end
end
import .Car
"""
macro moduleenum(name, content...)
    mname = esc(name)
    quote
        @eval module $name
            @enum t $(content...)
        end
        import .$name
    end
end

"""
     print_fieldnames(io::IO, obj)

 A better default printing for structs.
 """
 function print_fieldnames(io::IO, obj)
     println(io, typeof(obj), ":")
     for name in fieldnames(typeof(obj))
         println(io, " • ", name, " → ", getfield(obj, name))
     end
 end


"""
    check_kwargs_empty(kwargs, [allowed_kwargs])

Chack that the list of `kwargs` is empty. If not, print all unsupported keywords
with their arguments.
"""
function check_kwargs_empty(kwargs, allowed_kwargs=[])
    if !isempty(kwargs)
        msg = "Unexpected keyword argument(s): "
        first_el = true
        for kwarg in kwargs
            if !first_el
                msg *= ", "
            end
            msg *= "$(first(kwarg))=$(last(kwarg))"
            first_el = false
        end
        if !isempty(allowed_kwargs)
            msg *= "\nAllowed keywords are\n"
            msg *= join(allowed_kwargs, ", ")
        end
        throw(ErrorException(msg))
    end
end

"""
    isrealvector(v::AbstractVector, tol=1e-6)

Check whether the 2-norm of the imaginary part of `v` is at most `tol`.
."""
isrealvector(z::AbstractVector{<:Real}, tol=1e-6) = true
function isrealvector(z::AbstractVector{<:Complex}, tol=1e-6)
    total = zero(real(eltype(z)))
    for zᵢ in z
        total += abs2(imag(zᵢ))
    end
    sqrt(total) < tol
end

"""
    randseed(range=1_000:1_000_000)

Return a random seed in the range `range`.
"""
randseed(range=1_000:1_000_000) = rand(range)

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


start_solution_sample(xs) = first(xs) |> promote_start_solution
start_solution_sample(x::AbstractVector{<:Number}) = promote_start_solution(x )

promote_start_solution(x::AbstractVector{ComplexF64}) = x
function promote_start_solution(x)
    x_new =similar(x, promote_type(eltype(x), ComplexF64), length(x))
    copyto!(x_new, x)
    x_new
end

"""
    ComplexSegment(start, target)

Represents the line segment from `start` to `finish`.
Supports indexing of the values `t ∈ [0, length(target-start)]` in order to get the
corresponding point on the line segment.
"""
struct ComplexSegment
    start::ComplexF64
    target::ComplexF64
    # derived
    Δ_target_start::ComplexF64
    abs_target_start::Float64
end
function ComplexSegment(start, target)
    Δ_target_start = convert(ComplexF64, target) - convert(ComplexF64, start)
    abs_target_start = abs(Δ_target_start)

    ComplexSegment(start, target, Δ_target_start, abs_target_start)
end

function Base.getindex(segment::ComplexSegment, t::Real)
    Δ = t / segment.abs_target_start
    if 1.0 - Δ < 2eps()
        Δ = 1.0
    end
    segment.start + Δ * segment.Δ_target_start
end
Base.length(segment::ComplexSegment) = segment.abs_target_start

function Base.show(io::IO, segment::ComplexSegment)
    print(io, "ComplexSegment($(segment.start), $(segment.target))")
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", opts::ComplexSegment) = opts

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
