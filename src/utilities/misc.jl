export SymmetricGroup, infinity_norm, infinity_distance, fubini_study

"""
    SymmetricGroup(n)

Group action of the symmetric group S(n). This does not contain the identity element.
"""
struct SymmetricGroup{N,M}
    permutations::NTuple{M, SVector{N,Int}} # We should get rid of the necessary of a tuple
end
SymmetricGroup(N::Int) = SymmetricGroup(permutations(Val(N)))
function permutations(::Val{N}) where {N}
    s = StaticArrays.MVector{N}(1:N)
    perms = Vector{SVector{N, Int}}()
    while true
        i = N - 1
        while i>=1 && s[i] >= s[i+1]
            i -= 1
        end
        if i > 0
            j = N
            while j > i && s[i] >= s[j]
                j -= 1
            end
            s[i], s[j] = s[j], s[i]
            reverse!(s, i+1)
        else
            s[1] = N+1
        end

        s[1] > N && break
        push!(perms, SVector(s))
    end
    tuple(perms...)
end

Base.eltype(::Type{SymmetricGroup{N}}) where {N} = SVector{N, Int}
Base.length(p::SymmetricGroup{N}) where {N} = length(p.permutations)
Base.iterate(p::SymmetricGroup) = iterate(p.permutations)
Base.iterate(p::SymmetricGroup, s) = iterate(p.permutations, s)

"""
    deprecatekwarg(oldkw, newkw)

Logs a deprecation warning and assigns `oldkw` to `newkw` if `oldkw` is not `nothing`.
"""
macro deprecatekwarg(oldkw, newkw)
    quote
        if $(esc(oldkw)) !== nothing
            old = $(Expr(:quote, oldkw))
            new = $(Expr(:quote, newkw))
            @warn("`$(old)=$($(esc(oldkw)))` is deprecated, use `$(new)=$($(esc(oldkw)))` instead.")
            $(esc(newkw)) = $(esc(oldkw))
        end
    end
end

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
    log₁₀(x)

Same as `log10` but faster.
"""
log₁₀(x) = log(x) / 2.302585092994046


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
         if getfield(obj, name) !== nothing
             println(io, " • ", name, " → ", getfield(obj, name))
         end
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
isrealvector(z::NTuple{N, T}, tol=1e-6) where {N,T} = isrealvector(SVector{N}(z), tol)

"""
    randseed(range=1_000:1_000_000)

Return a random seed in the range `range`.
"""
randseed(range=1_000:1_000_000) = rand(range)

"""

    infinity_distance(u, v)

Compute the ∞-norm of `u-v`.
"""
infinity_distance(u, v) = infinity_norm(u, v)

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

Computes the Fubini-Study distance between `x` and `y`.
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


start_solution_sample(xs) = promote_start_solution(first(xs))
start_solution_sample(x::AbstractVector{<:Number}) = promote_start_solution(x)
promote_start_solution(x::AbstractVector{ComplexF64}) = x
promote_start_solution(x::SVector) = map(xᵢ -> first(promote(xᵢ, 0.0im)), x)
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
