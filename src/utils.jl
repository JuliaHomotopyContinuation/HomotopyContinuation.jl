export FiniteException, write_solutions, read_solutions, write_parameters, read_parameters

"""
    write_solutions(filename, solutions)

Stores the as a plain text file onto disk.
The storage format is as follows. The first line indicates the number of solutions stored
followed by a blank line.
Then the solutions are stored where each solution is separated by a blank line.
Note that the solutions are *always* considered as complex numbers.
See [`read_solutions`](@ref) for reading the solution file back.

Note that this is the same file format as used in [Bertini](https://bertini.nd.edu).

## Example

```julia
julia> write_solutions("solutions.txt", [[1, 1], [-1, 2]])

shell> cat solutions.txt
2

1 0
1 0

-1 0
2 0

julia> read_solutions("solutions.txt")
2-element Array{Array{Complex{Int64},1},1}:
 [1 + 0im, 1 + 0im]
 [-1 + 0im, 2 + 0im]
```
"""
function write_solutions(filename, S::AbstractVector{<:AbstractVector{<:Number}})
    open(filename, "w") do f
        write(f, "$(length(S))\n")
        for s in S
            write(f, "\n")
            DelimitedFiles.writedlm(f, [real.(s) imag.(s)], ' ')
        end
    end
end


"""
    read_solutions(filename)

Read the solutions stored with [`write_solutions`](@ref).

## Example

```julia
julia> write_solutions("solutions.txt", [[1, 1], [-1, 2]])

julia> read_solutions("solutions.txt")
2-element Array{Array{Complex{Int64},1},1}:
 [1 + 0im, 1 + 0im]
 [-1 + 0im, 2 + 0im]
```
"""
function read_solutions(filename)
    A = DelimitedFiles.readdlm(filename; skipblanks = false)
    n = convert(Int, A[1, 1])
    # figure out number of variables
    nvars = 0
    for k = 4:size(A, 1)
        if isempty(A[k, 1])
            nvars = k - 3
            break
        end
    end
    map(3:nvars+1:size(A, 1)) do i
        map(i:i+nvars-1) do k
            A[k, 1] + im * A[k, 2]

        end
    end
end

"""
    write_parameters(filename, parameters)

Stores the parameters as a plain text file onto disk.
The storage format is as follows. The first line indicates the number of parameter values
stored followed by a blank line. Then the parameter values are stored.

See [`read_parameters`](@ref)

Note that this is the same file format as used in [Bertini](https://bertini.nd.edu).

## Example

```julia
julia> write_parameters("parameters.txt", [2.0, -3.2 + 2im])

shell> cat parameters.txt
2

2.0 0.0
-3.2 2.0

julia> read_parameters("parameters.txt")
2-element Array{Complex{Float64},1}:
  2.0 + 0.0im
 -3.2 + 2.0im
```
"""
function write_parameters(filename, p::Vector{<:Number})
    open(filename, "w") do f
        write(f, "$(length(p))\n\n")
        DelimitedFiles.writedlm(f, [real.(p) imag.(p)], ' ')
    end
end

"""
    read_parameters(filename)

Read the parameters stored with [`write_parameters`](@ref).

## Example

```julia
julia> write_parameters("parameters.txt", [2.0, -3.2 + 2im])

julia> read_parameters("parameters.txt")
2-element Array{Complex{Float64},1}:
  2.0 + 0.0im
 -3.2 + 2.0im
```
"""
function read_parameters(filename)
    A = DelimitedFiles.readdlm(filename)
    n = A[1, 1]
    map(1:n) do i
        A[i+1, 1] + im * A[i+1, 2]
    end
end

###
### Exceptions
###
struct KeywordArgumentException <: Exception
    key::Symbol
    given::Any
    msg::String
end
function KeywordArgumentException(key, given)
    KeywordArgumentException(key, given, "")
end
function Base.showerror(io::IO, E::KeywordArgumentException)
    print(io, "Invalid argument $(E.given) for $(E.key). ", E.msg)
end

struct FiniteException <: Exception
    dim::Int
end
function Base.showerror(io::IO, E::FiniteException)
    print(
        io,
        "FiniteException: The solution set of the given system has at least dimension " *
        "$(E.dim) > 0. Consider intersecting with an (affine) subspace of codimension " *
        "$(E.dim) to reduce to (possibly) finitely many solutions.",
    )
end

@noinline function unsupported_kwargs(kwargs)
    if !(isempty(kwargs))
        msg = join(["$k = $v" for (k, v) in pairs(kwargs)], ", ")
        @warn "Ingored unsupported keyword arguments: $msg"
    end
end


# Perform type piracy here until https://github.com/JuliaArrays/StructArrays.jl/issues/131
# is fixed
const _C64SA = StructArrays.StructArray{
    Complex{Float64},
    2,
    NamedTuple{(:re, :im),Tuple{Array{Float64,2},Array{Float64,2}}},
    Int64,
}
Base.@propagate_inbounds function Base.setindex!(s::_C64SA, vals, I::Int)
    @boundscheck checkbounds(s, I)
    StructArrays.foreachfield(
        (col, val) -> (@inbounds col[I] = val),
        s,
        convert(ComplexF64, vals),
    )
    s
end

"""
    rand_unitary_matrix(n::Int, T=ComplexF64)

Samples a `n × n` unitary Matrix uniformly from the space of all unitary n × n matrices.

See https://arxiv.org/abs/math-ph/0609050 for a derivation.
"""
function rand_unitary_matrix(n::Int, T::Type = ComplexF64)
    Z = randn(T, n, n) ./ sqrt(2)
    Q, R = LA.qr(Z)
    Λ = LA.diagm(0 => [R[i, i] / abs(R[i, i]) for i = 1:n])
    Q * Λ
end

fast_abs(z::Complex) = sqrt(abs2(z))
fast_abs(x::Real) = abs(x)

"Like `min(a,b)`` but ignoring any `NaN` values."
nanmin(a, b) = isnan(a) ? b : (isnan(b) ? a : min(a, b))
nanmax(a, b) = isnan(a) ? b : (isnan(b) ? a : max(a, b))

always_false(::Any) = false

"""
    all2(f, a, b)

Checks if `f(a[i], b[i]) == true` for `i=1:length(a)`.
"""
function all2(f::F, a::AbstractVector, b::AbstractVector) where {F}
    length(a) == length(b) || return false
    for (ai, bi) in zip(a, b)
        f(ai, bi) || return false
    end
    true
end


"""
    unpack(a::Union{Nothing, T}, b::T)

Returns `a` if it is not `nothing`, otherwise `b`.
"""
unpack(a::Union{Nothing,T}, b::T) where {T} = a === nothing ? b : a

plural(singularstr, n) = n == 1 ? singularstr : singularstr * "s"

"""
     print_fieldnames(io::IO, obj)

 A better default printing for structs.
 """
function print_fieldnames(io::IO, obj)
    println(io, typeof(obj), ":")
    for name in fieldnames(typeof(obj))
        print_fieldname(io, obj, name)
    end
end
function print_fieldname(io::IO, obj, name)
    if getfield(obj, name) !== nothing
        val = getfield(obj, name)
        print(io, " • ", name, " → ")
        if val isa AbstractFloat
            println(io, round(val; sigdigits = 5))
        else
            println(io, val)
        end
    end
end

"""
    @tspawnat tid -> task
Mimics `Base.Threads.@spawn`, but assigns the task to thread `tid`.
# Example
```julia
julia> t = @tspawnat 4 Threads.threadid()
Task (runnable) @0x0000000010743c70
julia> fetch(t)
4
```
"""
macro tspawnat(thrdid, expr)
    thunk = esc(:(() -> ($expr)))
    var = esc(Base.sync_varname)
    tid = esc(thrdid)
    quote
        if $tid < 1 || $tid > Threads.nthreads()
            throw(
                AssertionError(
                    "@tspawnat thread assignment ($($tid)) must be between 1 and Threads.nthreads() (1:$(Threads.nthreads()))",
                ),
            )
        end
        local task = Task($thunk)
        task.sticky = false
        ccall(:jl_set_task_tid, Cvoid, (Any, Cint), task, $tid - 1)
        if $(Expr(:isdefined, var))
            push!($var, task)
        end
        schedule(task)
        task
    end
end


mutable struct SegmentStepper
    start::ComplexF64
    target::ComplexF64
    abs_Δ::Float64
    forward::Bool
    # current
    s::Float64
    # proposed
    s′::Float64
end
SegmentStepper(start::Number, target::Number) =
    SegmentStepper(ComplexF64(start), ComplexF64(target))
SegmentStepper(start::ComplexF64, target::ComplexF64) =
    init!(SegmentStepper(start, target, 0.0, true, 0.0, 0.0), start, target)

init!(S::SegmentStepper, start::Number, target::Number) =
    init!(S, ComplexF64(start), ComplexF64(target))
function init!(S::SegmentStepper, start::ComplexF64, target::ComplexF64)
    S.start = start
    S.target = target
    S.abs_Δ = abs(target - start)
    S.forward = abs(start) < abs(target)
    S.s = S.s′ = S.forward ? 0.0 : S.abs_Δ
    S
end

is_done(S::SegmentStepper) = S.forward ? S.s == S.abs_Δ : S.s == 0.0

step_success!(S::SegmentStepper) = (S.s = S.s′; S)
function propose_step!(S::SegmentStepper, Δs::Real)
    if S.forward
        S.s′ = min(S.s + Δs, S.abs_Δ)
    else
        S.s′ = max(S.s - Δs, 0.0)
    end
    S
end
dist_to_target(S::SegmentStepper) = S.forward ? S.abs_Δ - S.s : S.s

function Base.getproperty(S::SegmentStepper, sym::Symbol)
    if sym == :Δs
        s = getfield(S, :s)
        s′ = getfield(S, :s′)
        Δs = getfield(S, :forward) ? s′ - s : s - s′
        return Δs
    elseif sym === :t
        s = getfield(S, :s)
        start = getfield(S, :start)
        target = getfield(S, :target)
        abs_Δ = getfield(S, :abs_Δ)
        forward = getfield(S, :forward)
        return _t_helper(start, target, s, abs_Δ, forward)
    elseif sym == :t′
        s′ = getfield(S, :s′)
        start = getfield(S, :start)
        target = getfield(S, :target)
        abs_Δ = getfield(S, :abs_Δ)
        forward = getfield(S, :forward)
        return _t_helper(start, target, s′, abs_Δ, forward)
    elseif sym == :Δt
        s = getfield(S, :s)
        s′ = getfield(S, :s′)
        start = getfield(S, :start)
        target = getfield(S, :target)
        abs_Δ = getfield(S, :abs_Δ)
        forward = getfield(S, :forward)
        if forward
            Δt = ((s′ - s) / abs_Δ) * (target - start)
        else
            Δt = ((s - s′) / abs_Δ) * (target - start)
        end
        return Δt
    else # fallback to getfield
        return getfield(S, sym)
    end
end
function _t_helper(start, target, s, Δ, forward)
    if forward
        if s == 0.0
            return start
        elseif s == Δ
            return target
        else
            return start + (s / Δ) * (target - start)
        end
    else
        if s == Δ
            return start
        elseif s == 0.0
            return target
        else
            return target + (s / Δ) * (start - target)
        end
    end
end

function Base.show(io::IO, val::SegmentStepper)
    print(io, "SegmentStepper:")
    for field in [:start, :target, :t, :Δt]
        print(io, "\n • ", field, " → ", getproperty(val, field))
    end
end

"""
    nthroot(x, n)

Compute the `n`-th root of `x`.
"""
function nthroot(x::Real, N::Integer)
    if N == 4
        √(√(x))
    elseif N == 2
        √(x)
    elseif N == 3
        cbrt(x)
    elseif N == 1
        x
    elseif N == 0
        one(x)
    else
        x^(1 / N)
    end
end


import Base: MPFR
import Base.MPFR: MPFRRoundingMode


function set!(x::BigFloat, y::BigFloat, r::MPFRRoundingMode = MPFR.ROUNDING_MODE[])
    ccall(
        (:mpfr_set, :libmpfr),
        Int32,
        (Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode),
        x,
        y,
        r,
    )
    x
end


function set!(x::BigFloat, d::Float64, r::MPFRRoundingMode = MPFR.ROUNDING_MODE[])
    ccall(
        (:mpfr_set_d, :libmpfr),
        Int32,
        (Ref{BigFloat}, Float64, MPFRRoundingMode),
        x,
        d,
        r,
    )
    x
end


# Basic arithmetic without promotion
for (fJ, fC) in ((:add!, :add), (:mul!, :mul))
    @eval begin
        # BigFloat
        function ($fJ)(z::BigFloat, x::BigFloat, y::BigFloat)
            ccall(
                ($(string(:mpfr_, fC)), :libmpfr),
                Int32,
                (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode),
                z,
                x,
                y,
                MPFR.ROUNDING_MODE[],
            )
            return z
        end

        function ($fJ)(z::BigFloat, x::BigFloat, c::Float64)
            ccall(
                ($(string(:mpfr_, fC, :_d)), :libmpfr),
                Int32,
                (Ref{BigFloat}, Ref{BigFloat}, Cdouble, MPFRRoundingMode),
                z,
                x,
                c,
                MPFR.ROUNDING_MODE[],
            )
            return z
        end
        ($fJ)(c::Float64, x::BigFloat) = ($fJ)(x, c)

        function ($fJ)(z::BigFloat, x::BigFloat, c::Int64)
            ccall(
                ($(string(:mpfr_, fC, :_si)), :libmpfr),
                Int32,
                (Ref{BigFloat}, Ref{BigFloat}, Clong, MPFRRoundingMode),
                z,
                x,
                c,
                MPFR.ROUNDING_MODE[],
            )
            return z
        end
        ($fJ)(c::Int64, x::BigFloat) = ($fJ)(x, c)

        # BigInt
        function ($fJ)(z::BigFloat, x::BigFloat, c::BigInt)
            ccall(
                ($(string(:mpfr_, fC, :_z)), :libmpfr),
                Int32,
                (Ref{BigFloat}, Ref{BigFloat}, Ref{BigInt}, MPFRRoundingMode),
                z,
                x,
                c,
                MPFR.ROUNDING_MODE[],
            )
            return z
        end
        ($fJ)(c::BigInt, x::BigFloat) = ($fJ)(x, c)
    end
end

for (fJ, fC) in ((:sub!, :sub), (:div!, :div))
    @eval begin
        # BigFloat
        function ($fJ)(z::BigFloat, x::BigFloat, y::BigFloat)
            ccall(
                ($(string(:mpfr_, fC)), :libmpfr),
                Int32,
                (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode),
                z,
                x,
                y,
                MPFR.ROUNDING_MODE[],
            )
            return z
        end

        # Float32/Float64
        function ($fJ)(z::BigFloat, x::BigFloat, c::Float64)
            ccall(
                ($(string(:mpfr_, fC, :_d)), :libmpfr),
                Int32,
                (Ref{BigFloat}, Ref{BigFloat}, Cdouble, MPFRRoundingMode),
                z,
                x,
                c,
                MPFR.ROUNDING_MODE[],
            )
            return z
        end
        function ($fJ)(z::BigFloat, c::Float64, x::BigFloat)
            ccall(
                ($(string(:mpfr_, :d_, fC)), :libmpfr),
                Int32,
                (Ref{BigFloat}, Cdouble, Ref{BigFloat}, MPFRRoundingMode),
                z,
                c,
                x,
                MPFR.ROUNDING_MODE[],
            )
            return z
        end

        # BigInt
        function ($fJ)(z::BigFloat, x::BigFloat, c::BigInt)
            ccall(
                ($(string(:mpfr_, fC, :_z)), :libmpfr),
                Int32,
                (Ref{BigFloat}, Ref{BigFloat}, Ref{BigInt}, MPFRRoundingMode),
                z,
                x,
                c,
                MPFR.ROUNDING_MODE[],
            )
            return z
        end
        # no :mpfr_z_div function
    end
end

function sub!(z::BigFloat, c::BigInt, x::BigFloat)
    ccall(
        (:mpfr_z_sub, :libmpfr),
        Int32,
        (Ref{BigFloat}, Ref{BigInt}, Ref{BigFloat}, MPFRRoundingMode),
        z,
        c,
        x,
        MPFR.ROUNDING_MODE[],
    )
    return z
end

function rem!(z::BigFloat, x::BigFloat, y::BigFloat, ::RoundingMode{:Nearest})
    ccall(
        (:mpfr_remainder, :libmpfr),
        Int32,
        (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode),
        z,
        x,
        y,
        MPFR.MPFRRoundNearest,
    )
    return z
end
