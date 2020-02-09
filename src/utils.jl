fast_abs(z::Complex) = sqrt(abs2(z))
fast_abs(x::Real) = abs(x)

"""
     print_fieldnames(io::IO, obj)

 A better default printing for structs.
 """
function print_fieldnames(io::IO, obj)
    println(io, typeof(obj), ":")
    for name in fieldnames(typeof(obj))
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
end


"""
    ComplexLineSegment(a, b)

Models a straight line between two points `a` and `b` in the complex plane.
"""
struct ComplexLineSegment
    a::ComplexF64
    b::ComplexF64
    # derived
    Δ_b_a::ComplexF64
    abs_b_a::Float64
end

function ComplexLineSegment(start::Number, target::Number)
    a = ComplexF64(start)
    b = ComplexF64(target)
    Δ_b_a = b - a
    abs_b_a = abs(Δ_b_a)

    ComplexLineSegment(start, target, Δ_b_a, abs_b_a)
end

function Base.getindex(segment::ComplexLineSegment, t::Real)
    if t == segment.abs_b_a
        return segment.b
    elseif iszero(t)
        return segment.a
    else
        Δ = DoubleF64(t) / segment.abs_b_a
        ComplexF64(segment.a + Δ * segment.Δ_b_a)
    end
end
Base.length(segment::ComplexLineSegment) = segment.abs_b_a

function step_size(segment::ComplexLineSegment, Δs::Real)
    (Δs / segment.abs_b_a) * segment.Δ_b_a
end

function Base.show(io::IO, segment::ComplexLineSegment)
    print(io, "ComplexLineSegment($(segment.a), $(segment.b))")
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
