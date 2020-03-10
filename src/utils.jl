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
