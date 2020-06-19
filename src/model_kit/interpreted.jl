export InterpretedSystem, InterpretedHomotopy

"""
    InterpretedSystem <: AbstractSystem

An [`AbstractSystem`](@ref) which automatically generates a program for the
fast evaluation of `F` and its Jacobian. The program is however, not compiled
but rather interpreted. See also [`CompiledSystem`](@ref).

    InterpretedSystem(F::System; optimizations = true)

Construct an `InterpretedSystem` from the given [`System`](@ref) `F`.
If `optimizations = true` then [`optimize`](@ref) is called on `F` before compiling.
"""
struct InterpretedSystem{T} <: AbstractSystem
    system::System
    eval_interpreter::Interpreter{T,1}
    eval_interpreter_cache::InterpreterCache{ComplexF64}
    eval_interpreter_cache_ext::InterpreterCache{ComplexDF64}
    jac_interpreter::Interpreter{T,2}
    jac_interpreter_cache::InterpreterCache{ComplexF64}
    taylor_interpreters::Dict{
        NTuple{3,Int},
        Tuple{Interpreter{T,1},InterpreterCache{ComplexF64}},
    }
end

function InterpretedSystem(F::System; optimizations::Bool = true)
    F = optimizations ? optimize(F) : F
    eval_interpreter = evaluate_interpreter(F)
    jac_interpreter = jacobian_interpreter(F)
    T = promote_type(eltype(eval_interpreter), eltype(jac_interpreter))
    if !(T isa Rational)
        T = promote_type(T, Float64)
    end
    eval_interpreter = convert(Interpreter{T,1}, eval_interpreter)
    jac_interpreter = convert(Interpreter{T,2}, jac_interpreter)
    eval_interpreter_cache = InterpreterCache(ComplexF64, eval_interpreter)
    eval_interpreter_cache_ext = InterpreterCache(ComplexDF64, eval_interpreter)

    taylor_interpreters =
        Dict{NTuple{3,Int},Tuple{Interpreter{T,1},InterpreterCache{ComplexF64}}}()

    InterpretedSystem(
        F,
        eval_interpreter,
        eval_interpreter_cache,
        eval_interpreter_cache_ext,
        jac_interpreter,
        InterpreterCache(ComplexF64, jac_interpreter),
        taylor_interpreters,
    )
end

Base.size(F::InterpretedSystem) = size(F.system)
variables(F::InterpretedSystem) = variables(F.system)
parameters(F::InterpretedSystem) = parameters(F.system)
variable_groups(F::InterpretedSystem) = variable_groups(F.system)

function Base.show(io::IO, F::InterpretedSystem)
    print(io, "Interpreted: ")
    show(io, F.system)
end

(F::InterpretedSystem)(x, p = nothing) = F.system(x, p)
function evaluate!(u, F::InterpretedSystem, x, p = nothing)
    execute!(u, F.eval_interpreter, x, p, F.eval_interpreter_cache)
end
function evaluate!(u, F::InterpretedSystem, x::AbstractVector{ComplexDF64}, p = nothing)
    execute!(u, F.eval_interpreter, x, p, F.eval_interpreter_cache_ext)
end
function evaluate_and_jacobian!(u, U, F::InterpretedSystem, x, p = nothing)
    execute!(u, U, F.jac_interpreter, x, p, F.jac_interpreter_cache)
    nothing
end
function jacobian!(
    U,
    F::InterpretedSystem,
    x,
    p = nothing,
    cache = F.jac_interpreter_cache,
)
    execute!(U, F.jac_interpreter, x, p, cache)
    nothing
end

function taylor!(
    u::AbstractVector,
    v::Val{M},
    F::InterpretedSystem{T},
    x,
    p = nothing,
) where {T,M}
    order_x = _order(x)
    order_p = _order(p)
    if !haskey(F.taylor_interpreters, (M, order_x, order_p))
        I = taylor_interpreter(
            F.system;
            order_out = M,
            order_x = order_x,
            order_p = order_p,
        )
        C = InterpreterCache(ComplexF64, I)
        F.taylor_interpreters[(M, order_x, order_p)] = (I, C)
    else
        I, C = F.taylor_interpreters[(M, order_x, order_p)]
    end
    execute!(u, I, x, p, C)
    u
end

_order(::Nothing) = 0
_order(::AbstractArray) = 0
_order(::TaylorVector{N}) where {N} = N - 1

"""
    InterpretedHomotopy <: AbstractHomotopy

An [`AbstractHomotopy`](@ref) which automatically generates a program for the
fast evaluation of `H` and its Jacobian. The program is however, not compiled
but rather interpreted. See also [`CompiledHomotopy`](@ref).

    InterpretedHomotopy(H::Homotopy; optimizations = true)

Construct an `InterpretedHomotopy` from the given [`Homotopy`](@ref) `H`.
If `optimizations = true` then [`optimize`](@ref) is called on `H` before compiling.
"""
struct InterpretedHomotopy{T} <: AbstractHomotopy
    homotopy::Homotopy
    eval_interpreter::Interpreter{T,1}
    eval_interpreter_cache::InterpreterCache{ComplexF64}
    eval_interpreter_cache_ext::InterpreterCache{ComplexDF64}
    jac_interpreter::Interpreter{T,2}
    jac_interpreter_cache::InterpreterCache{ComplexF64}
    taylor_interpreters::Dict{
        NTuple{3,Int},
        Tuple{Interpreter{T,1},InterpreterCache{ComplexF64}},
    }
end

function InterpretedHomotopy(H::Homotopy; optimizations::Bool = true)
    H = optimizations ? optimize(H) : H
    eval_interpreter = evaluate_interpreter(H)
    jac_interpreter = jacobian_interpreter(H)
    T = promote_type(eltype(eval_interpreter), eltype(jac_interpreter))
    if !(T isa Rational)
        T = promote_type(T, Float64)
    end
    eval_interpreter = convert(Interpreter{T,1}, eval_interpreter)
    jac_interpreter = convert(Interpreter{T,2}, jac_interpreter)
    eval_interpreter_cache = InterpreterCache(ComplexF64, eval_interpreter)
    eval_interpreter_cache_ext = InterpreterCache(ComplexDF64, eval_interpreter)

    taylor_interpreters =
        Dict{NTuple{3,Int},Tuple{Interpreter{T,1},InterpreterCache{ComplexF64}}}()

    InterpretedHomotopy(
        H,
        eval_interpreter,
        eval_interpreter_cache,
        eval_interpreter_cache_ext,
        jac_interpreter,
        InterpreterCache(ComplexF64, jac_interpreter),
        taylor_interpreters,
    )
end


Base.size(H::InterpretedHomotopy) = size(H.homotopy)
variables(H::InterpretedHomotopy) = variables(H.homotopy)
parameters(H::InterpretedHomotopy) = parameters(H.homotopy)

function Base.show(io::IO, H::InterpretedHomotopy)
    print(io, "Interpreted: ")
    show(io, H.homotopy)
end

(H::InterpretedHomotopy)(x, t, p = nothing) = H.homotopy(x, t, p)
function evaluate!(
    u,
    H::InterpretedHomotopy,
    x,
    t,
    p = nothing,
    cache = H.eval_interpreter_cache,
)
    execute!(u, H.eval_interpreter, x, t, p, cache)
end
function evaluate!(
    u,
    H::InterpretedHomotopy,
    x::AbstractVector{ComplexDF64},
    t,
    p = nothing,
)
    execute!(u, H.eval_interpreter, x, t, p, H.eval_interpreter_cache_ext)
end
function evaluate_and_jacobian!(
    u,
    U,
    H::InterpretedHomotopy,
    x,
    t,
    p = nothing,
    cache = H.jac_interpreter_cache,
)
    execute!(u, U, H.jac_interpreter, x, t, p, cache)
    nothing
end

function jacobian!(
    U,
    H::InterpretedHomotopy,
    x,
    t,
    p = nothing,
    cache = H.jac_interpreter_cache,
)
    execute!(U, H.jac_interpreter, x, t, p, cache)
    nothing
end

function taylor!(
    u::AbstractVector,
    v::Val{M},
    H::InterpretedHomotopy{T},
    x,
    t,
    p::Union{AbstractArray,Nothing} = nothing,
) where {T,M}
    order_x = _order(x)
    order_p = _order(p)
    if !haskey(H.taylor_interpreters, (M, order_x, order_p))
        I = taylor_interpreter(
            H.homotopy;
            order_out = M,
            order_x = order_x,
            order_p = order_p,
        )
        C = InterpreterCache(ComplexF64, I)
        H.taylor_interpreters[(M, order_x, order_p)] = (I, C)
    else
        I, C = H.taylor_interpreters[(M, order_x, order_p)]
    end
    execute!(u, I, x, t, p, C)
    u
end
