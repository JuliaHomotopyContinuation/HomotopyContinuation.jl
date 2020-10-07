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
struct InterpretedSystem{T₁,T₂} <: AbstractSystem
    system::System
    eval_interpreter::Interpreter{T₁}
    eval_interpreter_cache::InterpreterCache{Vector{ComplexF64}}
    eval_interpreter_cache_ext::InterpreterCache{Vector{ComplexDF64}}
    jac_interpreter::Interpreter{T₂}
    jac_interpreter_cache::InterpreterCache{Vector{ComplexF64}}
    taylor_caches::Tuple{
        InterpreterCache{TaylorVector{2,ComplexF64}},
        InterpreterCache{TaylorVector{3,ComplexF64}},
        InterpreterCache{TaylorVector{4,ComplexF64}},
    }
end

function InterpretedSystem(F::System; optimizations::Bool = true)
    F = optimizations ? optimize(F) : F
    eval_interpreter = interpreter(F)
    jac_interpreter = jacobian_interpreter(F)
    N = cache_min_length(eval_interpreter)
    taylor_caches = (
        InterpreterCache(TaylorVector{2}(ComplexF64, N)),
        InterpreterCache(TaylorVector{3}(ComplexF64, N)),
        InterpreterCache(TaylorVector{4}(ComplexF64, N)),
    )

    InterpretedSystem(
        F,
        eval_interpreter,
        InterpreterCache(ComplexF64, eval_interpreter),
        InterpreterCache(ComplexDF64, eval_interpreter),
        jac_interpreter,
        InterpreterCache(ComplexF64, jac_interpreter),
        taylor_caches,
    )
end

Base.size(F::InterpretedSystem) = size(F.system)
variables(F::InterpretedSystem) = variables(F.system)
parameters(F::InterpretedSystem) = parameters(F.system)
variable_groups(F::InterpretedSystem) = variable_groups(F.system)
Base.:(==)(F::InterpretedSystem, G::InterpretedSystem) = F.system == G.system

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
function jacobian!(U, F::InterpretedSystem, x, p = nothing, cache = F.jac_interpreter_cache)
    execute!(nothing, U, F.jac_interpreter, x, p, cache)
    nothing
end

for M = 1:3
    @eval function taylor!(
        u::AbstractVector,
        Order::Val{$M},
        F::InterpretedSystem,
        x,
        p = nothing,
    )
        execute!(u, Order, F.eval_interpreter, x, p, F.taylor_caches[$M])
        u
    end
end

"""
    InterpretedHomotopy <: AbstractHomotopy

An [`AbstractHomotopy`](@ref) which automatically generates a program for the
fast evaluation of `H` and its Jacobian. The program is however, not compiled
but rather interpreted. See also [`CompiledHomotopy`](@ref).

    InterpretedHomotopy(H::Homotopy; optimizations = true)

Construct an `InterpretedHomotopy` from the given [`Homotopy`](@ref) `H`.
If `optimizations = true` then [`optimize`](@ref) is called on `H` before compiling.
"""
struct InterpretedHomotopy{T₁,T₂} <: AbstractHomotopy
    homotopy::Homotopy
    eval_interpreter::Interpreter{T₁}
    eval_interpreter_cache::InterpreterCache{Vector{ComplexF64}}
    eval_interpreter_cache_ext::InterpreterCache{Vector{ComplexDF64}}
    jac_interpreter::Interpreter{T₂}
    jac_interpreter_cache::InterpreterCache{Vector{ComplexF64}}
    taylor_caches::Tuple{
        InterpreterCache{TaylorVector{2,ComplexF64}},
        InterpreterCache{TaylorVector{3,ComplexF64}},
        InterpreterCache{TaylorVector{4,ComplexF64}},
    }
end

function InterpretedHomotopy(H::Homotopy; optimizations::Bool = true)
    H = optimizations ? optimize(H) : H
    eval_interpreter = interpreter(H)
    jac_interpreter = jacobian_interpreter(H)
    N = cache_min_length(eval_interpreter)
    taylor_caches = (
        InterpreterCache(TaylorVector{2}(ComplexF64, N)),
        InterpreterCache(TaylorVector{3}(ComplexF64, N)),
        InterpreterCache(TaylorVector{4}(ComplexF64, N)),
    )

    InterpretedHomotopy(
        H,
        eval_interpreter,
        InterpreterCache(ComplexF64, eval_interpreter),
        InterpreterCache(ComplexDF64, eval_interpreter),
        jac_interpreter,
        InterpreterCache(ComplexF64, jac_interpreter),
        taylor_caches,
    )
end


Base.size(H::InterpretedHomotopy) = size(H.homotopy)
variables(H::InterpretedHomotopy) = variables(H.homotopy)
parameters(H::InterpretedHomotopy) = parameters(H.homotopy)
Homotopy(H::InterpretedHomotopy) = H.homotopy

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
    execute!(nothing, U, H.jac_interpreter, x, t, p, cache)
    nothing
end

for M = 1:3
    @eval function taylor!(
        u::AbstractVecOrMat,
        Order::Val{$M},
        H::InterpretedHomotopy,
        x,
        t,
        p::Union{AbstractVector,Nothing} = nothing,
    )
        execute!(u, Order, H.eval_interpreter, x, t, p, H.taylor_caches[$M])
        u
    end
end
