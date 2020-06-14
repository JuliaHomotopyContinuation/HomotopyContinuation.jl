export InterpretedSystem

using .ModelKit: Interpreter, InterpreterCache

"""
    InterpretedSystem(F:System)

Construct a system from the given [`System`](@ref) `F`.
"""
mutable struct InterpretedSystem{T} <: AbstractSystem
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

function InterpretedSystem(F::ModelKit.System; optimizations::Bool = true)
    F = optimizations ? ModelKit.optimize(F) : F
    eval_interpreter = ModelKit.evaluate_interpreter(F)
    eval_interpreter_cache = ModelKit.InterpreterCache(ComplexF64, eval_interpreter)
    eval_interpreter_cache_ext = ModelKit.InterpreterCache(ComplexDF64, eval_interpreter)
    jac_interpreter = ModelKit.jacobian_interpreter(F)
    taylor_interpreters =
        Dict{NTuple{3,Int},Tuple{typeof(eval_interpreter),InterpreterCache{ComplexF64}}}()

    InterpretedSystem(
        F,
        eval_interpreter,
        eval_interpreter_cache,
        eval_interpreter_cache_ext,
        jac_interpreter,
        ModelKit.InterpreterCache(ComplexF64, jac_interpreter),
        taylor_interpreters,
    )
end


Base.size(F::InterpretedSystem) = size(F.system)
ModelKit.variables(F::InterpretedSystem) = variables(F.system)
ModelKit.parameters(F::InterpretedSystem) = parameters(F.system)
ModelKit.variable_groups(F::InterpretedSystem) = variable_groups(F.system)

function Base.show(io::IO, F::InterpretedSystem)
    println(io, typeof(F), ":")
    println(io, F.system)
end

(F::InterpretedSystem)(x, p = nothing) = F.system(x, p)
function evaluate!(u, F::InterpretedSystem, x, p = nothing)
    ModelKit.execute!(u, F.eval_interpreter, x, p, F.eval_interpreter_cache)
end
function evaluate!(u, F::InterpretedSystem, x::AbstractVector{ComplexDF64}, p = nothing)
    ModelKit.execute!(u, F.eval_interpreter, x, p, F.eval_interpreter_cache_ext)
end
function evaluate_and_jacobian!(u, U, F::InterpretedSystem, x, p = nothing)
    ModelKit.execute!(u, U, F.jac_interpreter, x, p, F.jac_interpreter_cache)
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
        I = ModelKit.taylor_interpreter(
            F.system;
            order_out = M,
            order_x = order_x,
            order_p = order_p,
        )
        C = ModelKit.InterpreterCache(ComplexF64, I)
        F.taylor_interpreters[(M, order_x, order_p)] = (I, C)
    else
        I, C = F.taylor_interpreters[(M, order_x, order_p)]
    end
    ModelKit.execute!(u, I, x, p, C)
    u
end

_order(::Nothing) = 0
_order(::AbstractArray) = 0
_order(::TaylorVector{N}) where {N} = N - 1
