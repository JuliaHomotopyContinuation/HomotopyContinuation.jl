const TSYSTEM_TABLE = Dict{UInt,Vector{System}}()

"""
    CompiledSystem <: AbstractSystem

An [`AbstractSystem`](@ref) which automatically compiles a straight line program for the
fast evaluation of a system `F` and its Jacobian.
For large systems the compilation can take some time and require a large amount of memory.
If this is a problem consider [`InterpretedSystem`](@ref).

    CompiledSystem(F::System; optimizations = true)

Construct a `CompiledSystem` from the given [`System`](@ref) `F`. If `optimizations = true`
then [`optimize`](@ref) is called on `F` before compiling.
"""
struct CompiledSystem{HI} <: AbstractSystem
    nexpressions::Int
    nvariables::Int
    nparameters::Int
    system::System
end

function CompiledSystem(F::System; optimizations::Bool = true)
    F = optimizations ? optimize(F) : F
    n = length(F)
    nvars = nvariables(F)
    nparams = nparameters(F)

    # We substitute all variable and parameter names away such that two systems
    # compile to the same code independently of the names
    @var _x_[1:nvars] _p_[1:nparams]
    D = ExpressionMap()
    for (x, y) in zip(variables(F), _x_)
        D[x] = y
    end
    for (x, y) in zip(parameters(F), _p_)
        D[x] = y
    end
    cleard_exprs = subs(expressions(F), D)
    sys = System(cleard_exprs, _x_, _p_)
    h = hash(cleard_exprs)
    k = 0
    if haskey(TSYSTEM_TABLE, h)
        # check that it is identical
        for (i, vi) in enumerate(TSYSTEM_TABLE[h])
            if vi == sys
                k = i
                break
            end
        end
        if k == 0
            push!(TSYSTEM_TABLE[h], sys)
            k = length(TSYSTEM_TABLE[h])
        end
    else
        k = 1
        TSYSTEM_TABLE[h] = [sys]
    end
    return CompiledSystem{(h, k)}(n, nvars, nparams, F)
end
function Base.show(io::IO, TS::CompiledSystem)
    print(io, "Compiled: ")
    show(io, TS.system)
end

(F::CompiledSystem)(x, p = nothing) = F.system(x, p)

interpret(TS::CompiledSystem) = TS.system
interpret(::Type{CompiledSystem{HI}}) where {HI} = TSYSTEM_TABLE[first(HI)][last(HI)]

Base.size(CS::CompiledSystem) = (CS.nexpressions, CS.nvariables)
Base.size(CS::CompiledSystem, i::Integer) = size(CS)[i]
Base.length(CS::CompiledSystem) = CS.nexpressions
nparameters(CS::CompiledSystem) = CS.nparameters
nvariables(CS::CompiledSystem) = CS.nvariables

variables(F::CompiledSystem) = variables(F.system)
parameters(F::CompiledSystem) = parameters(F.system)
variable_groups(F::CompiledSystem) = variable_groups(F.system)
System(F::CompiledSystem) = F.system

Base.:(==)(F::CompiledSystem{A}, G::CompiledSystem{B}) where {A,B} = A == B

######################
## CompiledHomotopy ##
######################

const THOMOTOPY_TABLE = Dict{UInt,Vector{Homotopy}}()

"""
    CompiledHomotopy <: AbstractHomotopy

An [`AbstractHomotopy`](@ref) which automatically compiles a straight line program for the
fast evaluation of a homotopy `H` and its Jacobian.
For large homotopies the compilation can take some time and require a large amount of
memory. If this is a problem consider [`InterpretedHomotopy`](@ref).

    CompiledSystem(F::System; optimizations = true)

Construct a `CompiledSystem` from the given [`System`](@ref) `F`. If `optimizations = true`
then [`optimize`](@ref) is called on `F` before compiling.
"""
struct CompiledHomotopy{HI} <: AbstractHomotopy
    nexpressions::Int
    nvariables::Int
    nparameters::Int
    homotopy::Homotopy
end

function CompiledHomotopy(H::Homotopy)
    n = length(H)
    nvars = nvariables(H)
    nparams = nparameters(H)

    # We substitute all variable and parameter names away such that two homotopies
    # compile to the same code independently of the names
    @var _x_[1:nvars] _t_ _p_[1:nparams]
    D = ExpressionMap()
    for (x, y) in zip(variables(H), _x_)
        D[x] = y
    end
    for (x, y) in zip(parameters(H), _p_)
        D[x] = y
    end
    D[H.t] = _t_
    cleard_exprs = subs(expressions(H), D)
    homotopy = Homotopy(cleard_exprs, _x_, _t_, _p_)
    h = hash(cleard_exprs)

    k = 0
    if haskey(THOMOTOPY_TABLE, h)
        # check that it is identical
        for (i, vi) in enumerate(THOMOTOPY_TABLE[h])
            if vi == homotopy
                k = 1
                break
            end
        end
        if k == 0
            push!(THOMOTOPY_TABLE[h], homotopy)
            k = 1
        end
    else
        k = 1
        THOMOTOPY_TABLE[h] = [homotopy]
    end
    return CompiledHomotopy{(h, 1)}(n, nvars, nparams, H)
end

Base.size(CH::CompiledHomotopy) = (CH.nexpressions, CH.nvariables)
Base.size(CH::CompiledHomotopy, i::Integer) = size(CH)[i]
Base.length(CH::CompiledHomotopy) = CH.nexpressions
nparameters(CH::CompiledHomotopy) = CH.nparameters
nvariables(CH::CompiledHomotopy) = CH.nvariables

function Base.show(io::IO, TH::CompiledHomotopy)
    print(io, "Compiled: ")
    show(io, TH.homotopy)
end

interpret(CH::CompiledHomotopy) = CH.homotopy
interpret(::Type{CompiledHomotopy{HI}}) where {HI} = THOMOTOPY_TABLE[first(HI)][last(HI)]


################
## EVALUATION ##
################
function _evaluate!_impl(T)
    I = interpret(T)
    list, out = instruction_list(expressions(I))
    generate_evaluate_impl(
        list,
        out;
        variables = Symbol.(variables(I)),
        parameters = Symbol.(parameters(I)),
        t = I isa Homotopy ? Symbol(I.t) : nothing,
    )
end

"""
    evaluate!(u, T::CompiledSystem, x, p = nothing)

Evaluate `T` for variables `x` and parameters `p` and store result in `u`.
"""
@generated function evaluate!(u, T::CompiledSystem, x, p = nothing)
    _evaluate!_impl(T)
end

"""
    evaluate!(u, T::CompiledHomotopy, x, t, p = nothing)

Evaluate `T` for variables `x`, `t` and parameters `p` and store result in `u`.
"""
@generated function evaluate!(u, T::CompiledHomotopy, x, t, p = nothing)
    _evaluate!_impl(T)
end

function _evaluate_and_jacobian!_impl(T)
    I = interpret(T)
    list, out = instruction_list(expressions(I))
    generate_evaluate_and_jacobian_impl(
        list,
        out;
        variables = Symbol.(variables(I)),
        parameters = Symbol.(parameters(I)),
        t = I isa Homotopy ? Symbol(I.t) : nothing,
    )
end

"""
    evaluate_and_jacobian!(u, U, T::CompiledSystem, x, p = nothing)

Evaluate `T` and its Jacobian for variables `x` and parameters `p` and
store result in `u`.
"""
@generated function evaluate_and_jacobian!(u, U, T::CompiledSystem, x, p = nothing)
    _evaluate_and_jacobian!_impl(T)
end
#
"""
    evaluate_and_jacobian!(u, U, T::CompiledHomotopy, x, t, p = nothing)

Evaluate `T` and its Jacobian for variables `x`, `t` and parameters `p` and
store result in `u`.
"""
@generated function evaluate_and_jacobian!(u, U, T::CompiledHomotopy, x, t, p = nothing)
    _evaluate_and_jacobian!_impl(T)
end

#
"""
    jacobian!(U, T::CompiledSystem, x, p = nothing)

Evaluate the Jacobian of `T` for variables `x`, `t` and parameters `p`
and store result in `u`.
"""
function jacobian!(U, T::CompiledSystem, x, p = nothing)
    evaluate_and_jacobian!(nothing, U, T, x, p)
end

"""
    jacobian!(U, T::CompiledHomotopy, x, t, p = nothing)

Evaluate the Jacobian of `T` for variables `x`, `t` and parameters `p` and
store result in `u`.
"""
function jacobian!(U, T::CompiledHomotopy, x, t, p = nothing)
    evaluate_and_jacobian!(nothing, U, T, x, t, p)
end


function _taylor!_impl(T)
    I = interpret(T)
    list, out = instruction_list(expressions(I))
    generate_taylor_impl(
        list,
        out;
        variables = Symbol.(variables(I)),
        parameters = Symbol.(parameters(I)),
        t = I isa Homotopy ? Symbol(I.t) : nothing,
    )
end

"""
    taylor!(
        u,
        ::Val{K},
        F::CompiledSystem,
        x,
        p = nothing,
    )

Compute the Taylor series order order `K` of ``u = F(x,p)``.
"""
@generated function taylor!(
    u::AbstractArray,
    Order::Val{K},
    T::CompiledSystem,
    x::AbstractArray,
    p::Union{Nothing,AbstractArray} = nothing,
) where {K}
    _taylor!_impl(T)
end

"""
    taylor!(
        u,
        ::Val{K},
        H::CompiledHomotopy,
        x,
        p = nothing,
    )

Compute the Taylor series order order `K` of ``u = H(x,t,p)``.
"""
@generated function taylor!(
    u::AbstractArray,
    Order::Val{K},
    T::CompiledHomotopy,
    x::AbstractArray,
    t,
    p::Union{Nothing,AbstractArray} = nothing,
) where {K}
    _taylor!_impl(T)
end
