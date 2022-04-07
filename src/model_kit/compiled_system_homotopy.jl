const TSYSTEM_TABLE = Dict{UInt,Vector{System}}()

"""
    CompiledSystem <: AbstractSystem

An [`AbstractSystem`](@ref) which automatically compiles a straight line program for the
fast evaluation of a system `F` and its Jacobian.
For large systems the compilation can take some time and require a large amount of memory.
If this is a problem consider [`InterpretedSystem`](@ref).

    CompiledSystem(F::System)

Construct a `CompiledSystem` from the given [`System`](@ref) `F`.
"""
struct CompiledSystem{HI} <: AbstractSystem
    nexpressions::Int
    nvariables::Int
    nparameters::Int
    system::System
end

function CompiledSystem(F::System)
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

Base.:(==)(::CompiledSystem{A}, ::CompiledSystem{B}) where {A,B} = A == B

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

    CompiledSystem(F::System)

Construct a `CompiledSystem` from the given [`System`](@ref) `F`.
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
                k = i
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
    return CompiledHomotopy{(h, k)}(n, nvars, nparams, H)
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
function boundschecks(; nexpressions::Int, nvariables, nparameters, has_second_output)
    checks = Expr[]
    if has_second_output
        push!(checks, :(@boundscheck checkbounds(U, 1:$nexpressions, 1:$nvariables)))
    end
    push!(checks, :(@boundscheck u === nothing || checkbounds(u, 1:$nexpressions)))
    push!(checks, :(@boundscheck checkbounds(x, 1:$nvariables)))
    push!(checks, :(@boundscheck p === nothing || checkbounds(p, 1:$nparameters)))

    Expr(:block, checks...)
end


function compiled_execute_impl(
    seq::InstructionSequence;
    has_second_output = false,
    taylor = false,
)
    if taylor
        expr, get_var_name = sequence_to_expr(seq; op_call = taylor_op_call, order = :Order)
    else
        expr, get_var_name = sequence_to_expr(seq; op_call = op_call)
    end
    assignments = if has_second_output
        quote
            zero!(U)
            if isnothing(u)
                $(map(seq.assignments) do (i, k)
                    if i > seq.output_dim
                        :(U[$(i - seq.output_dim)] = $(get_var_name(k)))
                    end
                end...)
            else
                zero!(u)
                idx = CartesianIndices(($(seq.output_dim), size(U, 2)))
                $(map(seq.assignments) do (i, k)
                    if i <= seq.output_dim
                        :(u[$(i)] = $(get_var_name(k)))
                    else
                        :(U[idx[$(i - seq.output_dim)]] = $(get_var_name(k)))
                    end
                end...)
            end
        end
    elseif taylor
        quote
            zero!(u)
            if (assign_highest_order_only)
                $(map(seq.assignments) do (i, k)
                    :(u[$(i)] = $(get_var_name(k))[K])
                end...)
            else
                $(map(seq.assignments) do (i, k)
                    :(u[$(i)] = $(get_var_name(k)))
                end...)
            end


        end

    else
        quote
            zero!(u)
            $(map(seq.assignments) do (i, k)
                :(u[$(i)] = $(get_var_name(k)))
            end...)
        end
    end

    quote
        $(boundschecks(
            nexpressions = seq.output_dim,
            has_second_output = has_second_output,
            nparameters = length(seq.parameters_range),
            nvariables = length(seq.variables_range),
        ))
        @inbounds begin
            $(expr)
            $(assignments)
        end
        u
    end
end

_evaluate!_impl(T) = compiled_execute_impl(instruction_sequence(interpret(T)))


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

_evaluate_and_jacobian!_impl(T) = compiled_execute_impl(
    jacobian_instruction_sequence(interpret(T));
    has_second_output = true,
)

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


_taylor!_impl(T) = compiled_execute_impl(instruction_sequence(interpret(T)); taylor = true)

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
    p::Union{Nothing,AbstractArray} = nothing;
    assign_highest_order_only::Bool = u isa Vector,
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
    t_,
    p::Union{Nothing,AbstractArray} = nothing;
    assign_highest_order_only::Bool = u isa Vector,
) where {K}
    quote
        t = (t_, 1)
        $(_taylor!_impl(T))
    end
end
