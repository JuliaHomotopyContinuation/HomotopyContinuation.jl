export compile

#####################
## CompiledSystem  ##
#####################
const TSYSTEM_TABLE = Dict{
    UInt,
    Vector{Tuple{Vector{Expression},Vector{Variable},Vector{Variable}}},
}()

struct CompiledSystem{HI}
    nexpressions::Int
    variables::Vector{Symbol}
    parameters::Vector{Symbol}
end

function CompiledSystem(
    exprs::Vector{Expression},
    var_order::Vector{Variable},
    param_order::Vector{Variable} = Variable[],
)
    val = (exprs, var_order, param_order)
    h = hash(val)

    n = length(exprs)
    vars = Symbol.(var_order)
    params = Symbol.(param_order)

    if haskey(TSYSTEM_TABLE, h)
        # check that it is identical
        for (i, vi) in enumerate(TSYSTEM_TABLE[h])
            if vi == val
                return CompiledSystem{(h, i)}(n, vars, params)
            end
        end
        push!(TSYSTEM_TABLE[h], val)
        return CompiledSystem{(h, length(TSYSTEM_TABLE[h]))}(n, vars, params)
    else
        TSYSTEM_TABLE[h] = [val]
        return CompiledSystem{(h, 1)}(n, vars, params)
    end
end

function Base.show(io::IO, TS::CompiledSystem)
    print(io, "Compiled: ")
    show(io, interpret(TS))
end

interpret(TS::CompiledSystem) = interpret(typeof(TS))
interpret(::Type{CompiledSystem{HI}}) where {HI} =
    System(TSYSTEM_TABLE[first(HI)][last(HI)]...)

Base.size(CS::CompiledSystem) = (CS.nexpressions, length(CS.variables))
Base.size(CS::CompiledSystem, i::Integer) = size(CS)[i]
Base.length(CS::CompiledSystem) = CS.nexpressions


######################
## CompiledHomotopy ##
######################

const THOMOTOPY_TABLE = Dict{
    UInt,
    Vector{Tuple{
        Vector{Expression},
        Vector{Variable},
        Variable,
        Vector{Variable},
    },},
}()

struct CompiledHomotopy{HI}
    nexpressions::Int
    variables::Vector{Symbol}
    parameters::Vector{Symbol}
end

function CompiledHomotopy(
    exprs::Vector{<:Expression},
    var_order::AbstractVector{<:Variable},
    t::Variable,
    param_order::AbstractVector{<:Variable} = Variable[],
)
    val = (exprs, var_order, t, param_order)
    h = hash(val)

    n = length(exprs)
    vars = Symbol.(var_order)
    params = Symbol.(param_order)

    if haskey(THOMOTOPY_TABLE, h)
    # check that it is identical
        for (i, vi) in enumerate(THOMOTOPY_TABLE[h])
            if vi == val
                return CompiledHomotopy{(h, i)}(n, vars, params)
            end
        end
        push!(THOMOTOPY_TABLE[h], val)
        return CompiledHomotopy{(h, length(TSYSTEM_TABLE[h]))}(n, vars, params)
    else
        THOMOTOPY_TABLE[h] = [val]
        return CompiledHomotopy{(h, 1)}(n, vars, params)
    end
end

Base.size(CH::CompiledHomotopy) = (CH.nexpressions, length(CH.variables))
Base.size(CS::CompiledHomotopy, i::Integer) = size(CS)[i]
Base.length(CH::CompiledHomotopy) = CH.nexpressions

function Base.show(io::IO, TH::CompiledHomotopy)
    print(io, "Compiled: ")
    show(io, interpret(TH))
end

interpret(TH::CompiledHomotopy) = interpret(typeof(TH))
interpret(::Type{CompiledHomotopy{HI}}) where {HI} =
    Homotopy(THOMOTOPY_TABLE[first(HI)][last(HI)]...)

"""
    compile(F::System)

Compile the given system. Returns a `CompiledSystem`.
    compile(H::Homotopy)

Compile the given homotopy. Returns a `CompiledHomotopy`.
"""
compile(F::System) = CompiledSystem(F.expressions, F.variables, F.parameters)
compile(H::Homotopy) =
    CompiledHomotopy(H.expressions, H.variables, H.t, H.parameters)

#############
## CODEGEN ##
#############

boundscheck_var_map(F::System; kwargs...) =
    boundscheck_var_map(F.expressions, F.variables, F.parameters; kwargs...)
boundscheck_var_map(H::Homotopy; kwargs...) = boundscheck_var_map(
    H.expressions,
    H.variables,
    H.parameters,
    H.t;
    kwargs...,
)

function boundscheck_var_map(
    exprs,
    vars,
    params,
    t = nothing;
    jacobian::Bool = false,
)
    n = length(exprs)
    m = length(vars)
    l = length(params)
    var_map = Dict{Symbol,Union{Symbol,Expr}}()
    for i = 1:m
        var_map[Symbol(vars[i])] = :(x[$i])
    end
    for i = 1:l
        var_map[Symbol(params[i])] = :(p[$i])
    end
    if t !== nothing
        var_map[Symbol(t)] = :(t)
    end

    checks = Expr[]
    if jacobian
        push!(checks, :(@boundscheck checkbounds(U, 1:$n, 1:$m)))
    else
        push!(checks, :(@boundscheck checkbounds(u, 1:$n)))
    end
    push!(checks, :(@boundscheck checkbounds(x, 1:$m)))
    push!(checks, :(@boundscheck p === nothing || checkbounds(p, 1:$l)))

    Expr(:block, checks...), var_map
end

function add_assignement!(D::Dict{Symbol,Vector{Expr}}, id::Symbol, e::Expr)
    if haskey(D, id)
        push!(D[id], e)
    else
        D[id] = [e]
    end
    D
end

function _evaluate!_impl(::Type{T},) where {T<:Union{
    CompiledSystem,
    CompiledHomotopy,
}}
    I = interpret(T)
    checks, var_map = boundscheck_var_map(I)
    slp = let
        list, ids = instruction_list(I.expressions)
        assignements = Dict{Symbol,Vector{Expr}}()
        for (i, id) in enumerate(ids)
            add_assignement!(assignements, id, :(u[$i] = $id))
        end
        to_expr(list, var_map, assignements)
    end

    quote
        $checks
        @inbounds $slp
        u
    end
end

function _jacobian!_impl(::Type{T},) where {T<:Union{
    CompiledSystem,
    CompiledHomotopy,
}}
    I = interpret(T)
    checks, var_map = boundscheck_var_map(I; jacobian = true)

    slp = let
        list, ids = instruction_list(I.expressions)
        vars = Symbol.(I.variables)
        params = Symbol.(I.parameters)
        dlist, J = diff(list, vars, ids)

        assignements = Dict{Symbol,Vector{Expr}}()

        U_constants = Expr[]
        for j = 1:size(J, 2), i = 1:size(J, 1)
            if J[i, j] isa Symbol
                if J[i, j] ∉ vars && J[i, j] ∉ params
                    add_assignement!(
                        assignements,
                        J[i, j],
                        :(U[$i, $j] = $(J[i, j])),
                    )
                else
                    push!(U_constants, :(U[$i, $j] = $(var_map[J[i, j]])))
                end
            elseif J[i, j] !== nothing
                push!(U_constants, :(U[$i, $j] = $(J[i, j])))
            end
        end
        expr = to_expr(dlist, var_map, assignements)
        append!(expr.args, U_constants)
        expr
    end
    quote
        $checks
        U .= zero(eltype(x))
        @inbounds $slp
        U
    end
end

function _evaluate_and_jacobian!_impl(::Type{T},) where {T<:Union{
    CompiledSystem,
    CompiledHomotopy,
}}
    I = interpret(T)
    checks, var_map = boundscheck_var_map(I; jacobian = true)

    slp = let
        list, ids = instruction_list(I.expressions)
        vars = Symbol.(I.variables)
        params = Symbol.(I.parameters)
        dlist, J = diff(list, vars, ids)

        assignements = Dict{Symbol,Vector{Expr}}()
        for (i, id) in enumerate(ids)
            add_assignement!(assignements, id, :(u[$i] = $id))
        end

        U_constants = Expr[]
        for j = 1:size(J, 2), i = 1:size(J, 1)
            if J[i, j] isa Symbol
                if J[i, j] ∉ vars && J[i, j] ∉ params
                    add_assignement!(
                        assignements,
                        J[i, j],
                        :(U[$i, $j] = $(J[i, j])),
                    )
                else
                    push!(U_constants, :(U[$i, $j] = $(var_map[J[i, j]])))
                end
            elseif J[i, j] !== nothing
                push!(U_constants, :(U[$i, $j] = $(J[i, j])))
            end
        end
        expr = to_expr(dlist, var_map, assignements)
        append!(expr.args, U_constants)
        expr
    end
    quote
        $checks
        U .= zero(eltype(x))
        @inbounds $slp
        nothing
    end
end

function _diff_t!_impl(T::Type{<:Union{CompiledHomotopy,CompiledSystem}}, d, DP)
    H = interpret(T)
    checks, var_map = boundscheck_var_map(H)

    list, ids = instruction_list(H.expressions)

    vars = Symbol.(H.variables)
    params = Symbol.(H.parameters)


    diff_map = DiffMap()
    for (i, v) in enumerate(vars)
        for k = 1:(d-1)
            diff_map[v, k] = :(dx[$k][$i])
        end
    end

    for (i, v) in enumerate(params)
        for k = 1:DP
            diff_map[v, k] = :(dp[$k][$i])
        end
    end

    if H isa Homotopy
        diff_map[Symbol(H.t), 1] = 1
    end
    dlist = univariate_diff!(list, d, diff_map)

    assignements = Dict{Symbol,Vector{Expr}}()
    u_constants = Expr[]
    for (i, id) in enumerate(ids)
        d_id = diff_map[id, d]
        if d_id isa Symbol
            if d_id ∉ vars && d_id ∉ params
                add_assignement!(assignements, d_id, :(u[$i] = $d_id))
            else
                push!(u_constants, :(u[$i] = $(var_map[d_id])))
            end
        elseif d_id isa Nothing
            push!(u_constants, :(u[$i] = zero(eltype(x))))
        else
            push!(u_constants, :(u[$i] = $d_id))
        end
    end
    slp = to_expr(dlist, var_map, assignements)
    append!(slp.args, u_constants)

    quote
        $checks
        @inbounds $slp
        u
    end
end


################
## EVALUATION ##
################

# inplace (generated)
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

"""
    jacobian!(U, T::CompiledHomotopy, x, p = nothing)

Evaluate the Jacobian of `T` for variables `x`, `t` and parameters `p`
and store result in `u`.
"""
@generated function jacobian!(U, T::CompiledSystem, x, p = nothing)
    _jacobian!_impl(T)
end

"""
    jacobian!(U, T::CompiledHomotopy, x, t, p = nothing)

Evaluate the Jacobian of `T` for variables `x`, `t` and parameters `p` and
store result in `u`.
"""
@generated function jacobian!(U, T::CompiledHomotopy, x, t, p = nothing)
    _jacobian!_impl(T)
end

"""
    evaluate_and_jacobian!(u, U, T::CompiledHomotopy, x, p = nothing)

Evaluate `T` and its Jacobian for variables `x` and parameters `p` and
store result in `u`.
"""
@generated function evaluate_and_jacobian!(
    u,
    U,
    T::CompiledSystem,
    x,
    p = nothing,
)
    _evaluate_and_jacobian!_impl(T)
end

"""
    evaluate_and_jacobian!(u, U, T::CompiledHomotopy, x, t, p = nothing)

Evaluate `T` and its Jacobian for variables `x`, `t` and parameters `p` and
store result in `u`.
"""
@generated function evaluate_and_jacobian!(
    u,
    U,
    T::CompiledHomotopy,
    x,
    t,
    p = nothing,
)
    _evaluate_and_jacobian!_impl(T)
end

"""
    diff_t!(u, F::CompiledSystem, x, (x₁,…,xᵣ₋₁) = (), p = nothing, (p₁,…,pⱼ) = ())

Evaluate the expression
```math
(1 / r!) dʳ/dλʳ \\, F(x + ∑_{k=1}^{r-1} xᵢλⁱ, p(t + λ))
```
at ``λ = 0``.
If ``j < r-1`` then ``pₐ = 0`` for ``a > j``.
"""
@generated function diff_t!(
    u,
    T::CompiledSystem,
    x::AbstractVector,
    dx::NTuple{D,<:AbstractVector} = (),
    p::Union{Nothing,AbstractVector} = nothing,
    dp::NTuple{DP,<:AbstractVector} = (),
) where {D,DP}
    _diff_t!_impl(T, D + 1, DP)
end

"""
    diff_t!(u, H::CompiledHomotopy, x, t, (x₁,…,xᵣ₋₁) = (), p = nothing, (p₁,…,pⱼ) = ())

Evaluate the expression
```math
(1 / r!) dʳ/dλʳ \\, H(x + ∑_{k=1}^{r-1} xᵢλⁱ,t + λ)
```
at ``λ = 0``.
If ``j < r-1`` then ``pₐ = 0`` for ``a > j``.
"""
@generated function diff_t!(
    u,
    T::CompiledHomotopy,
    x::AbstractVector,
    t,
    dx::NTuple{D,<:AbstractVector} = (),
    p::Union{Nothing,AbstractVector} = nothing,
    dp::NTuple{DP,<:AbstractVector} = (),
) where {D,DP}
    _diff_t!_impl(T, D + 1, DP)
end

# non-inplace

"""
    to_smallest_eltype(A::AbstractArray)

Convert an array to the smallest eltype such that all elements still fit.

## Example
```julia
typeof(to_smallest_elype(Any[2,3])) == Vector{Int}
```
"""
function to_smallest_eltype(A::AbstractArray)
    T = typeof(first(A))
    for a in A
        T = promote_type(T, typeof(a))
    end
    convert.(T, A)
end

function evaluate(T::CompiledSystem, x, p = nothing)
    to_smallest_eltype(evaluate!(Vector{Any}(undef, size(T, 1)), T, x, p))
end
function evaluate(T::CompiledHomotopy, x, t, p = nothing)
    to_smallest_eltype(evaluate!(Vector{Any}(undef, size(T, 1)), T, x, t, p))
end

(T::CompiledSystem)(x, p = nothing) = evaluate(x, p)
(T::CompiledHomotopy)(x, t, p = nothing) = evaluate(x, t, p)

function jacobian(T::CompiledSystem, x, p = nothing)
    n, m = size(T)
    U = Matrix{Any}(undef, n, m)
    to_smallest_eltype(jacobian!(U, T, x, p))
end
function jacobian(T::CompiledHomotopy, x, t, p = nothing)
    n, m = size(T)
    U = Matrix{Any}(undef, n, m)
    to_smallest_eltype(jacobian!(U, T, x, t, p))
end

function evaluate_and_jacobian(T::CompiledSystem, x, p = nothing)
    n, m = size(T)
    u = Vector{Any}(undef, n)
    U = Matrix{Any}(undef, n, m)
    evaluate_and_jacobian!(u, U, T, x, p)
    to_smallest_eltype(u), to_smallest_eltype(U)
end
function evaluate_and_jacobian(T::CompiledHomotopy, x, t, p = nothing)
    n, m = size(T)
    u = Vector{Any}(undef, n)
    U = Matrix{Any}(undef, n, m)
    evaluate_and_jacobian!(u, U, T, x, t, p)
    to_smallest_eltype(u), to_smallest_eltype(U)
end

function diff_t(H::CompiledSystem, x::AbstractVector, args...)
    u = Vector{Any}(undef, size(H, 1))
    to_smallest_eltype(diff_t!(u, H, x, args...))
end
function diff_t(H::CompiledHomotopy, x::AbstractVector, t, args...)
    u = Vector{Any}(undef, size(H, 1))
    to_smallest_eltype(diff_t!(u, H, x, t, args...))
end
