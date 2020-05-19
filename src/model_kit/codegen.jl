"""
    TaylorVector{N,T}

A data structure representing a vector of Taylor series with `N` terms.
Each element is an `NTuple{N,T}`.

    TaylorVector(T, N::Integer, n::Integer)

Create a vector of `n` `NTuple{N,T}`s.
"""
struct TaylorVector{N,T} <: AbstractVector{NTuple{N,T}}
    data::LinearAlgebra.Transpose{T,Matrix{T}}
    views::NTuple{N,SubArray{T,1,Matrix{T},Tuple{Int,UnitRange{Int}},true}}
end
TaylorVector{N}(TV::TaylorVector{M,T}) where {N,M,T} =
    TaylorVector{N,T}(TV.data, TV.views[1:N])
function TaylorVector{N}(data::Matrix) where {N}
    views = tuple((view(data, i, 1:size(data, 2)) for i = 1:N)...)
    TaylorVector(LinearAlgebra.transpose(data), views)
end
function TaylorVector{N}(T, n::Integer) where {N}
    TaylorVector{N}(zeros(T, N, n))
end

function Base.show(io::IO, ::MIME"text/plain", X::TaylorVector)
    summary(io, X)
    isempty(X) && return
    println(io, ":")
    Base.print_array(io, map(i -> X[i], 1:length(X)))
end

Base.length(TV::TaylorVector) = size(TV.data, 1)
Base.size(TV::TaylorVector) = (length(TV),)
Base.eltype(::Type{TaylorVector{N,T}}) where {N,T} = NTuple{N,T}
Base.IndexStyle(::TaylorVector) = Base.IndexLinear()

"""
    vectors(TV::TaylorVec{N})

Return the Taylor series as `N` seperate vectors.
"""
vectors(TV::TaylorVector) = TV.views

@generated function Base.getindex(TV::TaylorVector{N}, i::Integer) where {N}
    quote
        Base.@_propagate_inbounds_meta
        x = TV.data
        $(Expr(:tuple, (:(x[i, $k]) for k = 1:N)...))
    end
end
Base.getindex(TV::TaylorVector, i::Integer, j::Integer) = getindex(TV.data, i, j)

function Base.setindex!(TV::TaylorVector{N,T}, x, i::Integer) where {N,T}
    setindex!(TV, convert(NTuple{N,T}, x), i)
end
function Base.setindex!(TV::TaylorVector{1}, x::Number, i::Integer)
    TV.data[i] = x
    x
end
@generated function Base.setindex!(
    TV::TaylorVector{N,T},
    x::NTuple{N,T},
    i::Integer,
) where {N,T}
    quote
        Base.@_propagate_inbounds_meta
        d = TV.data
        $(Expr(:tuple, (:(d[i, $k]) for k = 1:N)...)) = x
        x
    end
end
Base.setindex!(TV::TaylorVector, x, i::Integer, j::Integer) = setindex!(TV.data, x, i, j)

#####################
## CompiledSystem  ##
#####################
const TSYSTEM_TABLE = Dict{UInt,Vector{System}}()

struct CompiledSystem{HI}
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

interpret(TS::CompiledSystem) = TS.system
interpret(::Type{CompiledSystem{HI}}) where {HI} = TSYSTEM_TABLE[first(HI)][last(HI)]

Base.size(CS::CompiledSystem) = (CS.nexpressions, CS.nvariables)
Base.size(CS::CompiledSystem, i::Integer) = size(CS)[i]
Base.length(CS::CompiledSystem) = CS.nexpressions
nparameters(CS::CompiledSystem) = CS.nparameters
nvariables(CS::CompiledSystem) = CS.nvariables

######################
## CompiledHomotopy ##
######################

const THOMOTOPY_TABLE = Dict{UInt,Vector{Homotopy}}()

struct CompiledHomotopy{HI}
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

"""
    compile(F::System; optimizations = true)

Compile the given system. Returns a `CompiledSystem`. If `optimizations == true` then
the given system is optimized for evaluation efficiency.

    compile(H::Homotopy)

Compile the given homotopy. Returns a `CompiledHomotopy`.
"""
compile(F::System; optimizations::Bool = true) =
    CompiledSystem(optimizations ? optimize(F) : F)
compile(H::Homotopy) = CompiledHomotopy(H)

#############
## CODEGEN ##
#############

boundscheck_var_map(F::System; kwargs...) =
    boundscheck_var_map(F.expressions, F.variables, F.parameters; kwargs...)
boundscheck_var_map(H::Homotopy; kwargs...) =
    boundscheck_var_map(H.expressions, H.variables, H.parameters, H.t; kwargs...)
function boundscheck_var_map(
    exprs,
    vars,
    params,
    t = nothing;
    taylor = false,
    jacobian::Bool = false,
)
    n = length(exprs)
    m = length(vars)
    l = length(params)
    var_map = Dict{Symbol,Union{Symbol,Expr}}()
    for i = 1:m
        if taylor
            var_map[Symbol(vars[i])] = :(x[$i, 1])
        else
            var_map[Symbol(vars[i])] = :(x[$i])
        end
    end
    for i = 1:l
        if taylor
            var_map[Symbol(params[i])] = :(p[$i, 1])
        else
            var_map[Symbol(params[i])] = :(p[$i])
        end
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

function _evaluate!_impl(::Type{T}) where {T<:Union{CompiledSystem,CompiledHomotopy}}
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

function _jacobian!_impl(::Type{T}) where {T<:Union{CompiledSystem,CompiledHomotopy}}
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
                    add_assignement!(assignements, J[i, j], :(U[$i, $j] = $(J[i, j])))
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

function _evaluate_and_jacobian!_impl(
    ::Type{T},
) where {T<:Union{CompiledSystem,CompiledHomotopy}}
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
                    add_assignement!(assignements, J[i, j], :(U[$i, $j] = $(J[i, j])))
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

Base.@propagate_inbounds function set_row!(u::AbstractMatrix, t::Tuple{A}, i) where {A}
    u[i, 1] = first(t)
end
Base.@propagate_inbounds function set_row!(u::AbstractMatrix, t::Tuple{A,B}, i) where {A,B}
    a, b = t
    u[i, 1] = a
    u[i, 2] = b
end
Base.@propagate_inbounds function set_row!(
    u::AbstractMatrix,
    t::Tuple{A,B,C},
    i,
) where {A,B,C}
    a, b, c = t
    u[i, 1] = a
    u[i, 2] = b
    u[i, 3] = c
end
Base.@propagate_inbounds function set_row!(
    u::AbstractMatrix,
    t::Tuple{A,B,C,D},
    i,
) where {A,B,C,D}
    a, b, c, d = t
    u[i, 1] = a
    u[i, 2] = b
    u[i, 3] = c
    u[i, 4] = d
end
Base.@propagate_inbounds function set_row!(
    u::AbstractMatrix,
    t::Tuple{A,B,C,D,E},
    i,
) where {A,B,C,D,E}
    a, b, c, d, e = t
    u[i, 1] = a
    u[i, 2] = b
    u[i, 3] = c
    u[i, 4] = d
    u[i, 5] = e
end

function _functions_taylor!_impl(
    ::Type{T},
    K::Int;
    highest_order_only::Bool,
) where {T<:Union{CompiledSystem,CompiledHomotopy}}
    I = interpret(T)
    checks, var_map = boundscheck_var_map(I)
    list, ids = instruction_list(I.expressions)
    assignements = Dict{Symbol,Vector{Expr}}()
    for (i, id) in enumerate(ids)
        if highest_order_only
            add_assignement!(assignements, id, :(u[$i] = last($id)))
        else
            add_assignement!(assignements, id, :(set_row!(u, $id, $i)))
        end
    end

    block = Expr(:block)
    exprs = block.args
    for (id, (op, arg1, arg2)) in list.instructions
        a = get(var_map, arg1, arg1)
        if op == :^
            r::Int = arg2
            if r == 2
                push!(exprs, :($id = taylor(Val{:sqr}, Val{$K}, $a)))
            else
                push!(exprs, :($id = taylor(Val{:^}, Val{$K}, $a, $r)))
            end
        elseif arg2 !== nothing
            b = get(var_map, arg2, arg2)
            push!(exprs, :($id = taylor(Val{$(QuoteNode(op))}, Val{$K}, $a, $b)))
        else
            push!(exprs, :($id = taylor(Val{$(QuoteNode(op))}, Val{$K}, $a)))
        end
        if haskey(assignements, id)
            append!(exprs, assignements[id])
        end
    end
    # TODO: let block only if homotopy....
    if I isa Homotopy
        quote
            $checks
            let t = (t, one(t))
                @inbounds $block
            end
            u
        end
    else
        quote
            $checks
            @inbounds $block
            u
        end
    end
end

function _inline_taylor!_impl(
    T::Type{<:Union{CompiledHomotopy,CompiledSystem}},
    K,
    dx,
    dp;
    highest_order_only::Bool,
)
    H = interpret(T)
    # @show H
    checks, var_map = boundscheck_var_map(H; taylor = true)

    list, ids = instruction_list(H.expressions)

    vars = Symbol.(H.variables)
    params = Symbol.(H.parameters)
    # @show vars, dx, dp
    # @show params
    diff_map = DiffMap()
    for (i, v) in enumerate(vars)
        for k = 1:dx
            diff_map[v, k] = :(x[$i, $(k + 1)])
        end
    end

    for (i, v) in enumerate(params)
        for k = 1:dp
            diff_map[v, k] = :(p[$i, $(k + 1)])
        end
    end

    if H isa Homotopy
        diff_map[Symbol(H.t), 1] = 1
    end
    # @show diff_map
    dlist = univariate_diff!(list, K, diff_map)

    assignements = Dict{Symbol,Vector{Expr}}()
    u_constants = Expr[]

    if highest_order_only
        for (i, id) in enumerate(ids)
            k = K
            k_id = diff_map[id, k]
            if k_id isa Symbol
                if k_id ∉ vars && k_id ∉ params
                    add_assignement!(assignements, k_id, :(u[$i] = $k_id))
                else
                    push!(u_constants, :(u[$i] = $(var_map[k_id])))
                end
            elseif k_id isa Nothing
                push!(u_constants, :(u[$i] = zero(eltype(u))))
            else
                push!(u_constants, :(u[$i] = $k_id))
            end
        end
    else
        for (i, id) in enumerate(ids)
            add_assignement!(assignements, id, :(u[$i, 1] = $id))
            for k = 1:K
                k_id = diff_map[id, k]
                if k_id isa Symbol
                    if k_id ∉ vars && k_id ∉ params
                        add_assignement!(assignements, k_id, :(u[$i, $(k + 1)] = $k_id))
                    else
                        push!(u_constants, :(u[$i, $(k + 1)] = $(var_map[k_id])))
                    end
                elseif k_id isa Nothing
                    push!(u_constants, :(u[$i, $(k + 1)] = zero(u[$i, $(k + 1)])))
                else
                    push!(u_constants, :(u[$i, $(k + 1)] = $k_id))
                end
            end
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

function _taylor!_impl(T, K, D, DP; kwargs...)
    # Experiments show that for K < 2 (i.e. u and u̇) the fully inlined version has better
    # compilation times, for K ≥ 2 the function version starts to become
    # significantly faster
    if K ≤ 1
        _inline_taylor!_impl(T, K, D, DP; kwargs...)
    else
        _functions_taylor!_impl(T, K; kwargs...)
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
@generated function evaluate_and_jacobian!(u, U, T::CompiledSystem, x, p = nothing)
    _evaluate_and_jacobian!_impl(T)
end

"""
    evaluate_and_jacobian!(u, U, T::CompiledHomotopy, x, t, p = nothing)

Evaluate `T` and its Jacobian for variables `x`, `t` and parameters `p` and
store result in `u`.
"""
@generated function evaluate_and_jacobian!(u, U, T::CompiledHomotopy, x, t, p = nothing)
    _evaluate_and_jacobian!_impl(T)
end

"""
    taylor!(
        u::AbstractMatrix,
        v::Val{M}
        F::CompiledSystem,
        x::TaylorVector{D},
        p::Union{Nothing,Vector,TaylorVector} = nothing,
    )
    taylor!(
        u::TaylorVector{M},
        v::Val{M}
        F::CompiledSystem,
        x::TaylorVector{D},
        p::Union{Nothing,Vector,TaylorVector} = nothing,
    )

Compute the Taylor series of ``u = F(x,p)`` given the Taylor series `x` and `p`.
"""
@generated function taylor!(
    u::AbstractMatrix,
    v::Val{M},
    T::CompiledSystem,
    x::TaylorVector{D},
    p::Nothing = nothing,
) where {M,D}
    _taylor!_impl(T, M, D - 1, -1; highest_order_only = false)
end
@generated function taylor!(
    u::AbstractMatrix,
    v::Val{M},
    T::CompiledSystem,
    x::TaylorVector{D},
    p::AbstractVector,
) where {M,D}
    _taylor!_impl(T, M, D - 1, 0; highest_order_only = false)
end
@generated function taylor!(
    u::AbstractMatrix,
    v::Val{M},
    T::CompiledSystem,
    x::TaylorVector{D},
    p::TaylorVector{DP},
) where {M,D,DP}
    _taylor!_impl(T, M, D - 1, DP - 1; highest_order_only = false)
end

function taylor!(
    u::TaylorVector{M},
    T::CompiledSystem,
    x::TaylorVector{D},
    p = nothing,
) where {M,D}
    taylor!(u.data, Val(M - 1), T, x, p)
    u
end


@generated function taylor!(
    u::AbstractVector,
    ::Val{M},
    T::CompiledSystem,
    x::TaylorVector{D},
    p::Nothing = nothing,
) where {M,D}
    _taylor!_impl(T, M, D - 1, -1; highest_order_only = true)
end
@generated function taylor!(
    u::AbstractVector,
    ::Val{M},
    T::CompiledSystem,
    x::TaylorVector{D},
    p::AbstractVector,
) where {M,D}
    _taylor!_impl(T, M, D - 1, 0; highest_order_only = true)
end
@generated function taylor!(
    u::AbstractVector,
    ::Val{M},
    T::CompiledSystem,
    x::TaylorVector{D},
    p::TaylorVector{DP},
) where {M,D,DP}
    _taylor!_impl(T, M, D - 1, DP - 1; highest_order_only = true)
end

@generated function taylor!(
    u::AbstractVector,
    ::Val{M},
    T::CompiledSystem,
    x::AbstractVector,
    p::Nothing = nothing,
) where {M}
    _taylor!_impl(T, M, 0, -1; highest_order_only = true)
end
@generated function taylor!(
    u::AbstractVector,
    ::Val{M},
    T::CompiledSystem,
    x::AbstractVector,
    p::AbstractVector,
) where {M}
    _taylor!_impl(T, M, 0, 0; highest_order_only = true)
end
@generated function taylor!(
    u::AbstractVector,
    ::Val{M},
    T::CompiledSystem,
    x::AbstractVector,
    p::TaylorVector{DP},
) where {M,DP}
    _taylor!_impl(T, M, 0, DP - 1; highest_order_only = true)
end

"""
    taylor!(
        u::AbstractVector,
        ::Val{M}
        H::Homotopy,
        x::TaylorVector{D},
        t::Number,
        p::Union{Nothing,Vector,TaylorVector} = nothing,
    )

Compute the `M`-th derivative of ``u=H(x,t,p)`` given the taylor series `x` and `p`.
"""
@generated function taylor!(
    u::AbstractVector,
    ::Val{M},
    T::CompiledHomotopy,
    x::TaylorVector{D},
    t::Number,
    p::Nothing = nothing,
) where {M,D}
    _taylor!_impl(T, M, D - 1, -1; highest_order_only = true)
end
@generated function taylor!(
    u::AbstractVector,
    ::Val{M},
    T::CompiledHomotopy,
    x::TaylorVector{D},
    t::Number,
    p::AbstractVector,
) where {M,D}
    _taylor!_impl(T, M, D - 1, 0; highest_order_only = true)
end
@generated function taylor!(
    u::AbstractVector,
    ::Val{M},
    T::CompiledHomotopy,
    x::TaylorVector{D},
    t::Number,
    p::TaylorVector{DP},
) where {M,D,DP}
    _taylor!_impl(T, M, D - 1, DP - 1; highest_order_only = true)
end
@generated function taylor!(
    u::AbstractVector,
    ::Val{M},
    T::CompiledHomotopy,
    x::AbstractVector,
    t::Number,
    p::Nothing = nothing,
) where {M}
    _taylor!_impl(T, M, 0, -1; highest_order_only = true)
end
@generated function taylor!(
    u::AbstractVector,
    ::Val{M},
    T::CompiledHomotopy,
    x::AbstractVector,
    t::Number,
    p::AbstractVector,
) where {M}
    _taylor!_impl(T, M, 0, 0; highest_order_only = true)
end
@generated function taylor!(
    u::AbstractVector,
    ::Val{1},
    T::CompiledHomotopy,
    x::AbstractVector,
    t::Number,
    p::TaylorVector{DP},
) where {M,DP}
    _taylor!_impl(T, M, 0, DP - 1; highest_order_only = true)
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
