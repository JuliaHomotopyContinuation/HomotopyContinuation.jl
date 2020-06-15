
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

function _evaluate!_impl(::Type{T}) where {T}
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

function _jacobian!_impl(::Type{T}) where {T}
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

function _evaluate_and_jacobian!_impl(::Type{T}) where {T}
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

function _functions_taylor!_impl(::Type{T}, K::Int; highest_order_only::Bool) where {T}
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
    T::Type,
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
