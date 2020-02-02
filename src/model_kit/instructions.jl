struct InstructionList
    instructions::Vector{Tuple{Symbol,Tuple{Symbol,Any,Any}}}
    var::Symbol
    n::Base.RefValue{Int}
end

function InstructionList(; var::Symbol = :_ι, n::Base.RefValue{Int} = Ref(0))
    InstructionList(Tuple{Symbol,Tuple{Symbol,Any,Any}}[], var, n)
end

function instruction_list(
    exprs::Vector{Expression};
    subexpressions::Dict{Expression,Expression} = Dict{Expression,Expression}(),
    perform_cse::Bool = true,
)
    v = InstructionList()
    PSE = Set{Symbol}()
    if perform_cse
        exprs, CSE = cse(exprs)
        merge!(CSE, subexpressions)
    else
        CSE = subexpressions
    end
    v, map(ex -> flat_expr!(v, ex, CSE, PSE), exprs)
end

Base.length(v::InstructionList) = v.n[]
function Base.show(io::IO, v::InstructionList)
    for (id, arg) in v.instructions
        println(io, :($id = $(Expr(:call, arg...))))
    end
end

function Base.push!(v::InstructionList, x)
    id = Symbol(v.var, v.n[] += 1)
    push!(v.instructions, (id, x))
    id
end
function Base.push!(v::InstructionList, x::Pair{Symbol,<:Tuple{Symbol,Any,Any}})
    push!(v.instructions, (first(x), last(x)))
    first(x)
end

function flat_expr!(
    v,
    ex::Basic,
    # common sub expression
    CSE::Dict{Expression,Expression},
    # processed sub expression
    PSE::Set{Symbol} = Set{Symbol}(),
)
    t = class(ex)
    if t == :Symbol
        s = Symbol(to_string(ex))
        if haskey(CSE, ex) && s ∉ PSE
            push!(PSE, s)
            val = flat_expr!(v, CSE[ex], CSE, PSE)
            if val isa Symbol
                v.instructions[end] = (s, last(v.instructions[end]))
            else
                push!(v, s => val)
            end
        end
        return s
    elseif t == :Pow
        x, k = ModelKit.args(ex)
        return push!(v, (:^, flat_expr!(v, x, CSE, PSE), convert(Int, k)))
    elseif t == :Add || t == :Mul
        op = t == :Add ? :+ : :*
        vec = ModelKit.args(ex)
        op_arg = flat_expr!(v, vec[1], CSE, PSE)
        for i = 2:length(vec)
            op_arg = push!(v, (op, op_arg, flat_expr!(v, vec[i], CSE, PSE)))
        end
        return op_arg
    elseif t == :Div
        x, y = ModelKit.args(ex)
        return push!(
            v,
            (:/, flat_expr!(v, x, CSE, PSE), flat_expr!(v, y, CSE, PSE)),
        )
    else
        return to_number(ex)
    end
end


function Base.diff(
    list::InstructionList,
    vars::Vector{Symbol},
    f::Vector{Symbol},
)
    seed = Dict{Tuple{Symbol,Int},Any}()
    for (i, v) in enumerate(vars)
        seed[(v, i)] = 1
    end
    dlist = diff!(list, length(vars), seed)

    J = Matrix{Union{Nothing,Symbol,Number}}(undef, length(f), length(vars))
    for (j, v) in enumerate(vars), (i, fi) in enumerate(f)
        if haskey(seed, (fi, j))
            J[i, j] = seed[(fi, j)]
        else
            J[i, j] = nothing
        end
    end
    dlist, J
end

"""
    diff!(list::InstructionList, N::Int, diff_map)

Returns an `InstructionList` list computing the derivatives of `N` variables.
The derivatives of the instructions are stored in `diff_map`.
"""
function diff!(list::InstructionList, N::Int, diff_map)
    n = length(list)
    v = InstructionList(n = Ref(n))
    for (id, el) in list.instructions
        (op, arg1, arg2) = el

        if op == :^
            let
                p1 = p2 = :NONE
                instr_added = false
                for ∂i = 1:N
                    exp::Int = arg2
                    if haskey(diff_map, (arg1, ∂i))
                        if p2 == :NONE
                            if exp == 2
                                p2 = push!(v, (:*, 2, arg1))
                            else
                                p1 = push!(v, (:^, arg1, exp - 1))
                                p2 = push!(v, (:*, exp, p1))
                            end
                        end
                        if !instr_added
                            if exp == 2
                                push!(v, id => (:^, arg1, 2))
                            else
                                push!(v, id => (:*, p1, arg1))
                            end
                            instr_added = true
                        end
                        ∂el = diff_map[(arg1, ∂i)]
                        if ∂el != 1
                            diff_map[(id, ∂i)] = push!(v, (:*, p2, ∂el))
                        else
                            diff_map[(id, ∂i)] = p2
                        end
                    elseif p2 != :NONE && !instr_added
                        push!(v, id => el)
                        instr_added = true
                    end

                end
                if !instr_added
                    push!(v, id => el)
                end
            end
        elseif op == :*
            let
                for ∂i = 1:N
                    ∂i == 1 && push!(v, id => el)

                    has_∂1 = haskey(diff_map, (arg1, ∂i))
                    has_∂2 = haskey(diff_map, (arg2, ∂i))

                    if has_∂2
                        a2::Symbol = arg2
                        ∂arg2 = diff_map[(a2, ∂i)]
                        if ∂arg2 != 1
                            e1 = push!(v, (:*, arg1, ∂arg2))
                        else
                            e1 = arg1
                        end
                    end

                    if has_∂1
                        a1::Symbol = arg1
                        ∂arg1 = diff_map[(a1, ∂i)]
                        if ∂arg1 != 1
                            e2 = push!(v, (:*, ∂arg1, arg2))
                        else
                            e2 = arg2
                        end
                    end

                    if has_∂1 && has_∂2
                        diff_map[(id, ∂i)] = push!(v, (:+, e1, e2))
                    elseif has_∂1
                        diff_map[(id, ∂i)] = e2
                    elseif has_∂2
                        diff_map[(id, ∂i)] = e1
                    end
                end
            end
        elseif op == :/
            let
                a1::Symbol
                a2::Symbol
                for ∂i = 1:N
                    ∂i == 1 && push!(v, id => el)

                    has_∂1 = haskey(diff_map, (arg1, ∂i))
                    has_∂2 = haskey(diff_map, (arg2, ∂i))

                    if has_∂1 && has_∂2
                        a1 = arg1
                        a2 = arg2
                        ∂arg1 = diff_map[(a1, ∂i)]
                        ∂arg2 = diff_map[(a2, ∂i)]
                        e1 = push!(v, (:*, ∂arg1, a2))
                        e2 = push!(v, (:*, a1, ∂arg2))
                        e3 = push!(v, (:-, e1, e2))
                        e4 = push!(v, :(:^, arg2, 2))
                        diff_map[(id, ∂i)] = push!(v, :(:/, e3, e4))
                    elseif has_∂1
                        a1 = arg1
                        ∂arg1 = diff_map[(a1, ∂i)]
                        diff_map[(id, ∂i)] = push!(v, (:/, ∂arg1, arg2))
                    elseif has_∂2
                        a2 = arg2
                        ∂arg2 = diff_map[(a2, ∂i)]
                        e1 = push!(v, (:*, arg1, ∂arg2))
                        e2 = push!(v, (:*, -1, e1))
                        e3 = push!(v, :(:^, arg2, 2))
                        diff_map[(id, ∂i)] = push!(v, :(:/, e2, e3))
                    end
                end
            end
        elseif op == :+
            let
                for ∂i = 1:N
                    ∂i == 1 && push!(v, id => el)

                    has_∂1 = haskey(diff_map, (arg1, ∂i))
                    has_∂2 = haskey(diff_map, (arg2, ∂i))

                    if has_∂1 && has_∂2
                        diff_map[(id, ∂i)] = push!(
                            v,
                            (:+, diff_map[(arg1, ∂i)], diff_map[(arg2, ∂i)]),
                        )
                    elseif has_∂1
                        diff_map[(id, ∂i)] = diff_map[(arg1, ∂i)]
                    elseif has_∂2
                        diff_map[(id, ∂i)] = diff_map[(arg2, ∂i)]
                    end
                end
            end
        end
    end
    v
end

## Higher order AD. This follows Chapter 13 of Griewank and Walther
"""
    univariate_diff!(list::InstructionList, K::Int, diff_map)

Returns an `InstructionList` list computing the first `K` derivatives.
The derivatives of the instructions are stored in `diff_map`.
"""
function univariate_diff!(list::InstructionList, K::Int, diff_map)
    n = length(list)
    v = InstructionList(n = Ref(n))
    for (id, el) in list.instructions
        (op, arg1, arg2) = el

        if op == :^
            r::Int = arg2
            if haskey(diff_map, (arg1, 1)) && r == 2
                push!(v, id => el)
                for k = 2:K
                    w_k = nothing
                    for j = 0:k
                        u_j = j == 0 ? arg1 : get(diff_map, (arg1, j), nothing)
                        u_kj = j == k ? arg1 :
                               get(diff_map, (arg1, k - j), nothing)
                        if j < k - j && u_j !== nothing && u_kj !== nothing
                            s = push!(v, (:*, u_j, u_kj))
                            if w_k === nothing
                                w_k = s
                            else
                                w_k = push!(v, (:+, w_k, s))
                            end
                        end
                        if j ≥ k - j && w_k !== nothing
                            w_k = push!(v, (:*, 2, w_k))
                        end
                        if j == k - j
                            s = push!(v, (:*, u_j, u_kj))
                            if w_k === nothing
                                w_k = s
                            else
                                w_k = push!(v, (:+, w_k, s))
                            end
                        end
                        if j ≥ k - j
                            break
                        end
                    end
                    if w_k !== nothing
                        diff_map[(id, k)] = w_k
                    end
                end
            elseif haskey(diff_map, (arg1, 1)) && r > 2
                w_0 = push!(v, id => (:^, arg1, r))
                u_0_inv = push!(v, (:inv_not_zero, arg1, nothing))
                # w_k is v_k in Griewank
                w = Symbol[]
                for k = 1:K
                    s = nothing
                    for j = 1:k
                        u_j = get(diff_map, (arg1, j), nothing)
                        if u_j !== nothing
                            ũ_j = j == 1 ? u_j : push!(v, (:*, j, u_j))
                        else
                            ũ_j = nothing
                        end

                        w_kj = k == j ? w_0 : w[k-j]
                        if ũ_j == 1
                            s_j = w_kj
                        elseif ũ_j !== nothing
                            s_j = push!(v, (:*, w_kj, ũ_j))
                        else
                            s_j = nothing
                        end
                        if s_j !== nothing
                            s = isnothing(s) ? s_j : push!(v, (:+, s, s_j))
                        end
                    end
                    if s !== nothing
                        s = push!(v, (:*, r, s))
                    end

                    t = nothing
                    for j = 1:k-1
                        u_kj = get(diff_map, (arg1, k - j), nothing)
                        w̃_j = j == 1 ? w[1] : push!(v, (:*, j, w[j]))

                        if u_kj == 1
                            t_j = w̃_j
                        elseif u_kj !== nothing
                            t_j = push!(v, (:*, u_kj, w̃_j))
                        else
                            t_j = nothing
                        end

                        if t_j !== nothing
                            t = isnothing(t) ? t_j : push!(v, (:+, t, t_j))
                        end
                    end
                    if t == nothing && s == nothing
                        continue
                    elseif t == nothing
                        w̃_k = push!(v, (:*, u_0_inv, s))
                    else # s cannot be nothing
                        w̃_k = push!(v, (:*, u_0_inv, push!(v, (:-, s, t))))
                    end

                    if k == 1
                        w_k = w̃_k
                    else
                        w_k = push!(v, (:/, w̃_k, k))
                    end

                    push!(w, w_k)
                    diff_map[(id, k)] = w_k
                end
            else
                push!(v, id => el)
            end
        elseif op == :*
            push!(v, id => el)

            for k = 1:K
                c_k = nothing
                for j = 0:k
                    if j == 0
                        a = arg1
                        b = get(diff_map, (arg2, k - j), nothing)
                    elseif j == k
                        a = get(diff_map, (arg1, k), nothing)
                        b = arg2
                    else
                        a = get(diff_map, (arg1, j), nothing)
                        b = get(diff_map, (arg2, k - j), nothing)
                    end

                    if a === nothing || b === nothing
                        continue
                    end

                    if c_k === nothing
                        c_k = push!(v, (:*, a, b))
                    else
                        c_k = push!(v, (:+, c_k, push!(v, (:*, a, b))))
                    end
                end

                if c_k === nothing
                    break
                else
                    diff_map[(id, k)] = c_k
                end
            end
        elseif op == :+
            push!(v, id => el)
            for k = 1:K
                a = get(diff_map, (arg1, k), nothing)
                b = get(diff_map, (arg2, k), nothing)
                if a !== nothing && b !== nothing
                    diff_map[(id, k)] = push!(v, (:+, a, b))
                elseif a !== nothing
                    diff_map[(id, k)] = a
                elseif b !== nothing
                    diff_map[(id, k)] = b
                end
            end
        elseif op == :/
            error("Not implemented")
        end
    end
    v
end

@inline sqr(x) = x^2
@inline function sqr(z::Complex)
    x, y = reim(z)
    Complex((x + y) * (x - y), (x + x) * y)
end

"""
    inv_not_zero(x)

Invert x unless it is 0, then return 0.
"""
@inline inv_not_zero(x) = ifelse(iszero(x), x, inv(x))

function unroll_pow(var, n)
    n == 0 && return :(one($var))
    n == 1 && return var
    n == 2 && return :(sqr($var))
    n == 3 && return :($var * sqr($var))
    n == 4 && return :(sqr(sqr($var)))
    n == 5 && return :($var * sqr(sqr($var)))

    # base to expansion shows up which power it is needed to compute
    d = digits(n, base = 2)
    exprs = map(2:length(d)) do i
        :($(Symbol(:x, 1 << (i - 1))) = sqr($(Symbol(:x, 1 << (i - 2)))))
    end
    prods = Symbol[]
    for (i, di) in enumerate(d)
        if !iszero(di)
            push!(prods, Symbol(:x, 1 << (i - 1)))
        end
    end
    if length(prods) > 1
        push!(exprs, Expr(:call, :*, prods...))
    end
    Expr(:let, :(x1 = $var), Expr(:block, exprs...))
end


function to_expr(
    list::InstructionList,
    var_map = Dict{Symbol,Union{Expr,Symbol}}(),
    assignements = Dict{Symbol,Expr}(),
)
    exprs = Expr[]
    for (id, (op, arg1, arg2)) in list.instructions
        if op == :^
            x::Symbol = arg1
            k::Int = arg2
            if k < 0
                push!(
                    exprs,
                    :($id = inv($(unroll_pow(get(var_map, x, x), -k)))),
                )
            else
                push!(exprs, :($id = $(unroll_pow(get(var_map, x, x), k))))
            end
        else
            a = get(var_map, arg1, arg1)
            if arg2 !== nothing
                b = get(var_map, arg2, arg2)
                push!(exprs, :($id = $(Expr(:call, op, a, b))))
            else
                push!(exprs, :($id = $(Expr(:call, op, a))))
            end
        end

        if haskey(assignements, id)
            push!(exprs, assignements[id])
        end
    end
    Expr(:block, exprs...)
end
