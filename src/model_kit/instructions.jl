@enum InstructionOp::UInt8 begin
    INSTR_ADD
    INSTR_SUB
    INSTR_NEG
    INSTR_MUL
    INSTR_DIV
    INSTR_MULADD
    INSTR_MULSUB
    INSTR_SUBMUL #-x*y+z
    INSTR_SQR
    INSTR_POW
    INSTR_SIN
    INSTR_COS
end

function call_op(op::InstructionOp)
    if op == INSTR_ADD
        :add_fast
    elseif op == INSTR_SUB
        :sub_fast
    elseif op == INSTR_NEG
        :neg_fast
    elseif op == INSTR_MUL
        :mul_fast
    elseif op == INSTR_DIV
        :div_fast
    elseif op == INSTR_MULADD
        :muladd_fast
    elseif op == INSTR_MULSUB
        :mulsub_fast
    elseif op == INSTR_SUBMUL
        :submul_fast
    elseif op == INSTR_SQR
        :sqr_fast
    elseif op == INSTR_POW
        :pow_fast
    elseif op == INSTR_SIN
        :sin
    elseif op == INSTR_COS
        :cos
    else
        error(op)
        :UNKNOWN
    end
end

function taylor_call_op(op::InstructionOp)
    if op == INSTR_ADD
        :taylor_add
    elseif op == INSTR_SUB
        :taylor_sub
    elseif op == INSTR_NEG
        :tayor_neg
    elseif op == INSTR_MUL
        :taylor_mul
    elseif op == INSTR_DIV
        :taylor_div
    elseif op == INSTR_MULADD
        :taylor_muladd
    elseif op == INSTR_MULSUB
        :taylor_mulsub
    elseif op == INSTR_SUBMUL
        :taylor_submul
    elseif op == INSTR_SQR
        :taylor_sqr
    elseif op == INSTR_POW
        :taylor_pow
    elseif op == INSTR_SIN
        :taylor_sin
    elseif op == INSTR_COS
        :taylor_cos
    else
        error(op)
        :UNKNOWN
    end
end

function is_univariate(op::InstructionOp)
    op == INSTR_SQR || op == INSTR_NEG
end

@inline neg_fast(x) = Base.FastMath.sub_fast(x)
@inline add_fast(x, y) = Base.FastMath.add_fast(x, y)
@inline sub_fast(x, y) = Base.FastMath.sub_fast(x, y)
@inline mul_fast(x, y) = Base.FastMath.mul_fast(x, y)
@inline div_fast(x, y) = Base.FastMath.div_fast(x, y)
@inline pow_fast(x, p::Integer) = Base.power_by_squaring(x, p)
@inline muladd_fast(x, y, z) = @fastmath x * y + z
@inline mulsub_fast(x, y, z) = @fastmath x * y - z
@inline submul_fast(x, y, z) = @fastmath z - x * y
@inline sqr_fast(x) = sqr(x)
@inline sqr(x) = @fastmath x * x
@inline function sqr(z::Complex)
    x, y = reim(z)
    @fastmath Complex((x + y) * (x - y), (x + x) * y)
end
"""
    inv_not_zero(x)

Invert x unless it is 0, then return 0.
"""
@inline inv_not_zero(x) = ifelse(is_zero(x), x, @fastmath inv(x))
@inline inv_fast(x) = @fastmath inv(x)

struct InstructionRef
    i::Int
end
Base.show(io::IO, i::InstructionRef) = print(io, "ι_", i.i)

const InstructionArg = Union{Nothing,Number,Symbol,InstructionRef}

struct Instruction
    op::InstructionOp
    args::Tuple{InstructionArg,InstructionArg,InstructionArg}
end

Instruction(op::InstructionOp, a::T) where {T>:Tuple} =
    Instruction(op, (a, nothing, nothing))
Instruction(op, a, b) = Instruction(op, (a, b, nothing))
Instruction(op, a, b, c) = Instruction(op, (a, b, c))
Base.getindex(I::Instruction, i) = I.args[i]


struct InstructionList
    instructions::Vector{Instruction}
    lifetime::Vector{Union{Nothing,Int}}
end

InstructionList() = InstructionList(Instruction[], Union{Nothing,Int}[])
Base.length(v::InstructionList) = length(v.instructions)
function Base.iterate(I::InstructionList, state = nothing)
    if isnothing(state)
        next_state = iterate(I.instructions)
    else
        next_state = iterate(I.instructions, state)
    end
    if isnothing(next_state)
        return nothing
    end
    next, state = next_state
    (InstructionRef(state - 1), next), state
end


function instruction_list(
    exprs::Vector{Expression};
    subexpressions::Dict{Expression,Expression} = Dict{Expression,Expression}(),
    perform_cse::Bool = true,
)
    list = InstructionList()
    if perform_cse
        exprs, CSE = cse(exprs)
        merge!(CSE, subexpressions)
    else
        CSE = subexpressions
    end
    symb_cse = Dict{Symbol,Expression}()
    for (k, v) in CSE
        push!(symb_cse, Symbol(to_string(k)) => v)
    end
    pse = Dict{Symbol,InstructionRef}()
    assignments = map(ex -> add_instructions!(list, ex, symb_cse, pse), exprs)

    # makesure to not overwrite assignments
    for iref in assignments
        update_lifetime!(list, iref, typemax(Int64))
    end
    list, assignments
end


function Base.show(io::IO, list::InstructionList)
    for (id, instr) in list
        println(
            io,
            :($id = $(Expr(:call, Symbol(instr.op), filter(!isnothing, instr.args)...))),
        )
    end
end

function update_lifetime!(list::InstructionList, x::InstructionRef, n::Int)
    list.lifetime[x.i] = n
    nothing
end
update_lifetime!(list::InstructionList, x::InstructionArg, n::Int) = nothing


function Base.push!(list::InstructionList, x::Instruction)
    push!(list.instructions, x)
    n = length(list)
    # by default does an element have infinite lifetime
    push!(list.lifetime, typemax(Int64))

    (a, b, c) = x.args
    update_lifetime!(list, a, n)
    update_lifetime!(list, b, n)
    update_lifetime!(list, c, n)

    InstructionRef(n)
end

# Helpers for InstructionOps
function add!(list::InstructionList, @nospecialize(a), @nospecialize(b))
    if a === nothing
        b === nothing ? nothing : b
    else
        b === nothing ? a : push!(list, Instruction(INSTR_ADD, a, b))
    end
end
function neg!(list::InstructionList, @nospecialize(a))
    a === nothing ? nothing : push!(list, Instruction(INSTR_NEG, a))
end
function sub!(list::InstructionList, @nospecialize(a), @nospecialize(b))
    if a === nothing
        b === nothing ? nothing : push!(list, Instruction(INSTR_NEG, b))
    else
        b === nothing ? a : push!(list, Instruction(INSTR_SUB, a, b))
    end
end

function mul!(list::InstructionList, @nospecialize(a), @nospecialize(b))
    (a === nothing || b === nothing) && return nothing
    if is_one(a)
        is_one(b) ? a : b
    elseif is_minus_one(a)
        push!(list, Instruction(INSTR_NEG, b))
    elseif is_minus_one(b)
        push!(list, Instruction(INSTR_NEG, a))
    elseif a isa Integer && a == 2
        push!(list, Instruction(INSTR_ADD, b, b))
    else
        is_one(b) ? a : push!(list, Instruction(INSTR_MUL, a, b))
    end
end

function muladd!(
    list::InstructionList,
    @nospecialize(a),
    @nospecialize(b),
    @nospecialize(c),
)
    (a === nothing || b === nothing) && return c
    c === nothing && return mul!(list, a, b)
    is_one(a) && return add!(list, b, c)
    is_one(b) && return add!(list, a, c)
    push!(list, Instruction(INSTR_MULADD, a, b, c))
end
function mulsub!(
    list::InstructionList,
    @nospecialize(a),
    @nospecialize(b),
    @nospecialize(c),
)
    (a === nothing || b === nothing) && return neg!(list, c)
    c === nothing && return mul!(list, a, b)
    is_one(a) && return sub!(list, b, c)
    is_one(b) && return sub!(list, a, c)
    push!(list, Instruction(INSTR_MULSUB, a, b, c))
end
function submul!(
    list::InstructionList,
    @nospecialize(a),
    @nospecialize(b),
    @nospecialize(c),
)
    # c - a * b
    (a === nothing || b === nothing) && return c
    c === nothing && return neg!(list, mul!(list, a, b))
    is_one(a) && return sub!(list, c, b)
    is_one(b) && return sub!(list, c, a)
    push!(list, Instruction(INSTR_SUBMUL, a, b, c))
end

function div!(list::InstructionList, @nospecialize(a), @nospecialize(b))
    a === nothing && return nothing
    is_one(b) ? a : push!(list, Instruction(INSTR_DIV, a, b))
end

function sqr!(list, @nospecialize(a))
    a === nothing && return nothing
    push!(list, Instruction(INSTR_SQR, a))
end
function pow!(list, @nospecialize(a), b::Integer)
    @assert b ≥ 0
    is_one(b) && return a
    b == 2 && return sqr!(list, a)
    is_zero(b) && return Int32(1)
    is_one(a) && return a
    a === nothing ? nothing : push!(list, Instruction(INSTR_POW, a, b))
end

function add_instructions!(
    list::InstructionList,
    ex::Basic,
    # common sub expression
    cse::Dict{Symbol,Expression},
    # processed sub expression
    # mapping of the cse var with the instruction used
    pse::Dict{Symbol,InstructionRef},
)
    t = class(ex)
    if t == :Symbol
        s = Symbol(to_string(ex))
        if haskey(pse, s)
            return pse[s]
        elseif haskey(cse, s)
            val = add_instructions!(list, cse[s], cse, pse)
            if val isa InstructionRef
                push!(pse, s => val)
            end
            val
        else
            s
        end
    elseif t == :Pow
        x, k_expr = args(ex)
        k = convert(Int32, k_expr)
        xinstr = add_instructions!(list, x, cse, pse)
        if k == -1
            div!(list, Int32(1), xinstr)
        elseif k < -1
            div!(list, Int32(1), pow!(list, xinstr, -k))
        else
            pow!(list, xinstr, k)
        end
    elseif t == :Mul
        _add_mul_instructions!(list, args(ex), cse, pse)
    elseif t == :Add
        # For sums we do work to find muladd, mulsub and negmuladd instructions
        vec = args(ex)
        neg = false
        op_arg = nothing
        for i = 1:length(vec)
            vᵢ = vec[i]
            if class(vᵢ) == :Mul
                vᵢ_args = args(vᵢ)
                nᵢ = length(vᵢ_args)
                if is_minus_one(vᵢ_args[1])
                    x = add_instructions!(list, vᵢ_args[2], cse, pse)
                    if length(vᵢ_args) > 2
                        y = _add_mul_instructions!(list, view(vᵢ_args, 3:nᵢ), cse, pse)
                        if isnothing(op_arg)
                            op_arg = mul!(list, x, y)
                            neg = true
                        elseif neg
                            op_arg = muladd!(list, x, y, op_arg)
                            neg = true
                        else
                            op_arg = submul!(list, x, y, op_arg)
                            neg = false
                        end
                    else
                        y = add_instructions!(list, vᵢ_args[2], cse, pse)
                        if isnothing(op_arg)
                            op_arg = add_instructions!(list, vᵢ_args[2], cse, pse)
                            neg = true
                        elseif neg
                            op_arg = add!(list, y, op_arg)
                            neg = true
                        else
                            op_arg = sub!(list, op_arg, y)
                            neg = false
                        end
                    end

                elseif isnothing(op_arg)
                    op_arg = add_instructions!(list, vᵢ, cse, pse)
                else
                    x = add_instructions!(list, vᵢ_args[1], cse, pse)
                    y = _add_mul_instructions!(list, view(vᵢ_args, 2:nᵢ), cse, pse)
                    if neg
                        op_arg = mulsub!(list, x, y, op_arg)
                        neg = false
                    else
                        op_arg = muladd!(list, x, y, op_arg)
                    end
                end
            elseif isnothing(op_arg)
                op_arg = add_instructions!(list, vᵢ, cse, pse)
            elseif neg
                op_arg = sub!(list, add_instructions!(list, vᵢ, cse, pse), op_arg)
                neg = false
            else
                op_arg = add!(list, add_instructions!(list, vᵢ, cse, pse), op_arg)
            end
        end
        neg ? neg!(list, op_arg) : op_arg
    elseif t == :Div
        x, y = ModelKit.args(ex)
        return div!(
            list,
            add_instructions!(list, x, cse, pse),
            add_instructions!(list, y, cse, pse),
        )
    elseif t == :Sin
        x = add_instructions!(list, ModelKit.args(ex)[1], cse, pse)
        return push!(list, Instruction(INSTR_SIN, x))
    elseif t == :Cos
        x = add_instructions!(list, ModelKit.args(ex)[1], cse, pse)
        return push!(list, Instruction(INSTR_COS, x))
    else
        return to_number(ex)
    end
end

function _add_mul_instructions!(list, vec, cse, pse)
    length(vec) == 1 && return add_instructions!(list, vec[1], cse, pse)
    op_arg = nothing
    is_neg = false
    divs = []
    for vec_arg in vec
        if class(vec_arg) == :Pow
            x, k_expr = args(vec_arg)
            k = convert(Int32, k_expr)
            xinstr = add_instructions!(list, x, cse, pse)
            if k == -1
                push!(divs, xinstr)
            elseif k < -1
                push!(divs, pow!(list, xinstr, -k))
            else
                if isnothing(op_arg)
                    op_arg = pow!(list, xinstr, k)
                else
                    op_arg = mul!(list, op_arg, pow!(list, xinstr, k))
                end
            end
        elseif is_minus_one(vec_arg)
            is_neg = !is_neg
        else
            vinstr = add_instructions!(list, vec_arg, cse, pse)
            if isnothing(op_arg)
                op_arg = vinstr
            else
                op_arg = mul!(list, op_arg, vinstr)
            end
        end
    end

    if !isempty(divs)
        denom = first(divs)
        for i = 2:length(divs)
            denom = mul!(list, denom, divs[i])
        end
        if isnothing(op_arg)
            op_arg = div!(list, Int32(1), denom)
        else
            op_arg = div!(list, op_arg, denom)
        end
    end
    is_neg ? neg!(list, op_arg) : op_arg
end

struct DiffMap
    D::Dict{Tuple{Any,Int},Any}
end
DiffMap() = DiffMap(Dict{Tuple{Any,Int},Any}())
#
Base.getindex(D::DiffMap, var, ∂i) = get(D.D, (var, ∂i), nothing)
Base.setindex!(D::DiffMap, val, var, ∂i) = D.D[(var, ∂i)] = val
Base.setindex!(D::DiffMap, ::Nothing, var, ∂i) = nothing

"""
    gradient(list::InstructionList, N::Int, diff_map)

Returns an `InstructionList` list computing the gradients all instructions.
"""
function gradient(list::InstructionList, assignments, variables::Vector{Symbol})
    D = DiffMap()
    for (∂i, var) in enumerate(variables)
        D[var, ∂i] = Int32(1)
    end
    N = length(variables)

    grad_list = InstructionList()
    list_ref_mapping = Dict{InstructionRef,InstructionRef}()
    for (instr_ref, instr) in list
        op = instr.op
        arg₁, arg₂, arg₃ = get.(Ref(list_ref_mapping), instr.args, instr.args)

        # handle pow first here since we change the original instruction
        if op == INSTR_POW
            base = pow!(grad_list, arg₁, arg₂ - 1)
            new_ref = mul!(grad_list, arg₁, base)
            for i = 1:N
                D[new_ref, i] = mul!(grad_list, arg₂, mul!(grad_list, base, D[arg₁, i]))
            end
        else
            new_ref = push!(grad_list, Instruction(op, (arg₁, arg₂, arg₃)))
        end
        list_ref_mapping[instr_ref] = new_ref

        if op == INSTR_ADD
            for i = 1:N
                D[new_ref, i] = add!(grad_list, D[arg₁, i], D[arg₂, i])
            end
        elseif op == INSTR_SUB
            for i = 1:N
                D[new_ref, i] = sub!(grad_list, D[arg₁, i], D[arg₂, i])
            end
        elseif op == INSTR_NEG
            for i = 1:N
                D[new_ref, i] = neg!(grad_list, D[arg₁, i])
            end
        elseif op == INSTR_MUL
            for i = 1:N
                D[new_ref, i] =
                    muladd!(grad_list, arg₁, D[arg₂, i], mul!(grad_list, D[arg₁, i], arg₂))
            end
        elseif op == INSTR_DIV
            denom = sqr!(grad_list, arg₂)
            for i = 1:N
                D[new_ref, i] = div!(
                    grad_list,
                    mulsub!(grad_list, D[arg₁, i], arg₂, mul!(grad_list, arg₁, D[arg₂, i])),
                    denom,
                )
            end
        elseif op == INSTR_MULADD
            for i = 1:N
                D[new_ref, i] = muladd!(
                    grad_list,
                    arg₁,
                    D[arg₂, i],
                    muladd!(grad_list, D[arg₁, i], arg₂, D[arg₃, i]),
                )
            end
        elseif op == INSTR_MULSUB
            for i = 1:N
                D[new_ref, i] = muladd!(
                    grad_list,
                    arg₁,
                    D[arg₂, i],
                    mulsub!(grad_list, D[arg₁, i], arg₂, D[arg₃, i]),
                )
            end
        elseif op == INSTR_SUBMUL
            for i = 1:N
                D[new_ref, i] = submul!(
                    grad_list,
                    arg₁,
                    D[arg₂, i],
                    submul!(grad_list, D[arg₁, i], arg₂, D[arg₃, i]),
                )
            end
        elseif op == INSTR_SQR
            for i = 1:N
                D[new_ref, i] = mul!(grad_list, mul!(grad_list, 2, arg₁), D[arg₁, i])
            end
        elseif op == INSTR_SIN
            for i = 1:N
                if !isnothing(D[arg₁, i])
                    D[new_ref, i] = mul!(
                        grad_list,
                        push!(grad_list, Instruction(INSTR_COS, arg₁)),
                        D[arg₁, i],
                    )
                end
            end
        elseif op == INSTR_COS
            for i = 1:N
                if !isnothing(D[arg₁, i])
                    D[new_ref, i] = neg!(
                        grad_list,
                        mul!(
                            grad_list,
                            push!(grad_list, Instruction(INSTR_SIN, arg₁)),
                            D[arg₁, i],
                        ),
                    )
                end
            end
        elseif op == INSTR_POW
            # handled above
        else
            error(op)
            :UNKNOWN
        end
    end

    new_assignments = [list_ref_mapping[instr] for instr in assignments]
    grad_assignments = map([D[instr, ∂i] for instr in new_assignments, ∂i = 1:N]) do x
        isnothing(x) ? 0 : x
    end
    grad_list, new_assignments, grad_assignments
end


"""
    allocation_alias(list::InstructionList)

This performs an analysis of the recorded instructions to figure to decrease the number of
allocations necessary.

## Example

```julia
julia> ex = (a * b * (h - u) - c * d) * (h - u);

julia> list, _ = ModelKit.instruction_list([ex]);

julia> ModelKit.to_expr(list)
quote
    ι1 = sub_fast(h, u)
    ι2 = mul_fast(c, d)
    ι3 = mul_fast(b, ι1)
    ι4 = mulsub_fast(a, ι3, ι2)
    ι5 = mul_fast(ι1, ι4)
end

julia> ModelKit.to_expr(list; allocation_alias = ModelKit.allocation_alias(list))
quote
    ι1 = sub_fast(h, u)
    ι2 = mul_fast(c, d)
    ι3 = mul_fast(b, ι1)
    ι3 = mulsub_fast(a, ι3, ι2)
    ι3 = mul_fast(ι1, ι3)
end

```
"""
function allocation_alias(list::InstructionList, assignments::Set{Int})
    n = length(list.lifetime)
    lifetimes = map((i, x) -> i ∈ assignments ? typemax(Int) : x, 1:n, list.lifetime)
    alias = Int[]
    open_slots = Int[]
    max_slots = 0
    deaths = sort!([(k, i) for (i, k) in enumerate(lifetimes)])
    death_index = 1

    for (k, Δ) in enumerate(lifetimes)
        # check for deaths
        n₁ = length(open_slots)
        # we use here < k and not ≤ k to avoid the potential of accidental overwrites due
        # to shared memory
        while first(deaths[death_index]) < k
            push!(open_slots, alias[last(deaths[death_index])])
            death_index += 1
        end

        if isempty(open_slots)
            max_slots += 1
            push!(alias, max_slots)
        else
            push!(alias, pop!(open_slots))
        end

        n₂ = length(open_slots)
        if n₁ != n₂
            sort!(open_slots, alg = InsertionSort)
        end
    end

    alias
end

to_expr_arg(i::InstructionRef, allocation_alias) = Symbol("ι", allocation_alias[i.i])
to_expr_arg(i::InstructionRef, ::Nothing = nothing) = Symbol("ι", i.i)
to_expr_arg(x, _allocation_alias = nothing) = x


function to_expr(
    f::F,
    list::InstructionList;
    variables = Dict{Symbol,Union{Expr,Symbol}}(),
    assignements = Dict{Symbol,Vector{Expr}}(),
    allocation_alias = nothing,
) where {F}
    block = Expr(:block)
    exprs = block.args
    for (instr_ref, el) in list
        id = to_expr_arg(instr_ref, allocation_alias)
        op, (arg1, arg2, arg3) = el.op, el.args
        a = to_expr_arg(get(variables, arg1, arg1), allocation_alias)
        if arg3 !== nothing && arg2 !== nothing
            c = to_expr_arg(get(variables, arg3, arg3), allocation_alias)
            b = to_expr_arg(get(variables, arg2, arg2), allocation_alias)
            push!(exprs, :($id = $(f(op, (a, b, c)))))
        elseif arg2 !== nothing
            if op == INSTR_POW && arg2 isa Int
                k::Int = arg2
                if k < 0
                    push!(exprs, :($id = inv_fast($(unroll_pow(get(variables, a, a), -k)))))
                else
                    push!(exprs, :($id = $(unroll_pow(get(variables, a, a), k))))
                end
            else
                b = to_expr_arg(get(variables, arg2, arg2), allocation_alias)
                push!(exprs, :($id = $(f(op, (a, b)))))
            end
        else
            push!(exprs, :($id = $(f(op, (a,)))))
        end

        if haskey(assignements, id)
            append!(exprs, assignements[id])
        end
    end
    block
end
to_expr(list::InstructionList; kwargs...) = to_expr(list; kwargs...) do op, args
    Expr(:call, call_op(op), args...)
end

# HIGHER ORDER DIFFERENTIATION
# This follows Chapter 13 of Griewank & Walther - Evaluating derivatives

# Taylor
function untuple(varname, dx)
    :($(Expr(:tuple, (Symbol(varname, k) for k = 0:dx)...)) = $varname)
end
taylor_tuple(ids) = Expr(:tuple, [isnothing(id) ? :(zero(x0)) : id for id in ids]...)

function taylor_impl(f!, K::Int, dx::Int)
    D = DiffMap()
    list = InstructionList()
    for k = 0:dx
        D[:x, k] = Symbol(:x, k)
    end
    ids = to_expr_arg.(f!(list, D), nothing)
    quote
        $(untuple(:x, dx))
        $(to_expr(list; variables = Dict(:x => :x)))
        $(taylor_tuple(ids))
    end
end
function taylor_impl(f!, K::Int, dx::Int, dy::Int)
    D = DiffMap()
    list = InstructionList()
    for k = 0:dx
        D[:x, k] = Symbol(:x, k)
    end
    for k = 0:dy
        D[:y, k] = Symbol(:y, k)
    end
    ids = to_expr_arg.(f!(list, D), nothing)
    quote
        $(untuple(:x, dx))
        $(untuple(:y, dy))
        $(to_expr(list; variables = Dict(:x => :x, :y => :y)))
        $(taylor_tuple(ids))
    end
end
function taylor_impl(f!, K::Int, dx::Int, dy::Int, dz::Int)
    D = DiffMap()
    list = InstructionList()
    for k = 0:dx
        D[:x, k] = Symbol(:x, k)
    end
    for k = 0:dy
        D[:y, k] = Symbol(:y, k)
    end
    for k = 0:dz
        D[:z, k] = Symbol(:z, k)
    end
    ids = to_expr_arg.(f!(list, D), nothing)
    quote
        $(untuple(:x, dx))
        $(untuple(:y, dy))
        $(untuple(:z, dz))
        $(to_expr(list; variables = Dict(:x => :x, :y => :y, :z => :z)))
        $(taylor_tuple(ids))
    end
end

make_ntuple(t::T) where {T>:Tuple} = (t,)
make_ntuple(t::Tuple{A}) where {A} = t
make_ntuple(t::Tuple{A,B}) where {A,B} = convert(NTuple{2,promote_type(A, B)}, t)
make_ntuple(t::Tuple{A,B,C}) where {A,B,C} = convert(NTuple{3,promote_type(A, B, C)}, t)
make_ntuple(t::Tuple{A,B,C,D}) where {A,B,C,D} =
    convert(NTuple{4,promote_type(A, B, C, D)}, t)
make_ntuple(t::Tuple{A,B,C,D,E}) where {A,B,C,D,E} =
    convert(NTuple{5,promote_type(A, B, C, D, E)}, t)


for f in [
    :taylor_add,
    :taylor_sub,
    :tayor_neg,
    :taylor_mul,
    :taylor_div,
    :taylor_muladd,
    :taylor_mulsub,
    :taylor_submul,
    :taylor_sqr,
    :taylor_pow,
    :taylor_sin,
    :taylor_cos,
]
    @eval begin
        $f(ord::O, x::X) where {V<:Val,O<:Val,X} = $f(ord, make_ntuple(x))
        $f(ord::O, x, y) where {V<:Val,O<:Val,X,Y} = $f(ord, make_ntuple(x), make_ntuple(y))
        $f(ord::O, x::X, y::Y, z::Z) where {O<:Val,X,Y,Z} =
            $f(ord, make_ntuple(x), make_ntuple(y), make_ntuple(z))
    end
end

@generated function taylor_neg(::Val{K}, x::NTuple{M}) where {K,M}
    taylor_impl(K, M - 1) do list, D
        [neg!(list, D[:x, k]) for k = 0:K]
    end
end

function taylor_sqr_impl(K::Int, dx::Int)
    taylor_impl(K, dx) do list, D
        map(0:K) do k
            # Compute ∑_{j=0}^k u_j U_{k - j} =
            #         2(∑_{j=0}^div(k - 1,2) u_j U_{k - j}) +
            #         iseven(k) * u_{k-1}^2
            k == 0 && return sqr!(list, D[:x, 0])
            w_k = nothing
            for j = 0:div(k - 1, 2)
                w_k = muladd!(list, D[:x, j], D[:x, k-j], w_k)
            end
            if iseven(k)
                muladd!(list, 2, w_k, sqr!(list, D[:x, div(k, 2)]))
            else
                add!(list, w_k, w_k)
            end
        end
    end
end
@generated function taylor_sqr(::Val{K}, x::NTuple{M}) where {K,M}
    taylor_sqr_impl(K, M - 1)
end

@generated function taylor_add(::Val{K},
    x::NTuple{M},
    y::NTuple{N},
) where {K,M,N}
    taylor_impl(K, M - 1, N - 1) do list, D
        [add!(list, D[:x, k], D[:y, k]) for k = 0:K]
    end
end

@generated function taylor_sub(
    ::Val{K},
    x::NTuple{M},
    y::NTuple{N},
) where {K,M,N}
    taylor_impl(K, M - 1, N - 1) do list, D
        [sub!(list, D[:x, k], D[:y, k]) for k = 0:K]
    end
end

@generated function taylor_mul(
    ::Val{K},
    x::NTuple{M},
    y::NTuple{N},
) where {K,M,N}
    taylor_impl(K, M - 1, N - 1) do list, D
        map(0:K) do k
            c_k = nothing
            for j = 0:k
                c_k = muladd!(list, D[:x, j], D[:y, k-j], c_k)
            end
            c_k
        end
    end
end

@generated function taylor_div(
    ::Val{K},
    x::NTuple{M},
    y::NTuple{N},
) where {K,M,N}
    taylor_impl(K, M - 1, N - 1) do list, D
        ids = []
        for k = 0:K
            s = nothing
            for j = 0:(k-1)
                s = muladd!(list, ids[j+1], D[:y, k-j], s)
            end
            push!(ids, div!(list, sub!(list, D[:x, k], s), D[:y, 0]))
        end
        ids
    end
end

function taylor_pow_impl(K, dx)
    D = DiffMap()
    list = InstructionList()
    for k = 0:dx
        D[:x, k] = Symbol(:x, k)
    end
    w = Any[:w₀]
    for k = 1:K
        s = nothing
        for j = 1:k
            ũ_j = mul!(list, j, D[:x, j])
            w_kj = k == j ? :w₀ : w[k-j+1]
            s = muladd!(list, w_kj, ũ_j, s)
        end
        s = mul!(list, :r, s)

        t = nothing
        for j = 1:k-1
            u_kj = D[:x, k-j]
            w̃_j = mul!(list, j, w[j+1])
            t = muladd!(list, u_kj, w̃_j, t)
        end
        w̃_k = mul!(list, :u₀_inv, sub!(list, s, t))
        w_k = div!(list, w̃_k, k)
        push!(w, w_k)
    end
    quote
        $(untuple(:x, dx))
        iszero(x0) && return $(taylor_tuple([nothing for _ = 0:K]))
        w₀ = pow_fast(x0, r)
        $((K > 0 ? (:(u₀_inv = inv_fast(x0)),) : ())...)
        $(to_expr(list; variables = Dict(:x => :x, :u₀_inv => :u₀_inv, :w₀ => :w₀)))
        $(taylor_tuple(to_expr_arg.(w, nothing)))
    end
end
@generated function taylor_pow(::Val{K}, x::NTuple{M}, r::Integer) where {K,M}
    taylor_pow_impl(K, M - 1)
end

@generated function taylor_muladd(
    ::Val{K},
    x::NTuple{M},
    y::NTuple{N},
    z::NTuple{L},
) where {K,M,N,L}
    taylor_impl(K, M - 1, N - 1, L - 1) do list, D
        map(0:K) do k
            c_k = D[:z, k]
            for j = 0:k
                c_k = muladd!(list, D[:x, j], D[:y, k-j], c_k)
            end
            c_k
        end
    end
end
@generated function taylor_mulsub(
    ::Val{K},
    x::NTuple{M},
    y::NTuple{N},
    z::NTuple{L},
) where {K,M,N,L}
    taylor_impl(K, M - 1, N - 1, L - 1) do list, D
        map(0:K) do k
            c_k = D[:z, k]
            neg = false
            for j = 0:k
                if neg
                    c_k = muladd!(list, D[:x, j], D[:y, k-j], c_k)
                else
                    c_k = mulsub!(list, D[:x, j], D[:y, k-j], c_k)
                    neg = true
                end
            end
            neg ? c_k : neg!(list, c_k)
        end
    end
end

@generated function taylor_submul(
    ::Val{K},
    x::NTuple{M},
    y::NTuple{N},
    z::NTuple{L},
) where {K,M,N,L}
    taylor_impl(K, M - 1, N - 1, L - 1) do list, D
        map(0:K) do k
            c_k = D[:z, k]
            for j = 0:k
                c_k = submul!(list, D[:x, j], D[:y, k-j], c_k)
            end
            c_k
        end
    end
end

function taylor_sin_helper!(list, D, k)
    k == 0 && return :s₀
    s_k = nothing
    for j = 1:k
        s_k =
            muladd!(list, mul!(list, j, D[:x, j]), taylor_cos_helper!(list, D, k - j), s_k)
    end
    div!(list, s_k, k)
end
function taylor_cos_helper!(list, D, k)
    k == 0 && return :c₀
    c_k = nothing
    for j = 1:k
        c_k =
            submul!(list, mul!(list, j, D[:x, j]), taylor_sin_helper!(list, D, k - j), c_k)
    end
    div!(list, c_k, k)
end
function taylor_sin_impl(K, M)
    D = DiffMap()
    list = InstructionList()
    for k = 0:M-1
        D[:x, k] = Symbol(:x, k)
    end

    ids = Any[:s₀]
    for k = 1:K
        push!(ids, taylor_sin_helper!(list, D, k))
    end

    quote
        $(untuple(:x, M - 1))
        $(K > 0 ? :((s₀, c₀) = sincos(x0)) : :(s₀ = sin(x0)))
        $(to_expr(list; variables = Dict(:x => :x, :s₀ => :s₀, :c₀ => :c₀)))
        $(taylor_tuple(to_expr_arg.(ids, nothing)))
    end
end
function taylor_cos_impl(K, M)
    D = DiffMap()
    list = InstructionList()
    for k = 0:M-1
        D[:x, k] = Symbol(:x, k)
    end

    ids = Any[:c₀]
    for k = 1:K
        push!(ids, taylor_cos_helper!(list, D, k))
    end

    quote
        $(untuple(:x, M - 1))
        $(K > 0 ? :((s₀, c₀) = sincos(x0)) : :(c₀ = cos(x0)))
        $(to_expr(list; variables = Dict(:x => :x, :s₀ => :s₀, :c₀ => :c₀)))
        $(taylor_tuple(to_expr_arg.(ids, nothing)))
    end
end
@generated function taylor_sin(::Val{K}, x::NTuple{M}) where {K,M}
    taylor_sin_impl(K, M)
end
@generated function taylor_cos(::Val{K}, x::NTuple{M}) where {K,M}
    taylor_cos_impl(K, M)
end
