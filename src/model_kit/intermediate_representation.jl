struct IRStatementRef
    i::Int
end

const IRStatementArg = Union{Nothing, Number, Symbol, IRStatementRef}

struct IRStatement
    op::OpType
    target::IRStatementRef
    args::NTuple{4, IRStatementArg}
end

struct IntermediateRepresentation
    statements::Vector{IRStatement}
    # tape_index, assignment_ref
    assignments::Vector{Tuple{Int, IRStatementArg}}
    output_dim::Int
end
IntermediateRepresentation() = IntermediateRepresentation([], [], 0)

## IRStatementRef
Base.show(io::IO, i::IRStatementRef) = print(io, "ι_", i.i)

# ## IRStatement
function IRStatement(op::OpType, target::Union{Symbol, IRStatementRef}, a::IRStatementArg)
    arity(op) == 1 || error("IRStatement has arity $(arity(op)) but called with 1 arg")
    return IRStatement(op, target, (a, nothing, nothing, nothing))
end
function IRStatement(
        op::OpType,
        target::Union{Symbol, IRStatementRef},
        a::IRStatementArg,
        b::IRStatementArg,
    )
    arity(op) == 2 || error("IRStatement has arity $(arity(op)) but called with 2 args")
    return IRStatement(op, target, (a, b, nothing, nothing))
end
function IRStatement(
        op::OpType,
        target::Union{Symbol, IRStatementRef},
        a::IRStatementArg,
        b::IRStatementArg,
        c::IRStatementArg,
    )
    arity(op) == 3 || error("IRStatement has arity $(arity(op)) but called with 3 args")
    return IRStatement(op, target, (a, b, c, nothing))
end
function IRStatement(
        op::OpType,
        target::Union{Symbol, IRStatementRef},
        a::IRStatementArg,
        b::IRStatementArg,
        c::IRStatementArg,
        d::IRStatementArg,
    )
    arity(op) == 4 || error("IRStatement has arity $(arity(op)) but called with 4 args")
    return IRStatement(op, target, (a, b, c, d))
end
Base.getindex(I::IRStatement, i) = I.args[i]


##  IntermediateRepresentation

IntermediateRepresentation(output_dim) =
    IntermediateRepresentation(IRStatement[], IRStatementRef[], output_dim)
Base.length(v::IntermediateRepresentation) = length(v.statements)
function Base.iterate(I::IntermediateRepresentation, state = nothing)
    return isnothing(state) ? iterate(I.statements) : iterate(I.statements, state)
end

function Base.show(io::IO, ir::IntermediateRepresentation)
    for stmt in ir
        println(
            io,
            "$(stmt.target) <- $(
                Expr(
                    :call,
                    Symbol(stmt.op),
                    filter(!isnothing, stmt.args)...,
                )
            )",
        )
    end
    return
end

function Base.push!(ir::IntermediateRepresentation, x::IRStatement)
    push!(ir.statements, x)
    return IRStatementRef(length(ir))
end

function add_op!(ir::IntermediateRepresentation, op::OpType, args::IRStatementArg...)
    stmt = n = length(ir) + 1
    ref = IRStatementRef(n)
    stmt = IRStatement(op, ref, args...)
    push!(ir.statements, stmt)

    return IRStatementRef(n)
end


## Convert to IntermediateRepresentation

intermediate_representation(ex::Expression; kwargs...) =
    intermediate_representation([ex]; kwargs...)
function intermediate_representation(
        exprs::Vector{Expression};
        output_dim = length(exprs),
        subexpressions::Dict{Expression, Expression} = Dict{Expression, Expression}(),
        perform_cse::Bool = true,
    )
    ir = IntermediateRepresentation(output_dim)
    if perform_cse
        exprs, CSE = cse(exprs)
        merge!(CSE, subexpressions)
    else
        CSE = subexpressions
    end
    symb_cse = Dict{Symbol, Expression}()
    for (k, v) in CSE
        push!(symb_cse, Symbol(to_string(k)) => v)
    end
    pse = Dict{Symbol, IRStatementRef}()
    for (k, ex) in enumerate(exprs)
        iszero(ex) && continue
        n = length(ir.statements)
        assignment = expr_to_ir_statements!(ir, ex, symb_cse, pse)
        # Add identity if result is a no-op to handle constants
        if (length(ir.statements) == n)
            ref = add_op!(ir, OP_IDENTITY, assignment)
            push!(ir.assignments, (k, ref))
        else
            push!(ir.assignments, (k, assignment))
        end
    end
    return ir
end


# Helpers
function add!(ir::IntermediateRepresentation, @nospecialize(a), @nospecialize(b))
    return if a === nothing
        b === nothing ? nothing : b
    else
        b === nothing ? a : add_op!(ir, OP_ADD, a, b)
    end
end
function neg!(ir::IntermediateRepresentation, @nospecialize(a))
    return a === nothing ? nothing : add_op!(ir, OP_NEG, a)
end
function sub!(ir::IntermediateRepresentation, @nospecialize(a), @nospecialize(b))
    return if a === nothing
        b === nothing ? nothing : add_op!(ir, OP_NEG, b)
    else
        b === nothing ? a : add_op!(ir, OP_SUB, a, b)
    end
end

function mul!(ir::IntermediateRepresentation, @nospecialize(a), @nospecialize(b))
    (a === nothing || b === nothing) && return nothing
    return if is_one(a)
        b
    elseif is_one(b)
        a
    elseif is_minus_one(a)
        add_op!(ir, OP_NEG, b)
    elseif is_minus_one(b)
        add_op!(ir, OP_NEG, a)
    elseif a isa Integer && a == 2
        add_op!(ir, OP_ADD, b, b)
    else
        add_op!(ir, OP_MUL, a, b)
    end
end

function muladd!(
        ir::IntermediateRepresentation,
        @nospecialize(a),
        @nospecialize(b),
        @nospecialize(c),
    )
    (a === nothing || b === nothing) && return c
    c === nothing && return mul!(ir, a, b)
    is_one(a) && return add!(ir, b, c)
    is_one(b) && return add!(ir, a, c)
    return add_op!(ir, OP_MULADD, a, b, c)
end
function mulsub!(
        ir::IntermediateRepresentation,
        @nospecialize(a),
        @nospecialize(b),
        @nospecialize(c),
    )
    (a === nothing || b === nothing) && return neg!(ir, c)
    c === nothing && return mul!(ir, a, b)
    is_one(a) && return sub!(ir, b, c)
    is_one(b) && return sub!(ir, a, c)
    return add_op!(ir, OP_MULSUB, a, b, c)
end
function submul!(
        ir::IntermediateRepresentation,
        @nospecialize(a),
        @nospecialize(b),
        @nospecialize(c),
    )
    # c - a * b
    (a === nothing || b === nothing) && return c
    c === nothing && return neg!(ir, mul!(ir, a, b))
    is_one(a) && return sub!(ir, c, b)
    is_one(b) && return sub!(ir, c, a)
    return add_op!(ir, OP_SUBMUL, a, b, c)
end

function mulmuladd!(ir::IntermediateRepresentation, a, b, c, d)
    return add_op!(ir, OP_MULMULADD, a, b, c, d)
end
function div!(ir::IntermediateRepresentation, @nospecialize(a), @nospecialize(b))
    a === nothing && return nothing
    return is_one(b) ? a : add_op!(ir, OP_DIV, a, b)
end

function sqr!(ir::IntermediateRepresentation, @nospecialize(a))
    a === nothing && return nothing
    return add_op!(ir, OP_SQR, a)
end
function pow!(ir::IntermediateRepresentation, @nospecialize(a), @nospecialize(b))
    is_zero(b) && return Int32(1)
    is_one(b) && return a
    k = to_number(b)
    if k isa Integer
        k == 2 && return sqr!(ir, a)
        k == 3 && return add_op!(ir, OP_CB, a)
        k == -1 && return add_op!(ir, OP_INV, a)
        k == -2 && return add_op!(ir, OP_INVSQR, a)
        return add_op!(ir, OP_POW_INT, a, k)
    elseif k == 1 // 2
        return add_op!(ir, OP_SQRT, a)
    else
        error("Cannot handle exponent: $b")
    end
end

function split_off_minus_one(ex)
    if class(ex) === :Mul
        v = args(ex)
        if is_minus_one(v[1])
            return -1, prod(@view v[2:end])
        end
        return 1, ex
    end
    return 1, ex
end

function split_into_positives_negatives(ex)
    positives = Basic[]
    negatives = Basic[]
    v = args(ex)
    for e in v
        (sign, val) = split_off_minus_one(e)
        if sign == -1
            push!(negatives, val)
        else
            push!(positives, val)
        end
    end
    return positives, negatives
end

function reduce_to_at_most_two_multiplicants!(ir, ex, cse, pse)
    if class(ex) === :Mul
        v = args(ex)
        if length(v) == 2
            return (
                expr_to_ir_statements!(ir, v[1], cse, pse),
                expr_to_ir_statements!(ir, v[2], cse, pse),
            )
        elseif length(v) == 1
            return (expr_to_ir_statements!(ir, v[1], cse, pse), nothing)
        elseif length(v) > 2
            v2 = expr_to_ir_statements!(ir, prod(@view v[1:(end - 1)]), cse, pse)
            return (expr_to_ir_statements!(ir, v[end], cse, pse), v2)
        end
    end
    return (expr_to_ir_statements!(ir, ex, cse, pse), nothing)
end

function sum_products!(ir, tuples::Vector)
    if isempty(tuples)
        return nothing
    end
    singles = convert(Vector{IRStatementArg}, first.(filter(t -> isnothing(t[2]), tuples)))
    pairs = filter(t -> !isnothing(t[2]), tuples)
    n = length(pairs)
    for k in 1:2:(n - 1)
        (a, b) = pairs[k]
        (c, d) = pairs[k + 1]
        ref = mulmuladd!(ir, a, b, c, d)
        push!(singles, ref)
    end


    if isodd(n)
        # we just had a single product
        if isempty(singles)
            return mul!(ir, pairs[n][1], pairs[n][2])
        end
        (a, b) = pairs[n]
        c = pop!(singles)
        ref = muladd!(ir, a, b, c)
        push!(singles, ref)
    end

    # Now add up singles
    while length(singles) > 1
        if length(singles) >= 4
            a = pop!(singles)
            b = pop!(singles)
            c = pop!(singles)
            d = pop!(singles)
            ref = add_op!(ir, OP_ADD4, a, b, c, d)
            push!(singles, ref)
        elseif length(singles) == 3
            a = pop!(singles)
            b = pop!(singles)
            c = pop!(singles)
            ref = add_op!(ir, OP_ADD3, a, b, c)
            push!(singles, ref)
        elseif length(singles) == 2
            a = pop!(singles)
            b = pop!(singles)
            ref = add_op!(ir, OP_ADD, a, b)
            push!(singles, ref)
        end
    end
    return singles[1]
end


function process_sum!(ir, ex::Basic, cse, pse)
    @assert class(ex) == :Add
    pos, neg = split_into_positives_negatives(ex)
    pos_reduced = map(e -> reduce_to_at_most_two_multiplicants!(ir, e, cse, pse), pos)
    neg_reduced = map(e -> reduce_to_at_most_two_multiplicants!(ir, e, cse, pse), neg)


    if length(pos_reduced) == 1 && length(neg_reduced) == 1
        (a, b) = pos_reduced[1]
        (c, d) = neg_reduced[1]
        if !isnothing(b) && !isnothing(d)
            return add_op!(ir, OP_MULMULSUB, a, b, c, d)
        elseif !isnothing(b)
            return add_op!(ir, OP_MULSUB, a, b, c)
        elseif !isnothing(d)
            # a - c * d
            return add_op!(ir, OP_SUBMUL, c, d, a)
        else
            return add_op!(ir, OP_SUB, a, c)
        end
    elseif length(neg_reduced) == 1
        a = sum_products!(ir, pos_reduced)
        (c, d) = neg_reduced[1]
        if !isnothing(d)
            # a - c * d
            return add_op!(ir, OP_SUBMUL, c, d, a)
        else
            return add_op!(ir, OP_SUB, a, c)
        end
    end

    pos_sum = sum_products!(ir, pos_reduced)
    neg_sum = sum_products!(ir, neg_reduced)
    return sub!(ir, pos_sum, neg_sum)
end


function split_into_num_denom!(ir, ex, cse, pse)
    nums = []
    denoms = []
    v = args(ex)
    for e in v
        if class(e) == :Pow
            x, k_ex = args(e)
            xref = expr_to_ir_statements!(ir, x, cse, pse)
            k = to_number(k_ex)
            if k isa Basic
                throw(ExprError("Cannot handle non-constant exponents"))
            end
            if k < 0
                push!(denoms, pow!(ir, xref, -k))
            else
                push!(nums, pow!(ir, xref, k))
            end
        else
            push!(nums, expr_to_ir_statements!(ir, e, cse, pse))
        end
    end
    return nums, denoms
end

function prod_parts!(ir, exs, cse, pse)
    isempty(exs) && return nothing
    # We multiply from the reverse to
    # possibly detect 2 * x (since constants are always in the front)
    parts = reverse(Any[(e) for e in exs])
    # Now add up parts
    while length(parts) > 1
        if length(parts) >= 4
            a = pop!(parts)
            b = pop!(parts)
            c = pop!(parts)
            d = pop!(parts)
            ref = add_op!(ir, OP_MUL4, a, b, c, d)
            push!(parts, ref)
        elseif length(parts) == 3
            a = pop!(parts)
            b = pop!(parts)
            c = pop!(parts)
            ref = add_op!(ir, OP_MUL3, a, b, c)
            push!(parts, ref)
        elseif length(parts) == 2
            a = pop!(parts)
            b = pop!(parts)
            ref = mul!(ir, a, b)
            push!(parts, ref)
        end
    end
    return parts[1]
end

function process_mul!(ir, ex, cse, pse)
    @assert class(ex) == :Mul
    (m, ex) = split_off_minus_one(ex)
    if class(ex) != :Mul
        return mul!(ir, m, expr_to_ir_statements!(ir, ex, cse, pse))
    end
    num, denom = split_into_num_denom!(ir, ex, cse, pse)
    num_prod = prod_parts!(ir, num, cse, pse)
    denom_prod = prod_parts!(ir, denom, cse, pse)
    if isnothing(num_prod)
        ref = add_op!(ir, OP_INV, denom_prod)
    elseif isnothing(denom_prod)
        ref = num_prod
    else
        ref = div!(ir, num_prod, denom_prod)
    end
    return mul!(ir, m, ref)
end


function expr_to_ir_statements!(
        ir::IntermediateRepresentation,
        ex::Basic,
        # common sub expression
        cse::Dict{Symbol, Expression},
        # processed sub expression
        # mapping of the cse var with the instruction used
        pse::Dict{Symbol, IRStatementRef},
    )
    t = class(ex)
    if t == :Symbol
        s = Symbol(to_string(ex))
        if haskey(pse, s)
            return pse[s]
        elseif haskey(cse, s)
            val = expr_to_ir_statements!(ir, cse[s], cse, pse)
            if val isa IRStatementRef
                push!(pse, s => val)
            end
            val
        else
            s
        end
    elseif t == :Pow
        x, k = args(ex)
        x_ref = expr_to_ir_statements!(ir, x, cse, pse)
        k_ref = expr_to_ir_statements!(ir, k, cse, pse)
        pow!(ir, x_ref, k_ref)
    elseif t == :Mul
        process_mul!(ir, ex, cse, pse)
    elseif t == :Add
        process_sum!(ir, ex, cse, pse)
    elseif t == :Div
        x, y = args(ex)
        return div!(
            ir,
            expr_to_ir_statements!(ir, x, cse, pse),
            expr_to_ir_statements!(ir, y, cse, pse),
        )
    elseif t == :Sin
        x = expr_to_ir_statements!(ir, args(ex)[1], cse, pse)
        return add_op!(ir, OP_SIN, x)
    elseif t == :Cos
        x = expr_to_ir_statements!(ir, args(ex)[1], cse, pse)
        return add_op!(ir, OP_COS, x)
    else
        n = to_number(ex)
        if (n === ex)
            error("Cannot match expression to an instruction.")
        end
        return n
    end
end


to_expr_arg(i::IRStatementRef) = Symbol("ι", i.i)
to_expr_arg(x) = x

function to_julia_expr(ir::IntermediateRepresentation)
    exprs = map(ir.statements) do stmt
        :(
            $(to_expr_arg(stmt.target)) = $(
                Expr(
                    :call,
                    op_call(stmt.op),
                    to_expr_arg.(filter(!isnothing, stmt.args))...,
                )
            )
        )
    end
    return quote
        $(exprs...)
    end
end
