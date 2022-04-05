struct Instruction
    input::NTuple{4,Int32}
    op::OpType
    output::Int32
end

struct InstructionSequence
    instructions::Vector{Instruction}
    constants::Vector{Number}
    constants_range::UnitRange{Int}
    parameters_range::UnitRange{Int}
    variables_range::UnitRange{Int}
    continuation_parameter_index::Union{Nothing,Int}
    # tape_index, assignment_index
    assignments::Vector{Tuple{Int,Int}}
    output_dim::Int
    tape_space_needed::Int
end

function sequence_to_expr(
    seq::InstructionSequence;
    op_call = op_call,
    instr_symbol = :τ,
    unique_variable_names = true,
    order = nothing,
)
    tape_index_symbol_map = Dict{Int32,Symbol}()

    r = (k) -> begin
        if k in seq.constants_range
            :($(seq.constants[k]))
        elseif k in seq.variables_range
            :(x[$(k - first(seq.variables_range) + 1)])
        elseif k in seq.parameters_range
            :(p[$(k - first(seq.parameters_range) + 1)])
        elseif k == seq.continuation_parameter_index
            :t
        else
            tape_index_symbol_map[k]
        end
    end

    order_arg = isnothing(order) ? () : (order,)

    exprs = map(seq.instructions) do instr
        op = instr.op
        arg₁, arg₂, arg₃, arg₄ = instr.input
        i = instr.output

        a = arity(instr.op)

        # Make sure that we have unique variable names (if desired)
        # We do this by remapping an index to a new variable name (the current one for this index appended by an !)
        if op != OP_STOP
            if unique_variable_names && haskey(tape_index_symbol_map, i)
                tape_index_symbol_map[i] = Symbol(tape_index_symbol_map[i], :!)
            else
                tape_index_symbol_map[i] = Symbol(instr_symbol, i)
            end
        end

        if instr.op == OP_POW_INT
            :($(r(i)) = $(op_call(OP_POW_INT))($(order_arg...), $(r(arg₁)), $(arg₂)))
        elseif a == 1
            :($(r(i)) = $(op_call(op))($(order_arg...), $(r(arg₁))))
        elseif a == 2
            :($(r(i)) = $(op_call(op))($(order_arg...), $(r(arg₁)), $(r(arg₂))))
        elseif a == 3
            :($(r(i)) = $(op_call(op))($(order_arg...), $(r(arg₁)), $(r(arg₂)), $(r(arg₃))))
        elseif a == 4
            :(
                $(r(i)) = $(op_call(op))(
                    $(order_arg...),
                    $(r(arg₁)),
                    $(r(arg₂)),
                    $(r(arg₃)),
                    $(r(arg₄)),
                )
            )
        end
    end
    quote
        $(exprs...)
    end, r
end

function Base.show(io::IO, instr::Instruction)
    print(io, "t[$(instr.output)] = ", string(Symbol(instr.op)), "(")
    args = map(enumerate(instr.input[1:arity(instr.op)])) do (k, arg)
        if (should_use_index_not_reference(instr.op, k))
            "$arg"
        else
            "t[$arg]"
        end
    end

    join(io, args, ",")
    print(io, ")")
end

function Base.show(io::IO, seq::InstructionSequence)
    print(io, first(sequence_to_expr(seq; unique_variable_names = false)))
end


function should_use_index_not_reference(instr::Instruction, index)
    should_use_index_not_reference(instr.op, index)
end


Base.length(seq::InstructionSequence) = length(seq.instructions)
Base.iterate(seq::InstructionSequence) = iterate(seq.instructions)
Base.iterate(seq::InstructionSequence, state) = iterate(seq.instructions, state)

number_or_nothing(a::Number) = a
number_or_nothing(a) = nothing
getconstants(i::IRStatement) = (
    number_or_nothing(i.args[1]),
    number_or_nothing(i.args[2]),
    number_or_nothing(i.args[3]),
    number_or_nothing(i.args[4]),
)

function instruction_sequence(F::Union{System,Homotopy}; output_dim::Integer = length(F))
    vars = Symbol.(variables(F))
    params = Symbol.(parameters(F))
    continuation_parameter = F isa Homotopy ? Symbol(F.t) : nothing
    instruction_sequence(
        intermediate_representation(expressions(F); output_dim = output_dim);
        variables = vars,
        parameters = params,
        continuation_parameter = continuation_parameter,
    )
end

function instruction_sequence(
    ir::IntermediateRepresentation;
    variables::Vector{Symbol},
    parameters::Vector{Symbol} = Symbol[],
    continuation_parameter::Union{Symbol,Nothing} = nothing,
)
    arg_index_map = Dict{IRStatementArg,Int32}()
    tape_index = 0
    # Find all constants
    constants = Number[]
    for stmt in ir
        for (k, arg) in enumerate(getconstants(stmt))
            isnothing(arg) && continue
            should_use_index_not_reference(stmt.op, k) && continue
            haskey(arg_index_map, arg) && continue

            push!(constants, arg)
            arg_index_map[arg] = (tape_index += 1)

        end
    end
    nconstants = length(constants)
    constants_range = 1:nconstants

    # Build up indices
    nparams = length(parameters)
    for v in parameters
        arg_index_map[v] = (tape_index += 1)
    end
    parameters_range = range(nconstants + 1, length = nparams)

    continuation_parameter_index = nothing
    if !isnothing(continuation_parameter)
        arg_index_map[continuation_parameter] = (tape_index += 1)
        continuation_parameter_index = tape_index
    end

    nvars = length(variables)
    for v in variables
        arg_index_map[v] = (tape_index += 1)
    end
    variables_range = range(tape_index - nvars + 1, length = nvars)

    input_block_size = tape_index
    initial_space = length(ir)

    for (k, (_, ref)) in enumerate(ir.assignments)
        arg_index_map[ref] = input_block_size + initial_space + k
    end

    instructions = map(ir) do stmt
        prev_stmt_arg = Int32(1)
        instr_op = stmt.op
        args = ntuple(Val(4)) do i
            arg = stmt.args[i]
            stmt_arg = if isnothing(arg)
                # dummy data
                Int32(prev_stmt_arg)
            elseif should_use_index_not_reference(instr_op, i)
                Int32(arg)
            elseif arg isa IRStatementRef
                get(arg_index_map, arg, Int32(arg.i + input_block_size))
            else
                arg_index_map[arg]
            end
            prev_stmt_arg = stmt_arg
            return prev_stmt_arg
        end
        target_index = get(arg_index_map, stmt.target, stmt.target.i + input_block_size)

        Instruction(args, instr_op, target_index)
    end
    # optimize sequence
    assignments_range =
        range(input_block_size + initial_space + 1; length = length(ir.assignments))

    instructions, tape_space_needed, updated_assignments_range =
        optimize(instructions, input_block_size, assignments_range)

    updated_assignments = map(ir.assignments, updated_assignments_range) do (k, _), idx
        (k, idx)
    end

    # Add stop instruction
    n = tape_space_needed
    push!(instructions, Instruction((n, n, n, n), OP_STOP, n))

    InstructionSequence(
        instructions,
        constants,
        constants_range,
        parameters_range,
        variables_range,
        continuation_parameter_index,
        updated_assignments,
        ir.output_dim,
        tape_space_needed,
    )
end

function jacobian_instruction_sequence(F::System, ; kwargs...)
    JF = System(
        [F.expressions; vec(jacobian(F))],
        variables = variables(F),
        parameters = parameters(F),
    )
    instruction_sequence(JF; output_dim = length(F))
end
function jacobian_instruction_sequence(H::Homotopy, ; kwargs...)
    JH = Homotopy(
        [H.expressions; vec(differentiate(H.expressions, H.variables))],
        variables(H),
        H.t,
        parameters(H),
    )
    instruction_sequence(JH; output_dim = length(H))
end


"""
    optimize(instructions::Vector{Instruction}, input_block_size::Int, assignment_indices)

Optimize the instruction order and needed space.
"""
function optimize(instructions::Vector{Instruction}, input_block_size::Int, assignments)
    instructions = optimize_instruction_order(instructions)
    instructions, space_needed, updated_assignments =
        reduce_space(instructions, input_block_size, assignments)
    instructions, space_needed, updated_assignments
end


# Optimize instruction order
function optimize_instruction_order(instructions::Vector{Instruction})
    # Build a directed acyclic graph (dag) from the instruction input - output
    # depencies. The graph is directed from the assignments to the inputs
    g, index_to_instr = dag_from_instructions(instructions; start_from_inputs = false)


    listing = Int32[]
    # Roots are all assignments that are not reused
    roots = filter(v -> SimpleGraphs.in_deg(g, v) == 0, SimpleGraphs.vlist(g))

    while !isempty(roots)
        root = pop!(roots)
        unlisted_nodes_without_parent = [root]
        while !isempty(unlisted_nodes_without_parent)
            u = pop!(unlisted_nodes_without_parent)
            push!(listing, u)
            children = SimpleGraphs.out_neighbors(g, u)
            delete!(g, u)
            for v in children
                if SimpleGraphs.in_deg(g, v) == 0 && SimpleGraphs.out_deg(g, v) > 0
                    push!(listing, v)
                    children = SimpleGraphs.out_neighbors(g, v)
                    delete!(g, v)
                end
            end

            if isempty(unlisted_nodes_without_parent)
                is_unlisted =
                    v ->
                        SimpleGraphs.in_deg(g, v) == 0 &&
                            SimpleGraphs.out_deg(g, v) > 0 &&
                            v ∉ roots
                unlisted_nodes_without_parent = filter(is_unlisted, SimpleGraphs.vlist(g))
            end
        end
    end
    reverse!(listing)

    return getindex.(Ref(index_to_instr), listing)
end

"""
Build a directed graph from the instructions.
"""
function dag_from_instructions(instructions::Vector{Instruction}; start_from_inputs = false)
    index_to_instr = Dict{Int32,Instruction}()
    g = SimpleGraphs.SimpleDigraph{Int32}()
    for instr in instructions
        index_to_instr[instr.output] = instr
        for k = 1:arity(instr.op)
            should_use_index_not_reference(instr.op, k) && continue
            if start_from_inputs
                SimpleGraphs.add_edges!(g, [(instr.input[k], instr.output)])
            else
                SimpleGraphs.add_edges!(g, [(instr.output, instr.input[k])])
            end
        end
    end
    g, index_to_instr
end



# Reduce tape space

function reduce_space(instructions::Vector{Instruction}, input_block_size::Int, assignments)
    index_map, space_needed, updated_assignments =
        index_compactification_mapping(instructions, input_block_size, assignments)

    instructions = map(instructions) do instr
        Instruction(
            ntuple(Val(4)) do k

                if should_use_index_not_reference(instr.op, k)
                    instr.input[k]
                else
                    get(index_map, instr.input[k], instr.input[k])
                end
            end,
            instr.op,
            index_map[instr.output],
        )
    end

    instructions, space_needed, updated_assignments
end

"""
    index_compactification_mapping(instructions::Vector{Instruction}, input_block_size::Int, ignored::Set{Int})

Assemble a mapping of instruction indices such that no instruction values are overwritten
but the number of different instruction output indices is reduced. `ignored` is a set or range
of indices whose output index will not be changed
"""
function index_compactification_mapping(
    instructions::Vector{Instruction},
    input_block_size::Int,
    assignments,
)
    used_indices = Set{Int32}()
    unused_indices = Vector{Int32}()

    index_map = Dict{Int32,Int32}()
    for idx = 1:input_block_size
        index_map[idx] = idx
    end

    get_register!() = begin
        if isempty(unused_indices)
            push!(unused_indices, input_block_size + length(used_indices) + 1)
        end
        r = pop!(unused_indices)
        push!(used_indices, r)
        r
    end

    # compute lifetime intervals of output
    output_lifetime_end = Dict{Int32,Int32}()
    for instr_idx = length(instructions):-1:1
        instr = instructions[instr_idx]
        for i = 1:arity(instr.op)
            should_use_index_not_reference(instr.op, i) && continue
            idx = instr.input[i]
            if !haskey(output_lifetime_end, idx)
                output_lifetime_end[idx] = instr_idx
            end
        end
    end

    for (instr_idx, instr) in enumerate(instructions)
        if instr.output ∉ assignments
            index_map[instr.output] = get_register!()
        end

        for i = 1:arity(instr.op)
            should_use_index_not_reference(instr.op, i) && continue

            input_idx = instr.input[i]
            if input_idx > input_block_size &&
               instr_idx == get(output_lifetime_end, input_idx, 0)

                mapped_idx = get(index_map, input_idx, input_idx)
                if (mapped_idx ∈ used_indices)
                    pop!(used_indices, mapped_idx)
                    push!(unused_indices, mapped_idx)
                end
            end
        end
    end

    max_indices_used = length(unused_indices)
    # update assignments and also add to index map
    updated_assignments =
        range(input_block_size + max_indices_used + 1; length = length(assignments))
    for (k, idx) in enumerate(assignments)
        index_map[idx] = input_block_size + max_indices_used + k
    end

    index_compactification_mapping = last(updated_assignments)

    index_map, index_compactification_mapping, updated_assignments
end
