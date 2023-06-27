mutable struct Interpreter{V<:AbstractVector}
    sequence::InstructionSequence
    tape::V
    variables::Vector{Symbol}
    parameters::Vector{Symbol}
end

variable_indices(I::Interpreter) = I.sequence.variables_range
input_block_size(I::Interpreter) = last(variable_indices(I))

function Base.show(io::IO, I::Interpreter{T}) where {T}
    print(io, "Interpreter{$T} for $(length(I.sequence)) instructions")
end

"""
    show_instructions([io::IO], I::Interpreter)

Show the instruction executed by the interpreter `I`.
"""
show_instructions(I::Interpreter) = show_instructions(stdout, I)
function show_instructions(io::IO, I::Interpreter)
    for i = 1:I.nconstants
        println(io, "t[$i] = $(I.tape[i])")
    end
    for (i, p) in enumerate(I.parameters)
        println(io, "t[", i + I.nconstants, "] = $(p)")
    end
    for (i, x) in enumerate(I.variables)
        println(io, "t[", i + I.nconstants + length(I.parameters), "] = $(x)")
    end
    for instr in I.sequence
        println(io, instr)
    end
end


interpreter(::Type{T}, F::Union{System,Homotopy}; kwargs...) where {T} =
    interpreter(Vector{T}, F; kwargs...)

function interpreter(
    ::Type{V},
    F::Union{System,Homotopy};
    output_dim::Integer = length(F),
) where {V<:AbstractVector}
    vars = Symbol.(variables(F))
    params = Symbol.(parameters(F))
    continuation_parameter = F isa Homotopy ? Symbol(F.t) : nothing

    interpreter(
        V,
        intermediate_representation(expressions(F); output_dim = output_dim);
        variables = vars,
        parameters = params,
        continuation_parameter = continuation_parameter,
    )
end

function interpreter(
    ::Type{V},
    ir::IntermediateRepresentation;
    variables::Vector{Symbol},
    parameters::Vector{Symbol} = Symbol[],
    continuation_parameter::Union{Symbol,Nothing} = nothing,
) where {V<:AbstractVector}
    sequence = instruction_sequence(
        ir;
        variables = variables,
        parameters = parameters,
        continuation_parameter = continuation_parameter,
    )
    interpreter(V, sequence; variables = variables, parameters = parameters)
end

function interpreter(
    ::Type{V},
    sequence::InstructionSequence;
    variables::Vector{Symbol},
    parameters::Vector{Symbol} = Symbol[],
) where {V<:AbstractVector}
    # build tape
    tape = create_tape(V, sequence.tape_space_needed)
    tape[sequence.constants_range] .= sequence.constants

    Interpreter(sequence, tape, variables, parameters)
end

function interpreter(::Type{V}, i::Interpreter) where {V<:AbstractVector}
    interpreter(V, i.sequence; variables = i.variables, parameters = i.parameters)
end
function interpreter(::Type{T}, i::Interpreter) where {T}
    interpreter(Vector{T}, i)
end

function jacobian_interpreter(::Type{T}, F::System, ; kwargs...) where {T}
    JF = System(
        [F.expressions; vec(jacobian(F))],
        variables = variables(F),
        parameters = parameters(F),
    )
    interpreter(T, JF; output_dim = length(F))
end
function jacobian_interpreter(::Type{T}, H::Homotopy, ; kwargs...) where {T}
    JH = Homotopy(
        [H.expressions; vec(differentiate(H.expressions, H.variables))],
        variables(H),
        H.t,
        parameters(H),
    )
    interpreter(T, JH; output_dim = length(H))
end

function setprecision!(i::Interpreter{<:AcbRefVector}, prec::Int)
    i.tape = Arblib.setprecision(i.tape, prec)
    return i
end

Base.@propagate_inbounds zero!(v::AbstractArray) = fill!(v, 0)
zero!(v::Arblib.AcbVectorLike) = Arblib.zero!(v)
zero!(v::Arblib.AcbMatrixLike) = Arblib.zero!(v)
create_tape(::Type{V}, len) where {V<:AbstractVector} = zero!(similar(V, len))
create_tape(::Type{Arblib.AcbRefVector}, len; prec::Int = 128) = Arblib.AcbRefVector(len)

function nested_ifs(cond_body, elsebranch = nothing)
    orig_expr = Expr(:if, cond_body[1][1], cond_body[1][2])
    expr = orig_expr
    for i = 2:length(cond_body)
        push!(expr.args, Expr(:elseif, cond_body[i][1], cond_body[i][2]))
        expr = expr.args[end]
    end
    if elsebranch !== nothing
        push!(expr.args, elsebranch)
    end
    return orig_expr
end

function execute_instructions_inner!_impl(level = 0)
    op_types = [instances(OpType)...]
    arity1 = filter(op -> arity(op) == 1, op_types)
    arity2 = filter(op -> arity(op) == 2 && op !== OP_POW_INT, op_types)
    arity3 = filter(op -> arity(op) == 3, op_types)
    arity4 = filter(op -> arity(op) == 4, op_types)
    # This should get compiled to a jump table by LLVM, so order doesn't matter
    branches = [
        map(arity1) do op
            (:(op == $(op)), :(tape[i] = $(op_call(op))(t₁)))
        end
        [(:(op == OP_POW_INT), :(tape[i] = $(op_call(OP_POW_INT))(t₁, arg₂)))]
        map(arity2) do op
            (:(op == $(op)), quote
                t₂ = tape[arg₂]
                tape[i] = $(op_call(op))(t₁, t₂)
            end)
        end
        map(arity3) do op
            (:(op == $(op)), quote
                t₂ = tape[arg₂]
                t₃ = tape[arg₃]
                tape[i] = $(op_call(op))(t₁, t₂, t₃)
            end)
        end
        map(arity4) do op
            (:(op == $(op)), quote
                t₂ = tape[arg₂]
                t₃ = tape[arg₃]
                t₄ = tape[arg₄]
                tape[i] = $(op_call(op))(t₁, t₂, t₃, t₄)
            end)
        end
    ]
    # Add one level of instruction recursion
    # to speed up things
    # More would help but are too prohibitiv in compile cost
    if level < 1
        branches = map(branches) do ((cond, code))
            (cond, :($code; $(execute_instructions_inner!_impl(level + 1))))
        end
    end


    quote
        Base.@_propagate_inbounds_meta
        instr = instructions[k+=1]
        op = instr.op
        arg₁, arg₂, arg₃, arg₄ = instr.input
        i = instr.output
        t₁ = tape[arg₁]
        $(nested_ifs([
            branches
            [(:(op == OP_STOP), :(break))]
        ]))
    end
end

function execute_acb_instructions_inner!_impl(level = 0)
    op_types = [instances(OpType)...]
    arity1 = filter(op -> arity(op) == 1, op_types)
    arity2 = filter(op -> arity(op) == 2 && op !== OP_POW_INT, op_types)
    arity3 = filter(op -> arity(op) == 3, op_types)
    arity4 = filter(op -> arity(op) == 4, op_types)
    # This should get compiled to a jump table by LLVM, so order doesn't matter
    branches = [
        map(arity1) do op
            (:(op == $(op)), :($(acb_op_call(op))(tape[i], t₁, m)))
        end
        [(:(op == OP_POW_INT), :($(acb_op_call(OP_POW_INT))(tape[i], t₁, arg₂, m)))]
        map(arity2) do op
            (:(op == $(op)), quote
                t₂ = tape[arg₂]
                $(acb_op_call(op))(tape[i], t₁, t₂, m)
            end)
        end
        map(arity3) do op
            (:(op == $(op)), quote
                t₂ = tape[arg₂]
                t₃ = tape[arg₃]
                $(acb_op_call(op))(tape[i], t₁, t₂, t₃, m)
            end)
        end
        map(arity4) do op
            (:(op == $(op)), quote
                t₂ = tape[arg₂]
                t₃ = tape[arg₃]
                t₄ = tape[arg₄]
                $(acb_op_call(op))(tape[i], t₁, t₂, t₃, t₄, m)
            end)
        end
    ]
    # Add one level of instruction recursion
    # to speed up things
    # More would help but are too prohibitiv in compile cost
    if level < 1
        branches = map(branches) do ((cond, code))
            (cond, :($code; $(execute_instructions_inner!_impl(level + 1))))
        end
    end


    quote
        Base.@_propagate_inbounds_meta
        instr = instructions[k+=1]
        op = instr.op
        arg₁, arg₂, arg₃, arg₄ = instr.input
        i = instr.output
        t₁ = tape[arg₁]
        $(nested_ifs([
            branches
            [(:(op == OP_STOP), :(break))]
        ]))
    end
end

@generated function execute_instructions!(tape, instructions)
    quote
        Base.@_propagate_inbounds_meta
        k = 0
        while true
            $(execute_instructions_inner!_impl())
        end
    end
end
@generated function execute_instructions!(tape::Arblib.AcbRefVector, instructions)
    quote
        Base.@_propagate_inbounds_meta
        k = 0
        # allocate some working memory
        m = (copy(tape[1]), copy(tape[1]))
        while true
            $(execute_acb_instructions_inner!_impl(1))
        end
    end
end

for has_parameters in [true, false],
    has_continuation_parameter in [true, false],
    has_second_output in [true, false]

    @eval Base.@propagate_inbounds function execute!(
        u::Union{Nothing,AbstractArray},
        $((has_second_output ? (:(U::AbstractArray),) : ())...),
        I::Interpreter,
        x::AbstractArray,
        $((has_continuation_parameter ? (:continuation_parameter,) : ())...),
        $(has_parameters ? :(parameters::AbstractArray) : :(parameters::Nothing)),
    )
        vars_range = I.sequence.variables_range
        params_range = I.sequence.parameters_range

        checkbounds(x, 1:length(vars_range))
        isnothing(parameters) || checkbounds(parameters, 1:length(params_range))

        @inbounds for (i, k) in enumerate(params_range)
            I.tape[k] = parameters[i]
        end
        $(
            has_continuation_parameter ?
            quote
                if !isnothing(I.sequence.continuation_parameter_index)
                    I.tape[I.sequence.continuation_parameter_index] =
                        continuation_parameter
                end
            end : :()
        )
        @inbounds for (i, k) in enumerate(vars_range)
            I.tape[k] = x[i]
        end
        @inbounds execute_instructions!(I.tape, I.sequence.instructions)
        $(
            has_second_output ?
            quote
                n = I.sequence.output_dim
                zero!(U)
                idx = CartesianIndices((I.sequence.output_dim, size(U, 2)))
                if isnothing(u)
                    for (i, k) in I.sequence.assignments
                        if i > n
                            U[idx[i-n]] = I.tape[k]
                        end
                    end
                else
                    zero!(u)
                    for (i, k) in I.sequence.assignments
                        if i <= n
                            u[i] = I.tape[k]
                        else
                            U[idx[i-n]] = I.tape[k]
                        end
                    end
                end
            end : quote
                @inbounds zero!(u)
                @inbounds for (i, k) in I.sequence.assignments
                    u[i] = I.tape[k]
                end
            end
        )

        u
    end
end

# TAYLOR
function execute_taylor_instructions_inner!_impl(K)
    op_types = [instances(OpType)...]
    arity1 = filter(op -> arity(op) == 1, op_types)
    arity2 = filter(op -> arity(op) == 2 && op !== OP_POW_INT, op_types)
    arity3 = filter(op -> arity(op) == 3, op_types)
    arity4 = filter(op -> arity(op) == 4, op_types)
    # This should get compiled to a jump table by LLVM, so order doesn't matter
    branches = [
        map(arity1) do op
            (:(op == $(op)), quote
                t₁ = tape[arg₁]
                tape[i] = $(taylor_op_call(op))(Val($K), t₁)
            end)
        end
        [(
            :(op == OP_POW_INT),
            quote
                t₁ = tape[arg₁]
                tape[i] = $(taylor_op_call(OP_POW_INT))(Val($K), t₁, arg₂)
            end,
        )]
        map(arity2) do op
            (
                :(op == $(op)),
                quote
                    if arg₁ > constants_params_end
                        if arg₂ > constants_params_end
                            t₁, t₂ = tape[arg₁], tape[arg₂]
                            tape[i] = $(taylor_op_call(op))(Val($K), t₁, t₂)
                        else
                            t₁, t₂ = tape[arg₁], tape[arg₂][0]
                            tape[i] = $(taylor_op_call(op))(Val($K), t₁, t₂)
                        end
                    else
                        t₁, t₂ = tape[arg₁][0], tape[arg₂]
                        tape[i] = $(taylor_op_call(op))(Val($K), t₁, t₂)
                    end
                end,
            )
        end
        map(arity3) do op
            (:(op == $(op)), quote
                    if arg₁ ≤ constants_params_end
                        t₁, t₂, t₃ = tape[arg₁][0], tape[arg₂], tape[arg₃]
                        tape[i] = $(taylor_op_call(op))(Val($K), t₁, t₂, t₃)
                    elseif arg₂ ≤ constants_params_end
                        t₁, t₂, t₃ = tape[arg₁], tape[arg₂][0], tape[arg₃]
                        tape[i] = $(taylor_op_call(op))(Val($K), t₁, t₂, t₃)
                    elseif arg₃ ≤ constants_params_end
                        t₁, t₂, t₃ = tape[arg₁], tape[arg₂], tape[arg₃][0]
                        tape[i] = $(taylor_op_call(op))(Val($K), t₁, t₂, t₃)
                    else
                        t₁, t₂ = tape[arg₁], tape[arg₂]
                        t₃ = tape[arg₃]
                        tape[i] = $(taylor_op_call(op))(Val($K), t₁, t₂, t₃)
                    end
                end)
        end
        map(arity4) do op
            (
                :(op == $(op)),
                quote
                    t₁, t₂, t₃, t₄ = tape[arg₁], tape[arg₂], tape[arg₃], tape[arg₄]
                    tape[i] = $(taylor_op_call(op))(Val($K), t₁, t₂, t₃, t₄)
                end,
            )
        end
    ]

    quote
        Base.@_propagate_inbounds_meta
        instr = instructions[k+=1]
        op = instr.op
        arg₁, arg₂, arg₃, arg₄ = instr.input
        i = instr.output
        $(nested_ifs([
            branches
            [(:(op == OP_STOP), :(break))]
        ]))
    end
end

function execute_taylor_instructions!_impl(K)
    quote
        Base.@_propagate_inbounds_meta
        instructions = sequence.instructions
        is_p_taylor = !isnothing(p) && !isempty(p) && length(p[1]) > 1
        is_x_taylor = !isempty(x) && length(x[1]) > 1
        if !is_p_taylor && !is_x_taylor && isnothing(sequence.continuation_parameter_index)
            constants_params_end = last(sequence.variables_range)
        elseif !is_p_taylor && isnothing(sequence.continuation_parameter_index)
            constants_params_end = last(sequence.parameters_range)
        elseif !is_p_taylor && !isnothing(sequence.continuation_parameter_index)
            constants_params_end = max(sequence.continuation_parameter_index - 1, 0)
        else
            constants_params_end = last(sequence.constants_range)
        end

        k = 0
        while true
            $(execute_taylor_instructions_inner!_impl(K))
        end
    end
end

@generated function execute_taylor_instructions!(::Val{K}, x, p, tape, sequence) where {K}
    execute_taylor_instructions!_impl(K)
end

for has_parameters in [true, false], has_continuation_parameter in [true, false]

    @eval Base.@propagate_inbounds function execute_taylor!(
        u::Union{Nothing,AbstractArray},
        V::Val{K},
        I::Interpreter,
        x::AbstractArray,
        $((has_continuation_parameter ? (:continuation_parameter,) : ())...),
        $(has_parameters ? :(parameters::AbstractArray) : :(parameters::Nothing));
        assign_highest_order_only::Bool = false,
    ) where {K}
        vars_range = I.sequence.variables_range
        params_range = I.sequence.parameters_range

        checkbounds(u, 1:I.sequence.output_dim)
        checkbounds(x, 1:length(vars_range))
        isnothing(parameters) || checkbounds(parameters, 1:length(params_range))

        @inbounds for (i, k) in enumerate(params_range)
            I.tape[k] = parameters[i]
        end
        $(
            has_continuation_parameter ?
            quote
                if !isnothing(I.sequence.continuation_parameter_index)
                    I.tape[I.sequence.continuation_parameter_index] =
                        continuation_parameter
                end
            end : :()
        )
        @inbounds for (i, k) in enumerate(vars_range)
            I.tape[k] = x[i]
        end
        @inbounds execute_taylor_instructions!(V, x, parameters, I.tape, I.sequence)
        @inbounds zero!(u)
        if (assign_highest_order_only)
            @inbounds for (i, k) in I.sequence.assignments
                u[i] = I.tape[k][K]
            end
        else
            @inbounds for (i, k) in I.sequence.assignments
                u[i] = I.tape[k]
            end
        end

        u
    end
end
