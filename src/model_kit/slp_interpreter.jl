@enum InterpreterDataSource::UInt8 begin
    DATA_TAPE # variables, t in homotopy, instruction results
    DATA_PARAMETERS
    DATA_CONSTANTS
end

struct InterpreterArg
    data::InterpreterDataSource
    ind::Int32
end
function Base.show(io::IO, arg::InterpreterArg)
    if arg.data == DATA_TAPE
        print(io, "tape_", arg.ind)
    elseif arg.data == DATA_PARAMETERS
        print(io, "p_", arg.ind)
    elseif arg.data == DATA_CONSTANTS
        print(io, "c_", arg.ind)
    end
end

struct InterpreterInstruction
    op::InstructionOp
    args::NTuple{3,InterpreterArg}
    tape_index::Int
end

function Base.show(io::IO, instr::InterpreterInstruction)
    print(
        io,
        :(
            tape[$(instr.tape_index)] =
                $(Expr(:call, Symbol(instr.op), Symbol.(sprint.(show, instr.args))...))
        ),
    )
end

struct Interpreter{T}
    instructions::Vector{InterpreterInstruction}
    constants::Vector{T}
    nexpressions::Int
    assignments::Vector{InterpreterArg}
    variables::Vector{Symbol}
    parameters::Vector{Symbol}
end

function Base.show(io::IO, I::Interpreter{T}) where {T}
    print(io, "Interpreter{$T} for $(length(I.instructions)) instructions")
end
struct InterpreterCache{A<:AbstractArray}
    tape::A
end
Base.length(C::InterpreterCache) = length(C.tape)
cache_min_length(I::Interpreter) = maximum(instr -> instr.tape_index, I.instructions)
InterpreterCache(T::Type, I::Interpreter) = InterpreterCache(zeros(T, cache_min_length(I)))

function InterpreterArg(instr_map, var_map, t, param_map, constants, arg)
    if haskey(instr_map, arg)
        return InterpreterArg(DATA_I, instr_map[arg])
    elseif haskey(var_map, arg)
        return InterpreterArg(DATA_X, var_map[arg])
    elseif haskey(param_map, arg)
        return InterpreterArg(DATA_P, param_map[arg])
    elseif t == arg
        return InterpreterArg(DATA_T, 0)
    else # constant
        c = findfirst(c -> c == arg, constants)
        if isnothing(c)
            push!(constants, arg)
            return InterpreterArg(DATA_C, length(constants))
        else
            return InterpreterArg(DATA_C, c)
        end
    end
end

function Interpreter(
    list::InstructionList,
    assignments;
    variables::Vector{Symbol},
    t::Union{Nothing,Symbol} = nothing,
    parameters::Vector{Symbol} = Symbol[],
    nexpressions::Int = length(assignments),
)
    # compress needed storage
    assignment_tape_indices = Set{Int}()
    for a in assignments
        if a isa InstructionRef
            push!(assignment_tape_indices, a.i)
        end
    end
    # alloc_alias = allocation_alias(list, assignment_tape_indices)
    alloc_alias = collect(1:length(list))
    # we shift all allocation alias by the number of variables (plus t if applicable)
    # since we will always put the variables at the start of the tape
    shift = length(variables) + !isnothing(t)
    alloc_alias .= alloc_alias .+ shift

    var_map = Dict(zip(variables, 1:length(variables)))
    if !isnothing(t)
        var_map[t] = shift
    end
    param_map = Dict(zip(parameters, 1:length(parameters)))

    constants_map = Dict{Any,Int}()
    nconstants = 0
    constants_type = Int32
    instructions = map(list) do (instr_ref, instr)
        args = ntuple(Val(3)) do i
            arg = instr.args[i]
            if isnothing(arg)
                # dummy data
                InterpreterArg(DATA_TAPE, Int32(shift + 1))
            elseif instr.op == INSTR_POW && i == 2
                # We store the power in the INDEX of the argument
                InterpreterArg(DATA_TAPE, Int32(arg))
            elseif haskey(var_map, arg)
                InterpreterArg(DATA_TAPE, Int32(var_map[arg]))
            elseif haskey(param_map, arg)
                InterpreterArg(DATA_PARAMETERS, Int32(param_map[arg]))
            elseif arg isa InstructionRef
                InterpreterArg(DATA_TAPE, Int32(alloc_alias[arg.i]))
            else
                if haskey(constants_map, arg)
                    InterpreterArg(DATA_CONSTANTS, constants_map[arg])
                else
                    nconstants += 1
                    push!(constants_map, arg => nconstants)
                    constants_type = promote_type(typeof(arg), constants_type)
                    InterpreterArg(DATA_CONSTANTS, nconstants)
                end
            end
        end
        InterpreterInstruction(instr.op, args, alloc_alias[instr_ref.i])
    end
    constants = Vector{constants_type}(undef, nconstants)
    for (k, ind) in constants_map
        constants[ind] = k
    end
    assignment_instructions = InterpreterArg[]
    for a in assignments
        if haskey(var_map, a)
            push!(assignment_instructions, InterpreterArg(DATA_TAPE, var_map[a]))
        elseif haskey(param_map, a)
            push!(assignment_instructions, InterpreterArg(DATA_PARAMETERS, param_map[a]))
        elseif haskey(constants_map, a)
            push!(assignment_instructions, InterpreterArg(DATA_CONSTANTS, constants_map[a]))
        elseif a isa InstructionRef
            push!(assignment_instructions, InterpreterArg(DATA_TAPE, alloc_alias[a.i]))
        elseif a isa Number
            nconstants += 1
            push!(constants_map, a => nconstants)
            push!(constants, a)
            push!(assignment_instructions, InterpreterArg(DATA_CONSTANTS, nconstants))
        else
            error(string("Unhandled assignment: ", sprint(show, a)))
        end
    end
    Interpreter(
        instructions,
        constants,
        nexpressions,
        assignment_instructions,
        variables,
        parameters,
    )
end

function interpreter(F::Union{System,Homotopy})
    vars = Symbol.(variables(F))
    params = Symbol.(parameters(F))
    t = F isa Homotopy ? Symbol(F.t) : nothing
    list, out = instruction_list(expressions(F))
    ModelKit.Interpreter(list, out; variables = vars, parameters = params, t = t)
end
function jacobian_interpreter(F::Union{System,Homotopy})
    vars = Symbol.(variables(F))
    params = Symbol.(parameters(F))
    t = F isa Homotopy ? Symbol(F.t) : nothing
    list, out = instruction_list(expressions(F))
    nexpressions = length(out)
    list, out, jac_out = gradient(list, out, vars)
    out′ = [out; vec(jac_out)]
    ModelKit.Interpreter(
        list,
        out′;
        variables = vars,
        parameters = params,
        t = t,
        nexpressions = nexpressions,
    )
end

function load_data(f, var_name, arg, has_parameters::Bool)
    v1, v2, v3 = gensym(var_name), gensym(var_name), gensym(var_name)
    if has_parameters
        quote
            if $arg.data == DATA_TAPE
                $v1 = tape[$arg.ind]
                $(f(v1))
            elseif $arg.data == DATA_PARAMETERS
                $v2 = parameters[$arg.ind]
                $(f(v2))
            else #if $arg.data == DATA_CONSTANTS
                $v3 = constants[$arg.ind]
                $(f(v3))
            end
        end
    else
        quote
            if $arg.data == DATA_TAPE
                $v1 = tape[$arg.ind]
                $(f(v1))
            else #if $arg.data == DATA_CONSTANTS
                $v3 = constants[$arg.ind]
                $(f(v3))
            end
        end
    end
end

Base.@propagate_inbounds neg_fast!(z, i, x) = (z[i] = neg_fast(x))
Base.@propagate_inbounds sqr_fast!(z, i, x) = (z[i] = sqr_fast(x))
Base.@propagate_inbounds sin_fast!(z, i, x) = (z[i] = sin(x))
Base.@propagate_inbounds cos_fast!(z, i, x) = (z[i] = cos(x))

Base.@propagate_inbounds add_fast!(z, i, x, y) = (z[i] = add_fast(x, y))
Base.@propagate_inbounds sub_fast!(z, i, x, y) = (z[i] = sub_fast(x, y))
Base.@propagate_inbounds mul_fast!(z, i, x, y) = (z[i] = mul_fast(x, y))
Base.@propagate_inbounds div_fast!(z, i, x, y) = (z[i] = div_fast(x, y))
Base.@propagate_inbounds pow_fast!(z, i, x, k::Integer) = (z[i] = pow_fast(x, k))

Base.@propagate_inbounds muladd_fast!(z, i, u, v, w) = (z[i] = muladd_fast(u, v, w))
Base.@propagate_inbounds mulsub_fast!(z, i, u, v, w) = (z[i] = mulsub_fast(u, v, w))
Base.@propagate_inbounds submul_fast!(z, i, u, v, w) = (z[i] = submul_fast(u, v, w))

function evaluate_block(; has_parameters::Bool = false)
    # This should get compiled to a jump table by LLVM, so order doesn't matter
    load_data(:v₁, :arg₁, has_parameters) do v₁
        quote
            if op == INSTR_NEG
                neg_fast!(tape, i, $v₁)
            elseif op == INSTR_SQR
                sqr_fast!(tape, i, $v₁)
            elseif op == INSTR_SIN
                sin_fast!(tape, i, $v₁)
            elseif op == INSTR_COS
                cos_fast!(tape, i, $v₁)
            elseif op == INSTR_POW
                pow_fast!(tape, i, $v₁, arg₂.ind)
            else
                $(
                    load_data(:v₂, :arg₂, has_parameters) do v₂
                        quote
                            if op == INSTR_ADD
                                add_fast!(tape, i, $v₁, $v₂)
                            elseif op == INSTR_SUB
                                sub_fast!(tape, i, $v₁, $v₂)
                            elseif op == INSTR_MUL
                                mul_fast!(tape, i, $v₁, $v₂)
                            elseif op == INSTR_DIV
                                div_fast!(tape, i, $v₁, $v₂)
                            else
                                $(
                                    load_data(:v₃, :arg₃, has_parameters) do v₃
                                        quote
                                            if op == INSTR_MULADD
                                                muladd_fast!(tape, i, $v₁, $v₂, $v₃)
                                            elseif op == INSTR_MULSUB
                                                mulsub_fast!(tape, i, $v₁, $v₂, $v₃)
                                            elseif op == INSTR_SUBMUL
                                                submul_fast!(tape, i, $v₁, $v₂, $v₃)
                                            else
                                                throw(error(string(
                                                    "Unexpected fall-through. instr: ",
                                                    Symbol(op),
                                                )))
                                            end
                                        end
                                    end
                                )
                            end
                        end
                    end
                )
            end
        end
    end
end

function taylor_block(; has_parameters::Bool = false)
    # This should get compiled to a jump table by LLVM, so order doesn't matter
    three_arg =
        (v₁, v₂) -> load_data(:v₃, :arg₃, has_parameters) do v₃
            quote
                if op == INSTR_MULADD
                    tape[i] = taylor_muladd(Order, $v₁, $v₂, $v₃)
                elseif op == INSTR_MULSUB
                    tape[i] = taylor_mulsub(Order, $v₁, $v₂, $v₃)
                elseif op == INSTR_SUBMUL
                    tape[i] = taylor_submul(Order, $v₁, $v₂, $v₃)
                else
                    throw(error(string("Unexpected fall-through. instr: ", Symbol(op))))
                end
            end
        end

    two_arg = v₁ -> load_data(:v₂, :arg₂, has_parameters) do v₂
        quote
            if op == INSTR_ADD
                tape[i] = taylor_add(Order, $v₁, $v₂)
            elseif op == INSTR_SUB
                tape[i] = taylor_sub(Order, $v₁, $v₂)
            elseif op == INSTR_MUL
                tape[i] = taylor_mul(Order, $v₁, $v₂)
            elseif op == INSTR_DIV
                tape[i] = taylor_div(Order, $v₁, $v₂)
            else
                $(three_arg(v₁, v₂))
            end
        end
    end
    load_data(:v₁, :arg₁, has_parameters) do v₁
        quote
            if op == INSTR_NEG
                tape[i] = taylor_neg(Order, $v₁)
            elseif op == INSTR_SQR
                tape[i] = taylor_sqr(Order, $v₁)
            elseif op == INSTR_SIN
                tape[i] = taylor_sin(Order, $v₁)
            elseif op == INSTR_COS
                tape[i] = taylor_cos(Order, $v₁)
            elseif op == INSTR_POW
                tape[i] = taylor_pow(Order, $v₁, arg₂.ind)
            else
                $(two_arg(v₁))
            end
        end
    end
end

for has_parameters in [true, false], has_t in [true, false]
    @eval function execute!(
        u::AbstractArray,
        I::Interpreter,
        x::AbstractArray,
        $((has_t ? (:t,) : ())...),
        $(has_parameters ? :(parameters::AbstractArray) : :(parameters::Nothing)),
        cache::InterpreterCache,
    )
        isnothing(parameters) || checkbounds(parameters, 1:length(I.parameters))

        constants = I.constants
        tape = cache.tape
        D = length(I.variables)
        for i = 1:D
            tape[i] = x[i]
        end
        $((has_t ? (:(tape[D+1] = t),) : ())...)
        @inbounds for instr in I.instructions
            i = instr.tape_index
            op = instr.op
            arg₁, arg₂, arg₃ = instr.args
            $(evaluate_block(has_parameters = has_parameters))
        end
        for (k, assign_arg) in enumerate(I.assignments)
            $(
                load_data(:v, :assign_arg, has_parameters) do v
                    :(u[k] = $v)
                end
            )
        end
        u
    end

    @eval function execute!(
        u::Union{Nothing,AbstractVector},
        U::AbstractMatrix,
        I::Interpreter,
        x::AbstractArray,
        $((has_t ? (:t,) : ())...),
        $(has_parameters ? :(parameters::AbstractArray) : :(parameters::Nothing)),
        cache::InterpreterCache,
    )
        isnothing(parameters) || checkbounds(parameters, 1:length(I.parameters))

        constants = I.constants
        tape = cache.tape
        D = length(I.variables)
        for i = 1:D
            tape[i] = x[i]
        end
        $((has_t ? (:(tape[D+1] = t),) : ())...)
        @inbounds for instr in I.instructions
            i = instr.tape_index
            op = instr.op
            arg₁, arg₂, arg₃ = instr.args
            $(evaluate_block(has_parameters = has_parameters))
        end
        if !(u isa Nothing)
            for k = 1:I.nexpressions
                a = I.assignments[k]
                $(
                    load_data(:v, :a, has_parameters) do v
                        :(u[k] = $v)
                    end
                )
            end
        end
        i = I.nexpressions
        for l = 1:D, k = 1:I.nexpressions
            i += 1
            a = I.assignments[i]
            $(
                load_data(:v, :a, has_parameters) do v
                    :(U[k, l] = $v)
                end
            )
        end
        nothing
    end

    @eval function execute!(
        u::A,
        Order::Val{K},
        I::Interpreter,
        x::AbstractArray,
        $((has_t ? (:t,) : ())...),
        $(has_parameters ? :(parameters::AbstractArray) : :(parameters::Nothing)),
        cache::InterpreterCache,
    ) where {A<:AbstractArray,K}
        isnothing(parameters) || checkbounds(parameters, 1:length(I.parameters))

        constants = I.constants
        tape = cache.tape
        D = length(I.variables)
        for i = 1:D
            tape[i] = x[i]
        end
        $((has_t ? (:(tape[D+1] = (t, one(t))),) : ())...)
        @inbounds for instr in I.instructions
            i = instr.tape_index
            op = instr.op
            arg₁, arg₂, arg₃ = instr.args
            $(taylor_block(has_parameters = has_parameters))
        end
        if u isa TaylorVector
            for (k, assign_arg) in enumerate(I.assignments)
                $(
                    load_data(:v, :assign_arg, has_parameters) do v
                        :(u[k] = $v)
                    end
                )
            end
        else
            for (k, assign_arg) in enumerate(I.assignments)
                $(
                    load_data(:v, :assign_arg, has_parameters) do v
                        :(u[k] = last($v))
                    end
                )
            end
        end
        u
    end
end
