# Predictor-only Taylor storage. Unlike the ordinary interpreter, every
# instruction retains its output, so lower coefficients survive the next pass.
mutable struct IncrementalTaylorInterpreter{T}
    interpreter::Interpreter{Vector{TruncatedTaylorSeries{4,T}}}
    order::Int
end

function IncrementalTaylorInterpreter(::Type{T}, I::Interpreter) where {T}
    s = I.sequence
    n = last(s.variables_range)
    slots = zeros(Int32, max(n, s.tape_space_needed))
    slots[1:n] .= 1:n
    instructions = Instruction[]
    for instr in s.instructions
        instr.op == OP_STOP && break
        args = ntuple(4) do j
            a = instr.input[j]
            j > arity(instr.op) || should_use_index_not_reference(instr.op, j) ? a :
            slots[a]
        end
        n += 1
        slots[instr.output] = n
        push!(instructions, Instruction(args, instr.op, n))
    end
    push!(instructions, Instruction((n, n, n, n), OP_STOP, n))
    assignments = [(i, Int(slots[k])) for (i, k) in s.assignments]
    ua = [(i, Int(slots[k])) for (i, k) in s.u_assignments]
    ja = [(i, Int(slots[k])) for (i, k) in s.U_assignments]
    seq = InstructionSequence(
        instructions,
        s.constants,
        s.constants_range,
        s.parameters_range,
        s.variables_range,
        s.continuation_parameter_index,
        assignments,
        s.output_dim,
        n,
        ua,
        ja,
        s.all_u_assigned,
        s.all_U_assigned,
    )
    tape = interpreter(
        Vector{TruncatedTaylorSeries{4,T}},
        seq;
        variables = I.variables,
        parameters = I.parameters,
    )
    IncrementalTaylorInterpreter{T}(tape, 0)
end

# Reuse the existing Taylor operation generator. Only store the newly affected
# coefficients; inlining lets LLVM discard the unused lower-order arithmetic.
function partial_taylor_stores(ex, K, L)
    ex isa Expr || return ex
    if ex.head == :(=) && ex.args[1] == :(tape[i])
        result = gensym(:series)
        values = [
            k < L ? :(tape[i][$k]) : k <= K ? :($result[$k]) : :(zero($result[0])) for
            k = 0:3
        ]
        return quote
            $result = $(ex.args[2])
            tape[i] = TruncatedTaylorSeries(($(values...),))
        end
    end
    Expr(ex.head, map(a -> partial_taylor_stores(a, K, L), ex.args)...)
end

@generated function execute_incremental_taylor!(
    ::Val{K},
    ::Val{L},
    x,
    p,
    tape,
    sequence,
) where {K,L}
    partial_taylor_stores(execute_taylor_instructions!_impl(K), K, L)
end

function predictor_taylor!(
    u,
    ::Val{K},
    cache::IncrementalTaylorInterpreter,
    x,
    t,
    p,
    incremental,
) where {K}
    1 <= K <= 3 || throw(ArgumentError("predictor Taylor order must be 1, 2, or 3"))
    I = cache.interpreter
    s, tape = I.sequence, I.tape
    checkbounds(u, 1:s.output_dim)
    checkbounds(x, 1:length(s.variables_range))
    p === nothing || checkbounds(p, 1:length(s.parameters_range))
    # Only coefficients below the new order survive incremental execution.
    # Compare values before overwriting the input slots, including mutations of
    # an existing TaylorVector; object identity alone cannot establish validity.
    reuse = incremental && K > 1 && cache.order == K - 1
    @inbounds for (i, k) in enumerate(s.variables_range)
        previous = tape[k]
        tape[k] = x[i]
        reuse &= all(j -> isequal(previous[j], tape[k][j]), 0:K-2)
    end
    @inbounds for (i, k) in enumerate(s.parameters_range)
        previous = tape[k]
        tape[k] = p[i]
        reuse &= all(j -> isequal(previous[j], tape[k][j]), 0:K-2)
    end
    if s.continuation_parameter_index !== nothing
        k = s.continuation_parameter_index
        previous = tape[k]
        tape[k] = t
        reuse &= all(j -> isequal(previous[j], tape[k][j]), 0:K-2)
    end
    if reuse
        execute_incremental_taylor!(Val(K), Val(K - 1), x, p, tape, s)
    else
        execute_incremental_taylor!(Val(K), Val(0), x, p, tape, s)
    end
    cache.order = K
    fill!(u, 0)
    @inbounds for (i, k) in s.u_assignments
        u[i] = tape[k][K]
    end
    u
end
