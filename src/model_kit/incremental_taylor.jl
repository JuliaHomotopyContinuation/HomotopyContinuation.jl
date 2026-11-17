# Incremental Taylor evaluation for the predictor.
#
# The predictor evaluates a homotopy's Taylor series at orders k = 1, 2, 3 at the same point, 
# where the call of order k only changes coefficient k-1 of the input.
# The standard interpreter overwrites its tape each time, so order 3 recomputes orders 0-2. 
# Here every instruction gets its own tape slot (no slot is reused), so the lower-order coefficients computed by the previous call survive and only the coefficients L..K are computed by `execute_incremental_taylor!`.
# Only orders K = 1, 2, 3 are supported.

mutable struct IncrementalTaylorInterpreter{T}
    interpreter::Interpreter{Vector{TruncatedTaylorSeries{4,T}}}
    order::Int
end

function IncrementalTaylorInterpreter(::Type{T}, I::Interpreter) where {T}
    s = I.sequence
    n = last(s.variables_range)
    # slots[a] is the new tape position of the value that was at position a.
    # Variables, parameters and constants keep their positions.
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
        # every instruction gets a fresh tape slot
        n += 1
        slots[instr.output] = n
        push!(instructions, Instruction(args, instr.op, n))
    end
    # terminator expected by the instruction loop
    push!(instructions, Instruction((n, n, n, n), OP_STOP, n))
    # the output assignments point to the new slots
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

# Reuse the existing Taylor operation generator. 
# Only store the newly affected coefficients.
# K is the target order and L the first coefficient to (re)compute.
# The coefficients k < L are kept from the old tape, L ≤ k ≤ K are computed and k > K are set to zero.
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

# computes Taylor coefficients L..K
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

# Whether the coefficients 0..K-2 of two series agree. `isequal`, so NaN equals NaN.
# This is a function instead of an `all(j -> ..., 0:K-2)` closure on purpose: the closure captures variables that are mutated, which boxes them and allocates on each call.
@inline function same_lower_coefficients(a, b, ::Val{K}) where {K}
    for j = 0:K-2
        isequal(a[j], b[j]) || return false
    end
    true
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
    # Reuse the tape of the previous call only if (a) that call had order K-1 and (b) the variables, parameters and t are unchanged in the coefficients 0..K-2.
    # Only coefficients below the new order survive incremental execution.
    # Compare values before overwriting the input slots, including mutations of an existing TaylorVector.
    reuse = incremental && K > 1 && cache.order == K - 1
    @inbounds for (i, k) in enumerate(s.variables_range)
        previous = tape[k]
        tape[k] = x[i]
        reuse &= same_lower_coefficients(previous, tape[k], Val(K))
    end
    @inbounds for (i, k) in enumerate(s.parameters_range)
        previous = tape[k]
        tape[k] = p[i]
        reuse &= same_lower_coefficients(previous, tape[k], Val(K))
    end
    if s.continuation_parameter_index !== nothing
        k = s.continuation_parameter_index
        previous = tape[k]
        tape[k] = t
        reuse &= same_lower_coefficients(previous, tape[k], Val(K))
    end
    # Val(K-1): recompute the coefficients K-1 and K only; Val(0): recompute all.
    if reuse
        execute_incremental_taylor!(Val(K), Val(K - 1), x, p, tape, s)
    else
        execute_incremental_taylor!(Val(K), Val(0), x, p, tape, s)
    end
    cache.order = K
    # output entries without an assignment have coefficient K equal to zero
    fill!(u, 0)
    @inbounds for (i, k) in s.u_assignments
        u[i] = tape[k][K]
    end
    u
end
