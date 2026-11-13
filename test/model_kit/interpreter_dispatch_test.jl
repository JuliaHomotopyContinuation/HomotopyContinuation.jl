@testset "Interpreter instruction dispatch" begin
    instruction = ModelKit.Instruction
    stop = instruction((0, 0, 0, 0), ModelKit.OP_STOP, 0)
    # STOP must not read or write the tape, even when no tape is needed.
    @test isnothing(ModelKit.execute_instructions!(Float64[], [stop]))
    for T in (Float64, ComplexF64), op in instances(ModelKit.OpType)
        op == ModelKit.OP_STOP && continue
        n = ModelKit.arity(op)
        f = getfield(ModelKit, ModelKit.op_call(op))
        values = T[0.4, 0.3, 0.2, 0.1]
        # Invalid indices in unused operands catch speculative operand loads.
        inputs = ntuple(i -> Int32(i <= n ? i : -1), 4)
        if op == ModelKit.OP_POW_INT
            inputs = (Int32(1), Int32(-3), Int32(-1), Int32(-1))
        end
        expected = op == ModelKit.OP_POW_INT ? f(values[1], -3) : f(values[1:n]...)
        for output in (1, 5) # both overwriting an input and writing a fresh slot
            tape = [values; zero(T)]
            ModelKit.execute_instructions!(tape, [instruction(inputs, op, output), stop])
            @test isequal(tape[output], expected)
        end
    end
    # Consecutive instructions must see previous writes, at either loop parity.
    for count in (1, 2, 3, 8, 9)
        tape = [1.0, 2.0]
        instructions = fill(instruction((1, 2, -1, -1), ModelKit.OP_ADD, 1), count)
        push!(instructions, stop)
        ModelKit.execute_instructions!(tape, instructions)
        @test tape == [1.0 + 2count, 2.0]
    end
end
