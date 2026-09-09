@testset "test operations against taylor operations (N=$N, K=$K)" for (K, N) in [
    (3, 4),
    (3, 2),
    (2, 1),
]
    @var x[0:5] y[0:5] z[0:5] w[0:5]

    tx = ModelKit.TruncatedTaylorSeries(ModelKit.Expression.(x[1:N+1]))
    ty = ModelKit.TruncatedTaylorSeries(ModelKit.Expression.(y[1:N+1]))
    tz = ModelKit.TruncatedTaylorSeries(ModelKit.Expression.(z[1:N+1]))
    tw = ModelKit.TruncatedTaylorSeries(ModelKit.Expression.(w[1:N+1]))

    function normalized_taylor_term(expr, order, ε = Variable(:ε))
        isnothing(expr) && return nothing
        expand(subs(differentiate(expr, ε, order), ε => 0) / factorial(order))
    end

    function check_correctness(op)
        op_f = getfield(ModelKit, ModelKit.op_call(op))
        taylor_op_f = getfield(ModelKit, ModelKit.taylor_op_call(op))
        f =
            (val, args...) -> begin
                taylor_res = taylor_op_f(val, args...)
                expected_res = ModelKit.TruncatedTaylorSeries([
                    normalized_taylor_term(op_f(ModelKit.expression.(args)...), k)
                    for k = 0:K
                ])
                (expand.(taylor_res), expand.(expected_res))
            end

        (taylor_res, expected_res) = if ModelKit.arity(op) === 0

            taylor_op_f(Val(K)), nothing
        elseif ModelKit.arity(op) === 1
            f(Val(K), tx)
        elseif ModelKit.arity(op) === 2
            if op === ModelKit.OP_POW_INT
                f(Val(K), tx, 5)
            elseif op === ModelKit.OP_POW
                f(Val(K), tx, 1.5)
            else
                f(Val(K), tx, ty)
            end
        elseif ModelKit.arity(op) === 3
            f(Val(K), tx, ty, tz)
        elseif ModelKit.arity(op) === 4
            f(Val(K), tx, ty, tz, tw)
        end

        @test taylor_res == expected_res
    end

    ops = instances(ModelKit.OpType)

    @testset "op: $(op)" for op in ops
        check_correctness(op)
    end

end

@testset "Integer Taylor powers with zero constant coefficient" begin
    for pow in (ModelKit.taylor_op_pow_int, ModelKit.taylor_op_pow)
        for r in (0, 1, 2, 5, 8)
            expected = ntuple(k -> k == r + 1 ? 1.0 : 0.0, 6)
            @test Tuple(pow(Val(5), (0.0, 1.0), r)) == expected
        end
        # (2ε + 3ε²)^3 = 8ε³ + 36ε⁴ + 54ε⁵ + 27ε⁶
        @test Tuple(pow(Val(6), (0.0, 2.0, 3.0), 3)) ==
              (0.0, 0.0, 0.0, 8.0, 36.0, 54.0, 27.0)
        @test Tuple(pow(Val(2), (0.0im, 1.0im), 2)) == (0.0, 0.0, -1.0)
        @test Tuple(pow(Val(0), (0.0, 1.0), 0)) == (1.0,)
        @test Tuple(pow(Val(2), (0.0,), 2)) == (0.0, 0.0, 0.0)
    end
end
