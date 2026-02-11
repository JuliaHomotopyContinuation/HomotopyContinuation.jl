@testset "test operations against taylor operations (N=$N, K=$K)" for (K, N) in [
        (3, 4),
        (3, 2),
        (2, 1),
    ]
    @var x[0:5] y[0:5] z[0:5] w[0:5]

    tx = ModelKit.TruncatedTaylorSeries(ModelKit.Expression.(x[1:(N + 1)]))
    ty = ModelKit.TruncatedTaylorSeries(ModelKit.Expression.(y[1:(N + 1)]))
    tz = ModelKit.TruncatedTaylorSeries(ModelKit.Expression.(z[1:(N + 1)]))
    tw = ModelKit.TruncatedTaylorSeries(ModelKit.Expression.(w[1:(N + 1)]))

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
            expected_res = ModelKit.TruncatedTaylorSeries(
                [
                    normalized_taylor_term(op_f(ModelKit.expression.(args)...), k)
                        for k in 0:K
                ]
            )
            (expand.(taylor_res), expand.(expected_res))
        end

        (taylor_res, expected_res) = if ModelKit.arity(op) === 0

            taylor_op_f(Val(K)), nothing
        elseif ModelKit.arity(op) === 1
            f(Val(K), tx)
        elseif ModelKit.arity(op) === 2
            if op === ModelKit.OP_POW_INT
                f(Val(K), tx, 5)
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
