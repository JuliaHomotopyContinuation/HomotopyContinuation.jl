@testset "ModelKit - Instructions" begin
    @testset "evaluate" begin
        function test_instr_interpreter(instr, vars, expected)
            list = ModelKit.InstructionList()
            out = push!(list, ModelKit.Instruction(instr, vars...))
            I = ModelKit.Interpreter(list, [out]; variables = vars)
            cache =
                ModelKit.InterpreterCache(zeros(Expression, ModelKit.cache_min_length(I)))
            u = Expression[0]
            ModelKit.execute!(u, I, Variable.(vars), nothing, cache)
            @test expand(u[1] - expected) == 0
        end
        function test_instr_pow_interpreter(instr, r, expected)
            @var x
            list = ModelKit.InstructionList()
            out = push!(list, ModelKit.Instruction(instr, :x, r))
            I = ModelKit.Interpreter(list, [out]; variables = [:x])
            cache =
                ModelKit.InterpreterCache(zeros(Expression, ModelKit.cache_min_length(I)))
            u = Expression[0]
            ModelKit.execute!(u, I, [x], nothing, cache)
            @test expand(u[1] - expected) == 0
        end

        @var x y z
        test_instr_interpreter(ModelKit.INSTR_NEG, [:x], -x)
        test_instr_interpreter(ModelKit.INSTR_SQR, [:x], x^2)
        test_instr_interpreter(ModelKit.INSTR_SIN, [:x], sin(x))
        test_instr_interpreter(ModelKit.INSTR_COS, [:x], cos(x))

        test_instr_pow_interpreter(ModelKit.INSTR_POW, 2, x^2)
        test_instr_pow_interpreter(ModelKit.INSTR_POW, 5, x^5)
        test_instr_pow_interpreter(ModelKit.INSTR_POW, 6, x^6)
        test_instr_interpreter(ModelKit.INSTR_ADD, [:x, :y], x + y)
        test_instr_interpreter(ModelKit.INSTR_MUL, [:x, :y], x * y)
        test_instr_interpreter(ModelKit.INSTR_SUB, [:x, :y], x - y)
        test_instr_interpreter(ModelKit.INSTR_DIV, [:x, :y], x / y)

        test_instr_interpreter(ModelKit.INSTR_MULADD, [:x, :y, :z], x * y + z)
        test_instr_interpreter(ModelKit.INSTR_MULSUB, [:x, :y, :z], x * y - z)
        test_instr_interpreter(ModelKit.INSTR_SUBMUL, [:x, :y, :z], z - x * y)
    end

    @testset "evaluate with arb" begin
        function test_instr_interpreter(instr, vars, values, expected)
            list = ModelKit.InstructionList()
            out = push!(list, ModelKit.Instruction(instr, vars...))
            I = ModelKit.Interpreter(list, [out]; variables = vars)
            cache = ModelKit.InterpreterCache(Arblib.Acb, I; prec = 256)
            u = Arblib.AcbVector(1; prec = 256)
            ModelKit.execute!(u, I, values, nothing, cache)
            @test Float64(Arblib.get!(Arblib.Mag(), u[1] - expected)) < 1e-14
        end
        function test_instr_pow_interpreter(instr, r, x, expected)
            list = ModelKit.InstructionList()
            out = push!(list, ModelKit.Instruction(instr, :x, r))
            I = ModelKit.Interpreter(list, [out]; variables = [:x])
            cache = ModelKit.InterpreterCache(Arblib.Acb, I; prec = 256)
            u = Arblib.AcbVector(1; prec = 256)
            ModelKit.execute!(u, I, [x], nothing, cache)
            @test Float64(Arblib.get!(Arblib.Mag(), u[1] - expected)) < 1e-14
        end

        x, y, z = Arblib.Acb.(randn(ComplexF64, 3); prec = 256)
        test_instr_interpreter(ModelKit.INSTR_NEG, [:x], [x], -x)
        test_instr_interpreter(ModelKit.INSTR_SQR, [:x], [x], x^2)
        test_instr_interpreter(ModelKit.INSTR_SIN, [:x], [x], sin(x))
        test_instr_interpreter(ModelKit.INSTR_COS, [:x], [x], cos(x))

        test_instr_pow_interpreter(ModelKit.INSTR_POW, 2, x, x^2)
        test_instr_pow_interpreter(ModelKit.INSTR_POW, 5, x, x^5)
        test_instr_pow_interpreter(ModelKit.INSTR_POW, 6, x, x^6)
        test_instr_interpreter(ModelKit.INSTR_ADD, [:x, :y], [x, y], x + y)
        test_instr_interpreter(ModelKit.INSTR_MUL, [:x, :y], [x, y], x * y)
        test_instr_interpreter(ModelKit.INSTR_SUB, [:x, :y], [x, y], x - y)
        test_instr_interpreter(ModelKit.INSTR_DIV, [:x, :y], [x, y], x / y)

        test_instr_interpreter(ModelKit.INSTR_MULADD, [:x, :y, :z], [x, y, z], x * y + z)
        test_instr_interpreter(ModelKit.INSTR_MULSUB, [:x, :y, :z], [x, y, z], x * y - z)
        test_instr_interpreter(ModelKit.INSTR_SUBMUL, [:x, :y, :z], [x, y, z], z - x * y)
    end

    @testset "gradient" begin
        function test_instr_interpreter_grad(instr, vars, expected)
            list = ModelKit.InstructionList()
            out = push!(list, ModelKit.Instruction(instr, vars...))
            grad_list, ev_out, jac_out = ModelKit.gradient(list, [out], vars)
            I = ModelKit.Interpreter(
                grad_list,
                [ev_out; vec(jac_out)];
                variables = vars,
                nexpressions = 1,
            )
            cache =
                ModelKit.InterpreterCache(zeros(Expression, ModelKit.cache_min_length(I)))
            u = zeros(Expression, 1 + length(vars))
            ModelKit.execute!(u, I, Variable.(vars), nothing, cache)
            @test expand.(u - expected) == zeros(Expression, 1 + length(vars))
        end
        function test_instr_pow_interpreter_grad(instr, r, expected)
            @var x
            list = ModelKit.InstructionList()
            out = push!(list, ModelKit.Instruction(instr, :x, r))
            grad_list, ev_out, jac_out = ModelKit.gradient(list, [out], [:x])
            I = ModelKit.Interpreter(
                grad_list,
                [ev_out; vec(jac_out)];
                variables = [:x],
                nexpressions = 1,
            )
            cache =
                ModelKit.InterpreterCache(zeros(Expression, ModelKit.cache_min_length(I)))
            u = Expression[0, 0]
            ModelKit.execute!(u, I, [x], nothing, cache)
            @test expand.(u - expected) == Expression[0, 0]
        end
        @var x y z

        test_instr_interpreter_grad(ModelKit.INSTR_NEG, [:x], [-x, -1])
        test_instr_interpreter_grad(ModelKit.INSTR_SQR, [:x], [x^2, 2x])
        test_instr_interpreter_grad(ModelKit.INSTR_SIN, [:x], [sin(x), cos(x)])
        test_instr_interpreter_grad(ModelKit.INSTR_COS, [:x], [cos(x), -sin(x)])

        test_instr_pow_interpreter_grad(ModelKit.INSTR_POW, 2, [x^2, 2x])
        test_instr_pow_interpreter_grad(ModelKit.INSTR_POW, 5, [x^5, 5x^4])
        test_instr_pow_interpreter_grad(ModelKit.INSTR_POW, 6, [x^6, 6x^5])
        test_instr_interpreter_grad(ModelKit.INSTR_ADD, [:x, :y], [x + y, 1, 1])
        test_instr_interpreter_grad(ModelKit.INSTR_MUL, [:x, :y], [x * y, y, x])
        test_instr_interpreter_grad(ModelKit.INSTR_SUB, [:x, :y], [x - y, 1, -1])
        test_instr_interpreter_grad(ModelKit.INSTR_DIV, [:x, :y], [x / y, 1 / y, -x / y^2])

        test_instr_interpreter_grad(
            ModelKit.INSTR_MULADD,
            [:x, :y, :z],
            [x * y + z, y, x, 1],
        )
        test_instr_interpreter_grad(
            ModelKit.INSTR_MULSUB,
            [:x, :y, :z],
            [x * y - z, y, x, -1],
        )
        test_instr_interpreter_grad(
            ModelKit.INSTR_SUBMUL,
            [:x, :y, :z],
            [z - x * y, -y, -x, 1],
        )
    end

    @testset "taylor" begin
        function test_instr_taylor1(instr, op, K, M)
            vars = [:x]
            list = ModelKit.InstructionList()
            out = push!(list, ModelKit.Instruction(instr, vars...))
            I = ModelKit.Interpreter(list, [out]; variables = vars)
            tape = ModelKit.TaylorVector{K + 1}(Expression, ModelKit.cache_min_length(I))
            cache = ModelKit.InterpreterCache(tape)
            u = ModelKit.TaylorVector{K + 1}(Expression, 1)
            @var x[0:M]
            X = ModelKit.TaylorVector{M + 1}(Expression, 1)
            X[1] = tuple(x...)
            ModelKit.execute!(u, Val(K), I, X, nothing, cache)

            @var λ
            tx = sum(λ .^ (k - 1) .* xk for (k, xk) in enumerate(x))
            true_value = [differentiate(op(tx), λ, k)(λ => 0) / factorial(k) for k = 0:K]
            @test expand.(u[1] .- true_value) == zeros(Expression, K + 1)
        end

        function test_instr_taylor_pow(r, K, M)
            vars = [:x]
            list = ModelKit.InstructionList()
            out = push!(list, ModelKit.Instruction(ModelKit.INSTR_POW, :x, r))
            I = ModelKit.Interpreter(list, [out]; variables = vars)
            tape = ModelKit.TaylorVector{K + 1}(Expression, ModelKit.cache_min_length(I))
            cache = ModelKit.InterpreterCache(tape)
            u = ModelKit.TaylorVector{K + 1}(Expression, 1)
            @var x[0:M]
            X = ModelKit.TaylorVector{M + 1}(Expression, 1)
            X[1] = tuple(x...)
            ModelKit.execute!(u, Val(K), I, X, nothing, cache)

            @var λ
            tx = sum(λ .^ (k - 1) .* xk for (k, xk) in enumerate(x))
            true_value = [differentiate(tx^r, λ, k)(λ => 0) / factorial(k) for k = 0:K]
            @test expand.(u[1] .- true_value) == zeros(Expression, K + 1)
        end
        function test_instr_taylor2(instr, op, K, M)
            vars = [:x, :y]
            list = ModelKit.InstructionList()
            out = push!(list, ModelKit.Instruction(instr, vars...))
            I = ModelKit.Interpreter(list, [out]; variables = vars)
            tape = ModelKit.TaylorVector{K + 1}(Expression, ModelKit.cache_min_length(I))
            cache = ModelKit.InterpreterCache(tape)
            u = ModelKit.TaylorVector{K + 1}(Expression, 1)
            @var x[0:M] y[0:M]
            X = ModelKit.TaylorVector{M + 1}(Expression, 2)
            X[1] = tuple(x...)
            X[2] = tuple(y...)
            ModelKit.execute!(u, Val(K), I, X, nothing, cache)

            @var λ
            tx = sum(λ .^ (k - 1) .* xk for (k, xk) in enumerate(x))
            ty = sum(λ .^ (k - 1) .* yk for (k, yk) in enumerate(y))
            true_value =
                [differentiate(op(tx, ty), λ, k)(λ => 0) / factorial(k) for k = 0:K]
            @test expand.(u[1] .- true_value) == zeros(Expression, K + 1)
        end

        function test_instr_taylor3(instr, op, K, M)
            vars = [:x, :y, :z]
            list = ModelKit.InstructionList()
            out = push!(list, ModelKit.Instruction(instr, vars...))
            I = ModelKit.Interpreter(list, [out]; variables = vars)
            tape = ModelKit.TaylorVector{K + 1}(Expression, ModelKit.cache_min_length(I))
            cache = ModelKit.InterpreterCache(tape)
            u = ModelKit.TaylorVector{K + 1}(Expression, 1)
            @var x[0:M] y[0:M] z[0:M]
            X = ModelKit.TaylorVector{M + 1}(Expression, 3)
            X[1] = tuple(x...)
            X[2] = tuple(y...)
            X[3] = tuple(z...)
            ModelKit.execute!(u, Val(K), I, X, nothing, cache)

            @var λ
            tx = sum(λ .^ (k - 1) .* xk for (k, xk) in enumerate(x))
            ty = sum(λ .^ (k - 1) .* yk for (k, yk) in enumerate(y))
            tz = sum(λ .^ (k - 1) .* yk for (k, yk) in enumerate(z))
            true_value =
                [differentiate(op(tx, ty, tz), λ, k)(λ => 0) / factorial(k) for k = 0:K]
            @test expand.(u[1] .- true_value) == zeros(Expression, K + 1)
        end


        for K = 0:4, M = 0:K
            test_instr_taylor1(ModelKit.INSTR_NEG, -, K, M)
            test_instr_taylor1(ModelKit.INSTR_SQR, a -> a^2, K, M)
            test_instr_taylor1(ModelKit.INSTR_SIN, sin, K, M)
            test_instr_taylor1(ModelKit.INSTR_COS, cos, K, M)
            test_instr_taylor2(ModelKit.INSTR_ADD, +, K, M)
            test_instr_taylor2(ModelKit.INSTR_SUB, -, K, M)
            test_instr_taylor2(ModelKit.INSTR_MUL, *, K, M)
            test_instr_taylor2(ModelKit.INSTR_DIV, /, K, M)
            test_instr_taylor_pow(2, K, M)
            test_instr_taylor_pow(4, K, M)
            test_instr_taylor_pow(7, K, M)
            test_instr_taylor3(ModelKit.INSTR_MULADD, muladd, K, M)
            test_instr_taylor3(ModelKit.INSTR_MULSUB, (x, y, z) -> x * y - z, K, M)
            test_instr_taylor3(ModelKit.INSTR_SUBMUL, (x, y, z) -> z - x * y, K, M)
        end
    end
end
