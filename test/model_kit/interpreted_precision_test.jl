@testset "Interpreted arithmetic follows input precision" begin
    D = HC.DoubleDouble.DoubleF64
    C = HC.DoubleDouble.ComplexDF64
    δ = D(2.0^-60)
    @var x p t
    F = InterpretedSystem(System([(x - 1)^2]; variables = [x]))
    u, J = zeros(C, 1), zeros(C, 1, 1)
    z = C[1+δ]
    @test isnothing(F.jac_ComplexDF64)
    evaluate_and_jacobian!(u, J, F, z)
    @test u[1] == δ^2
    @test J[1] == 2δ
    J .= 0
    jacobian!(J, F, z)
    @test J[1] == 2δ
    cache = F.jac_ComplexDF64
    evaluate_and_jacobian!(u, J, F, z)
    @test F.jac_ComplexDF64 === cache

    # A Float64 destination must not force Float64 intermediate arithmetic.
    uf, Jf = zeros(ComplexF64, 1), zeros(ComplexF64, 1, 1)
    evaluate_and_jacobian!(uf, Jf, F, z)
    @test uf[1] == Float64(δ)^2
    @test Jf[1] == 2Float64(δ)
    tx = TaylorVector{2}(reshape(C[1+δ, 1], 2, 1))
    taylor!(u, Val(1), F, tx)
    @test u[1] == 2δ
    taylor!(u, Val(2), F, tx)
    @test u[1] == 1

    # Parameters alone can require double-double, with ordinary Float64 x.
    G = InterpretedSystem(System([(p - 1) * x]; variables = [x], parameters = [p]))
    for params in (C[1+δ], D[1+δ])
        evaluate!(u, G, [2.0], params)
        @test u[1] == 2δ
        evaluate_and_jacobian!(u, J, G, [2.0], params)
        @test J[1] == δ
        jacobian!(J, G, [2.0], params)
        @test J[1] == δ
        tx = TaylorVector{4}(reshape([2.0, 3.0, 5.0, 7.0], 4, 1))
        for k = 1:3
            taylor!(u, Val(k), G, tx, params)
            @test u[1] == (3, 5, 7)[k] * δ
        end
    end

    H = InterpretedHomotopy(Homotopy([(t - 1) * x], [x], t))
    evaluate!(u, H, [2.0], 1 + δ)
    @test u[1] == 2δ
    evaluate_and_jacobian!(u, J, H, [2.0], 1 + δ)
    @test u[1] == 2δ
    @test J[1] == δ
    jacobian!(J, H, [2.0], 1 + δ)
    @test J[1] == δ
    # At x(s)=s, t(s)=1+δ+s, the first coefficient is δ.
    tx = TaylorVector{2}(reshape([0.0, 1.0], 2, 1))
    taylor!(u, Val(1), H, tx, 1 + δ)
    @test u[1] == δ
    taylor!(u, Val(2), H, tx, 1 + δ)
    @test u[1] == 1

    # Changing precision between calls must select independent cached tapes.
    evaluate_and_jacobian!(uf, Jf, F, [2.0])
    @test uf == [1.0]
    @test Jf[1] == 2.0
    evaluate_and_jacobian!(u, J, F, z)
    @test J[1] == 2δ
end

@testset "Compiled double-double Jacobians and Taylor series are interpreted" begin
    D = HC.DoubleDouble.DoubleF64
    C = HC.DoubleDouble.ComplexDF64
    δ = D(2.0^-60)
    @var x p t
    F = CompiledSystem(System([(p - 1) * x^2]; variables = [x], parameters = [p]))
    @test isnothing(F.extended[])
    u, J = zeros(C, 1), zeros(C, 1, 1)
    evaluate_and_jacobian!(u, J, F, C[2], C[1+δ])
    @test u[1] == 4δ
    @test J[1] == 4δ
    I = F.extended[]
    @test I isa InterpretedSystem
    tx = TaylorVector{4}(reshape(C[2, 1, 0, 0], 4, 1))
    taylor!(u, Val(1), F, tx, C[1+δ])
    @test u[1] == 4δ
    @test F.extended[] === I
    # Float64 Jacobians stay compiled and do not build the interpreter.
    G = CompiledSystem(System([(p - 1) * x^2]; variables = [x], parameters = [p]))
    evaluate_and_jacobian!(zeros(ComplexF64, 1), zeros(ComplexF64, 1, 1), G, [2.0], [3.0])
    @test isnothing(G.extended[])

    H = CompiledHomotopy(Homotopy([(t - 1) * x], [x], t))
    evaluate_and_jacobian!(u, J, H, C[2], 1 + δ)
    @test u[1] == 2δ
    @test J[1] == δ
    @test H.extended[] isa InterpretedHomotopy
end
