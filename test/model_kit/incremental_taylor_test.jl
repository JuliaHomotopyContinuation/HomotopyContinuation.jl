@testset "Incremental predictor Taylor evaluation" begin
    @var x y p t
    F = System(
        [(x * y + p)^3, x^3 - y^2 + p * x, inv(x + p)];
        variables = [x, y],
        parameters = [p],
    )
    for T in (ComplexF64, ComplexDF64)
        G = InterpretedSystem(F)
        u, ref = zeros(T, 3), zeros(T, 3)
        for base in (1, 2, 1)
            tx = TaylorVector{4}(T.([base 2; 0.3 0.7; 0.2 0.1; 0 0]))
            tp = TaylorVector{4}(reshape(T[1, 0.2, 0.1, 0.3], 4, 1))
            for k = 1:3
                partial =
                    k == 1 ? T[base, 2] : k == 2 ? TaylorVector{2}(tx) : TaylorVector{3}(tx)
                ModelKit.predictor_taylor!(u, Val(k), G, partial, tp, k > 1)
                taylor!(ref, Val(k), G, partial, tp)
                @test u ≈ ref rtol = 1e-13
            end
            # Missing the preceding order takes the independent-evaluation path.
            ModelKit.predictor_taylor!(u, Val(3), G, TaylorVector{3}(tx), tp, true)
            @test u ≈ ref rtol = 1e-13
        end
        # Parameters may be mutated in place between successive orders too.
        tx = TaylorVector{3}(T.([2 3; 0.3 0.7; 0.2 0.1]))
        tp = TaylorVector{3}(reshape(T[1, 0.2, 0.1], 3, 1))
        for coefficient in (1, 2)
            ModelKit.predictor_taylor!(u, Val(1), G, T[2, 3], tp, false)
            ModelKit.predictor_taylor!(u, Val(2), G, TaylorVector{2}(tx), tp, true)
            tp[1, coefficient] += 1
            ModelKit.predictor_taylor!(u, Val(3), G, tx, tp, true)
            taylor!(ref, Val(3), G, tx, tp)
            @test u ≈ ref rtol = 1e-13
        end
    end
    H = InterpretedHomotopy(Homotopy([(x + t)^3], [x], t))
    u, ref = zeros(ComplexF64, 1), zeros(ComplexF64, 1)
    tx = TaylorVector{3}(reshape(ComplexF64[2, 3, 4], 3, 1))
    for k = 1:3
        partial = k == 1 ? [2.0] : k == 2 ? TaylorVector{2}(tx) : tx
        taylor!(u, Val(k), H, partial, 0.4, k > 1)
        taylor!(ref, Val(k), H, partial, 0.4)
        @test u ≈ ref
    end
    # Reuse is safe only while the coefficients already consumed are unchanged.
    for change in (:point, :time, :coefficient)
        taylor!(u, Val(1), H, [2.0], 0.4, false)
        taylor!(u, Val(2), H, TaylorVector{2}(tx), 0.4, true)
        changed = deepcopy(tx)
        changed[1, change == :coefficient ? 2 : 1] += change == :time ? 0 : 1
        time = change == :time ? 0.7 : 0.4
        taylor!(u, Val(3), H, changed, time, true)
        taylor!(ref, Val(3), H, changed, time)
        @test u ≈ ref
    end
end
