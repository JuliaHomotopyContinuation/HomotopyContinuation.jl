@testset "utils.jl" begin
    Q = HC.rand_unitary_matrix(8, ComplexF64)
    @test norm(I - Q' * Q, Inf) < 1e-14

    @test HC.fast_abs(-2.0) ≈ abs(-2.0)
    @test HC.fast_abs(2.0 + 3im) ≈ abs(2.0 + 3im)

    @test HC.nthroot(42.2, 0) == 1.0
    @test HC.nthroot(42.2, 1) == 42.2
    @test HC.nthroot(42.1, 2) == sqrt(42.1)
    @test HC.nthroot(42.1, 3) == cbrt(42.1)
    @test HC.nthroot(42.1, 4) == sqrt(sqrt(42.1))
    @test HC.nthroot(42.1, 5) == 42.1^(1 / 5)

    @testset "SegmentStepper" begin
        seg = HC.SegmentStepper(0, 1)
        @test seg.start == 0
        @test seg.target == 1
        @test seg.abs_Δ == 1.0
        @test seg.forward
        @test !isempty(sprint(show, seg))

        HC.propose_step!(seg, exp2(-60))
        @test seg.Δs == exp2(-60)
        @test seg.s′ == exp2(-60)
        @test seg.t′ == exp2(-60)
        @test seg.t == 0
        @test seg.Δt == exp2(-60)
        HC.step_success!(seg)
        @test seg.s == exp2(-60)
        @test seg.t == exp2(-60)
        @test seg.t′ == exp2(-60)
        HC.propose_step!(seg, exp2(-50))
        @test seg.t == exp2(-60)
        @test seg.s′ == exp2(-60) + exp2(-50)
        @test seg.t′ == exp2(-60) + exp2(-50)

        HC.propose_step!(seg, 4)
        @test seg.s′ == 1.0
        @test seg.t′ == 1.0
        @test !HC.is_done(seg)
        HC.step_success!(seg)
        @test HC.is_done(seg)

        seg = HC.SegmentStepper(im, 0)
        @test seg.start == im
        @test seg.target == 0
        @test seg.abs_Δ == 1.0
        @test !seg.forward

        HC.propose_step!(seg, exp2(-20))
        @test seg.Δs == exp2(-20)
        @test seg.s′ == 1 - exp2(-20)
        @test seg.t′ == im * (1 - exp2(-20))
        @test seg.t == im
        @test seg.Δt == -im * exp2(-20)
        HC.step_success!(seg)
        @test seg.s == 1 - exp2(-20)
        @test seg.t == im * (1 - exp2(-20))
        @test seg.t′ == im * (1 - exp2(-20))
        HC.propose_step!(seg, exp2(-10))
        @test seg.t == im * (1 - exp2(-20))
        @test seg.s′ == 1 - (exp2(-20) + exp2(-10))
        @test seg.t′ == im * (1 - (exp2(-20) + exp2(-10)))

        HC.propose_step!(seg, 2)
        @test seg.s′ == 0.0
        @test seg.t′ == 0.0
        @test !HC.is_done(seg)
        HC.step_success!(seg)
        @test HC.is_done(seg)
    end

    struct Foo_print_fieldnames
        a::Int
        b::Float64
    end

    @test sprint(HC.print_fieldnames, Foo_print_fieldnames(2, 2.2)) ==
          "Foo_print_fieldnames:\n • a → 2\n • b → 2.2\n"
end
