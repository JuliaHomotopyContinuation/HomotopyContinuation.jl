@testset "utils.jl" begin
    @test HC2.fast_abs(-2.0) ≈ abs(-2.0)
    @test HC2.fast_abs(2.0 + 3im) ≈ abs(2.0 + 3im)

    @test HC2.nthroot(42.2, 0) == 1.0
    @test HC2.nthroot(42.2, 1) == 42.2
    @test HC2.nthroot(42.1, 2) == sqrt(42.1)
    @test HC2.nthroot(42.1, 3) == cbrt(42.1)
    @test HC2.nthroot(42.1, 4) == sqrt(sqrt(42.1))
    @test HC2.nthroot(42.1, 5) == 42.1^(1 / 5)

    line = HC2.ComplexLineSegment(1.0, 1.0 + 1im)
    @test line isa HC2.ComplexLineSegment
    @test length(line) ≈ 1
    @test HC2.step_size(line, 0.5) ≈ 0.5im
    @test sprint(show, line) == "ComplexLineSegment(1.0 + 0.0im, 1.0 + 1.0im)"

    struct Foo_print_fieldnames
        a::Int
        b::Float64
    end

    @test sprint(HC2.print_fieldnames, Foo_print_fieldnames(2, 2.2)) ==
          "Foo_print_fieldnames:\n • a → 2\n • b → 2.2\n"
end
