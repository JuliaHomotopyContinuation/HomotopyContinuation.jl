@testset "Norms" begin
    x = [2.0im, 3.0 - 1im, 5.0 + 2.0im]
    y = [-2.0im, 3.0 - 1im, 5.0 + 2.0im]
    norm = HC2.InfNorm()
    @test norm isa HC2.AbstractNorm
    @test norm(x) ≈ abs(5.0 + 2.0im)
    @test norm(x, x) ≈ 0.0
    @test norm(x, y) ≈ abs(4im)

    # test overflow check
    huge_x = exp2(700) .* x
    huge_y = exp2(700) .* y
    @test norm(huge_x) ≈ exp2(700) * abs(5.0 + 2.0im)
    @test norm(huge_x, huge_x) ≈ 0.0
    @test norm(huge_x, huge_y) ≈ exp2(700) * abs(4im)

    weighted_norm = HC2.WeightedNorm(norm, x)
    @test HC2.weights(weighted_norm) == ones(3)
    weighted_norm .= 2.0
    @test HC2.weights(weighted_norm) == [2.0, 2, 2]
    weighted_norm[1] = 4.0
    @test weighted_norm isa HC2.AbstractNorm
    @test weighted_norm(x) ≈ 0.5 * abs(5 + 2.0im)
    @test weighted_norm(x, x) ≈ 0.0
    @test weighted_norm(x, y) ≈ 0.25 * abs(4im)
    @test !isempty(sprint(show, weighted_norm))


    weighted_norm = HC2.WeightedNorm(norm, huge_x)
    weighted_norm .= 2.0
    weighted_norm[1] = 4.0
    @test weighted_norm isa HC2.AbstractNorm
    @test weighted_norm(huge_x) ≈ exp2(700) * 0.5 * abs(5 + 2.0im)
    @test weighted_norm(huge_x, huge_x) ≈ exp2(700) * 0.0
    @test weighted_norm(huge_x, huge_y) ≈ exp2(700) * 0.25 * abs(4im)
end
