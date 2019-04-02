@testset "PathTracker" begin

    @testset "EarlyAtInfinity" begin
        f = equations(griewank_osborne())

        tracker, startsolutions = pathtracker_startsolutions(f, seed=130793)
        S = collect(startsolutions)
        # This path has the special case that we track towards a non-singular solution
        # at infinity and the valuation doesn't stabilize fast enough
        result = track(tracker, S[3])
        @test result.return_code == :at_infinity
    end


    @testset "Affine tracking" begin
        @polyvar x a y b
        F = [x^2-a, x*y-a+b]
        p = [a, b]

        tracker = pathtracker(F, parameters=p, p₁=[1, 0], p₀=[2, 4], affine=true)

        @test affine_tracking(tracker.core_tracker) == true
        @test HC.type_of_x(tracker) == Vector{ComplexF64}
        res = track(tracker, [1, 1])
        @test res.return_code == :success
        @test isa(solution(res), Vector{ComplexF64})
        @test length(solution(res)) == 2
        @test isprojective(res) == false
        @test isaffine(res)

        x = @SVector [1.0, 1.0 + 0.0*im]
        tracker, starts = pathtracker_startsolutions(F, x;
                                parameters=p, p₁=[1, 0], p₀=[2, 4], affine=true)
        @test length(starts) == 1
        res = track(tracker, starts[1])
        @test isa(solution(res), Vector{ComplexF64})
        @test length(solution(res)) == 2

    end

    @testset "Details Level" begin
        @polyvar x a y b
        F = [x^2-a, x*y-a+b]
        p = [a, b]

        tracker = pathtracker(F, parameters=p, p₁=[1, 0], p₀=[2, 4], affine=true)

        res = track(tracker, [1, 1]; details=:minimal)
        test_show_juno(res)
        @test res.condition_jacobian === nothing

        res = track(tracker, [1, 1]; details=:default)
        @test res.condition_jacobian !== nothing

        res = track(tracker, [1, 1]; details=:extensive)
        @test res.valuation !== nothing

        @test isnonsingular(res)
    end
end
