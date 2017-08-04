using HomotopyContinuation
using Base.Test

@testset "TestSystems" begin
    @testset "cyclic5 validate solutions" begin
        cyclic5 = TestSystems.cyclic5()

        solutions = TestSystems.cyclic5Solutions()

        @test all(sol -> norm(evaluate(cyclic5, sol)) < 1e-8, solutions) == true
    end

    @testset "cyclic7 validate solutions" begin
        cyclic7 = TestSystems.cyclic7()

        solutions = TestSystems.cyclic7Solutions()

        @test all(sol -> norm(evaluate(cyclic7, sol)) < 1e-8, solutions) == true
    end
end
