@testset "TestSystems" begin
    @testset "cyclic5 validate solutions" begin
        @TP.polyvar x[1:5]
        cyclic5 = TestSystems.cyclic(x...)

        solutions = TestSystems.cyclic5Solutions()

        @test all(sol -> norm(evaluate(cyclic5, sol)) < 1e-8, solutions) == true
    end

    @testset "cyclic7 validate solutions" begin
        @TP.polyvar x[1:7]
        cyclic7 = TestSystems.cyclic(x...)

        solutions = TestSystems.cyclic7Solutions()

        @test all(sol -> norm(evaluate(cyclic7, sol)) < 1e-8, solutions) == true
    end
end
