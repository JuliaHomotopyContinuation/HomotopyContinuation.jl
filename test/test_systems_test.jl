import FixedPolynomials: evaluate

@testset "TestSystems" begin
    @testset "cyclic5 validate solutions" begin
        cyclic5 = TestSystems.cyclic5()
        cyclic5 isa Vector{FixedPolynomials.Polynomial{Complex128}}

        solutions = TestSystems.cyclic5Solutions()

        @test all(sol -> norm(map(p -> p(sol), cyclic5)) < 1e-8, solutions)
    end

    @testset "cyclic7 validate solutions" begin
        cyclic7 = TestSystems.cyclic7()
        cyclic7 isa Vector{FixedPolynomials.Polynomial{Complex128}}

        solutions = TestSystems.cyclic7Solutions()

        @test all(sol -> norm(map(p -> p(sol), cyclic7)) < 1e-8, solutions)
    end
end
