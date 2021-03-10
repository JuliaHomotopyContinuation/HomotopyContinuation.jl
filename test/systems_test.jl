@testset "System Tests" begin
    @testset "Composition System" begin
        @var x y a b

        f = System([y^2 + 2x + 3, x - 1])
        g = System([x + y * a, x - 1 * b], parameters = [a, b])
        @test parameters(g ∘ f) == [a, b]
        @test nparameters(g ∘ f) == 2

        @test parameters(g ∘ g) == [a, b]
        @test nparameters(g ∘ g) == 2

        @test parameters(f ∘ g) == [a, b]
        @test nparameters(f ∘ g) == 2

        @test parameters(f ∘ f) == []
        @test nparameters(f ∘ f) == 0

        @test parameters(f ∘ g ∘ f) == [a, b]
        @test nparameters(f ∘ g ∘ f) == 2
    end
end
