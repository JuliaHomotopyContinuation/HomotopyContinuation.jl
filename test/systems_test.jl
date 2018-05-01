@testset "SPSystem" begin
    Systems = HomotopyContinuation.Systems

    fs = TestSystems.equations(katsura5())

    F = Systems.SPSystem(fs)

    @test size(F) == (6, 6)
    x = rand(Complex{Float64}, 6)

    u = zeros(Complex{Float64}, 6)
    Systems.evaluate!(u, F, x)
    @test Systems.evaluate(F, x) == u

    U = zeros(Complex{Float64}, 6, 6)
    Systems.jacobian!(U, F, x)
    @test Systems.jacobian(F, x) == U

    Systems.evaluate_and_jacobian!(u, U, F, x)
    @test Systems.evaluate_and_jacobian(F, x) == (u, U)


    cache = Systems.cache(F, x)
    @test cache isa Systems.NullCache
    u = zeros(Complex{Float64}, 6)
    Systems.evaluate!(u, F, x, cache)
    @test Systems.evaluate(F, x, cache) == u

    U = zeros(Complex{Float64}, 6, 6)
    Systems.jacobian!(U, F, x, cache)
    @test Systems.jacobian(F, x, cache) == U

    Systems.evaluate_and_jacobian!(u, U, F, x, cache)
    @test Systems.evaluate_and_jacobian(F, x, cache) == (u, U)
end
