@testset "StraightLineHomotopy" begin
    Hom = HomotopyContinuation.Homotopies
    Systems = HomotopyContinuation.Systems

    fs = TestSystems.equations(katsura5())

    F = Systems.SPSystem(fs)
    G = Systems.SPSystem(3.52 .* fs)

    x = rand(Complex{Float64}, 6)
    t = rand()
    u = zeros(Complex{Float64}, 6)
    U = zeros(Complex{Float64}, 6, 6)

    H = Hom.StraightLineHomotopy(F, G)
    @test H isa Hom.AbstractStartTargetHomotopy
    @test size(H) == (6, 6)

    Hom.evaluate!(u, H, x, t)
    @test Hom.evaluate(H, x, t) == u

    Hom.dt!(u, H, x, t)
    @test Hom.dt(H, x, t) == u

    Hom.jacobian!(U, H, x, t)
    @test Hom.jacobian(H, x, t) == U

    Hom.evaluate_and_jacobian!(u, U, H, x, t)
    @test Hom.evaluate_and_jacobian(H, x, t) == (u, U)

    Hom.jacobian_and_dt!(U, u, H, x, t)
    @test Hom.jacobian_and_dt(H, x, t) == (U, u)

    cache = Hom.cache(H, x, t)
    @test cache isa Hom.StartTargetHomotopyCache

    Hom.evaluate!(u, H, x, t, cache)
    @test Hom.evaluate(H, x, t, cache) == u

    Hom.dt!(u, H, x, t, cache)
    @test Hom.dt(H, x, t, cache) == u

    Hom.jacobian!(U, H, x, t, cache)
    @test Hom.jacobian(H, x, t, cache) == U

    Hom.evaluate_and_jacobian!(u, U, H, x, t, cache)
    @test Hom.evaluate_and_jacobian(H, x, t, cache) == (u, U)

    Hom.jacobian_and_dt!(U, u, H, x, t, cache)
    @test Hom.jacobian_and_dt(H, x, t, cache) == (U, u)
end
