@testset "StraightLineHomotopy" begin
    fs = TestSystems.equations(katsura5())

    F = Systems.SPSystem(fs)
    G = Systems.SPSystem(3.52 .* fs)

    x = rand(Complex{Float64}, 6)
    t = rand()
    u = zeros(Complex{Float64}, 6)
    U = zeros(Complex{Float64}, 6, 6)

    H = Homotopies.StraightLineHomotopy(F, G)
    @test H isa Homotopies.AbstractHomotopy
    @test size(H) == (6, 6)

    Homotopies.evaluate!(u, H, x, t)
    @test Homotopies.evaluate(H, x, t) == u

    Homotopies.dt!(u, H, x, t)
    @test Homotopies.dt(H, x, t) == u

    Homotopies.jacobian!(U, H, x, t)
    @test Homotopies.jacobian(H, x, t) == U

    Homotopies.evaluate_and_jacobian!(u, U, H, x, t)
    @test Homotopies.evaluate_and_jacobian(H, x, t) == (u, U)

    Homotopies.jacobian_and_dt!(U, u, H, x, t)
    @test Homotopies.jacobian_and_dt(H, x, t) == (U, u)

    cache = Homotopies.cache(H, x, t)
    @test cache isa Homotopies.StraightLineHomotopyCache

    Homotopies.evaluate!(u, H, x, t, cache)
    @test Homotopies.evaluate(H, x, t, cache) == u

    Homotopies.dt!(u, H, x, t, cache)
    @test Homotopies.dt(H, x, t, cache) == u

    Homotopies.jacobian!(U, H, x, t, cache)
    @test Homotopies.jacobian(H, x, t, cache) == U

    Homotopies.evaluate_and_jacobian!(u, U, H, x, t, cache)
    @test Homotopies.evaluate_and_jacobian(H, x, t, cache) == (u, U)

    Homotopies.jacobian_and_dt!(U, u, H, x, t, cache)
    @test Homotopies.jacobian_and_dt(H, x, t, cache) == (U, u)
end

@testset "FixedPointHomotopy" begin
    fs = TestSystems.equations(katsura5())

    F = Systems.SPSystem(fs)
    G = Systems.SPSystem(3.52 .* fs)

    x = rand(Complex{Float64}, 6)
    t = rand()
    u = zeros(Complex{Float64}, 6)
    U = zeros(Complex{Float64}, 6, 6)

    H = Homotopies.StraightLineHomotopy(F, G)
    @test H isa Homotopies.AbstractHomotopy
    @test size(H) == (6, 6)

    Homotopies.evaluate!(u, H, x, t)
    @test Homotopies.evaluate(H, x, t) == u

    Homotopies.dt!(u, H, x, t)
    @test Homotopies.dt(H, x, t) == u

    Homotopies.jacobian!(U, H, x, t)
    @test Homotopies.jacobian(H, x, t) == U

    Homotopies.evaluate_and_jacobian!(u, U, H, x, t)
    @test Homotopies.evaluate_and_jacobian(H, x, t) == (u, U)

    Homotopies.jacobian_and_dt!(U, u, H, x, t)
    @test Homotopies.jacobian_and_dt(H, x, t) == (U, u)

    cache = Homotopies.cache(H, x, t)
    @test cache isa Homotopies.StraightLineHomotopyCache

    Homotopies.evaluate!(u, H, x, t, cache)
    @test Homotopies.evaluate(H, x, t, cache) == u

    Homotopies.dt!(u, H, x, t, cache)
    @test Homotopies.dt(H, x, t, cache) == u

    Homotopies.jacobian!(U, H, x, t, cache)
    @test Homotopies.jacobian(H, x, t, cache) == U

    Homotopies.evaluate_and_jacobian!(u, U, H, x, t, cache)
    @test Homotopies.evaluate_and_jacobian(H, x, t, cache) == (u, U)

    Homotopies.jacobian_and_dt!(U, u, H, x, t, cache)
    @test Homotopies.jacobian_and_dt(H, x, t, cache) == (U, u)
end
