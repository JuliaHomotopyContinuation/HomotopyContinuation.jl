function testevaluations(H, x)
    m, n = size(H)
    t = rand()
    u = zeros(Complex{Float64}, m)
    U = zeros(Complex{Float64}, m, n)

    cache = Homotopies.cache(H, x, t)

    Homotopies.evaluate!(u, H, x, t, cache)
    @test Homotopies.evaluate(H, x, t, cache) == u
    @test Homotopies.evaluate(H, x, t) == u

    Homotopies.dt!(u, H, x, t, cache)
    @test Homotopies.dt(H, x, t, cache) == u
    @test Homotopies.dt(H, x, t) == u

    Homotopies.jacobian!(U, H, x, t, cache)
    @test Homotopies.jacobian(H, x, t, cache) == U
    @test Homotopies.jacobian(H, x, t) == U

    Homotopies.evaluate_and_jacobian!(u, U, H, x, t, cache)
    @test Homotopies.evaluate_and_jacobian(H, x, t, cache) == (u, U)
    @test (u, U) == (Homotopies.evaluate(H, x, t), Homotopies.jacobian(H, x, t))

    Homotopies.jacobian_and_dt!(U, u, H, x, t, cache)
    @test Homotopies.jacobian_and_dt(H, x, t, cache) == (U, u)
    @test (U, u) == (Homotopies.jacobian(H, x, t), Homotopies.dt(H, x, t))
end

@testset "StraightLineHomotopy" begin
    F = Systems.SPSystem(TestSystems.equations(katsura5()))
    G = Systems.SPSystem(TestSystems.equations(cyclic6()))
    H = Homotopies.StraightLineHomotopy(F, G)
    @test H isa Homotopies.AbstractHomotopy
    @test size(H) == (6, 6)


    testevaluations(H, rand(Complex{Float64}, 6))
end

@testset "FixedPointHomotopy" begin
    F = Systems.SPSystem(TestSystems.equations(katsura5()))
    H = Homotopies.FixedPointHomotopy(F, rand(Complex128, 6))
    @test H isa Homotopies.AbstractHomotopy
    @test size(H) == (6, 6)

    testevaluations(H, rand(Complex{Float64}, 6))
end

@testset "PatchedHomotopy" begin
    F = Systems.SPSystem(TestSystems.equations(katsura5()))
    G = Systems.SPSystem(TestSystems.equations(cyclic6()))
    x = ProjectiveVectors.PVector(rand(Complex{Float64}, 6), 1)
    H = Homotopies.PatchedHomotopy(Homotopies.StraightLineHomotopy(F, G),
        AffinePatches.OrthogonalPatch(),
        x)
    @test H isa Homotopies.AbstractHomotopy
    @test size(H) == (7, 6)

    testevaluations(H, x)
end
