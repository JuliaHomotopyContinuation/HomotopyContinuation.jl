module InterfaceTest

import ..Homotopies
import Base.Test: @test


function homotopy(H::Homotopies.AbstractHomotopy, x=rand(Complex{Float64}, size(H, 2)) )
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
    @test Homotopies.jacobian(H, x, t, cache) ≈ U
    @test Homotopies.jacobian(H, x, t) ≈ U

    Homotopies.evaluate_and_jacobian!(u, U, H, x, t, cache)
    (v, V) = Homotopies.evaluate_and_jacobian(H, x, t, cache)
    @test v ≈ u
    @test V ≈ U
    @test Homotopies.evaluate(H, x, t) ≈ u
    @test Homotopies.jacobian(H, x, t) ≈ U

    Homotopies.jacobian_and_dt!(U, u, H, x, t, cache)
    (V, v) = Homotopies.jacobian_and_dt(H, x, t, cache)
    @test V ≈ U
    @test v ≈ u
    @test Homotopies.jacobian(H, x, t) ≈ U
    @test Homotopies.dt(H, x, t) ≈ u
end

end
