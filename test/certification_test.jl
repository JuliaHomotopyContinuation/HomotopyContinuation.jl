@testset "Certification" begin
    @testset "Simple" begin
        @var x y
        # define the polynomials
        f₁ = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        f₂ = x^2 + 2x * y^2 - 2 * y^2 - 1 / 2
        result = solve([f₁, f₂])
        cert = certify(System([f₁, f₂]), result)
        @test nnonsingular(result) == 18
        @test ncertified(cert) == 18
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 4
        @test ndistinct_real_certified(cert) == 4

        S = solutions(result)
        cert = certify(System([f₁, f₂]), [S; S])
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 4
    end

    @testset "EDD" begin
        @var x y
        f = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        # define new variables u₁, u₂ and λ₁
        @var λ[1:1] u[1:2]
        # define the jacobian of F
        J = differentiate([f], [x, y])
        # J' defines the transpose of J
        C = System([[x, y] - u - J' * λ; f], parameters = u)
        u₀ = [-0.32, -0.1]
        res = solve(C, target_parameters = u₀)
        cert = certify(C, res; target_parameters = u₀)
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 36
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 8
    end
end
