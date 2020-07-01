@testset "Certification" begin
    @testset "Simple 1: Input = Vector{Expression}" begin
        @var x y
        # define the polynomials
        f₁ = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        f₂ = x^2 + 2x * y^2 - 2 * y^2 - 1 / 2
        F = [f₁, f₂]
        result = solve(F)

        # Input: System(F)
        cert = certify(System(F), result)
        @test nnonsingular(result) == 18
        @test ncertified(cert) == 18
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 4
        @test ndistinct_real_certified(cert) == 4
        save("tmp_cert.txt", cert)
        @test !isempty(read("tmp_cert.txt", String))

        # Input: F
        cert = certify(F, result)
        @test nnonsingular(result) == 18
        @test ncertified(cert) == 18
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 4
        @test ndistinct_real_certified(cert) == 4

        # Control Display
        cert = certify(F, result, show_progress = false)

        # Double solutions
        S = solutions(result)
        cert = certify(F, [S; S])
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 4
    end

    @testset "Simple 2: Input = Vector{MP.Polynomial}" begin
        @polyvar x y
        # define the polynomials
        f₁ = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        f₂ = x^2 + 2x * y^2 - 2 * y^2 - 1 / 2
        F = [f₁, f₂]
        result = solve(F)

        # Input: System(F)
        cert = certify(System(F), result)
        @test nnonsingular(result) == 18
        @test ncertified(cert) == 18
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 4
        @test ndistinct_real_certified(cert) == 4

        # Input: F
        cert = certify(F, result)
        @test nnonsingular(result) == 18
        @test ncertified(cert) == 18
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 4
        @test ndistinct_real_certified(cert) == 4

        # Control Display
        cert = certify(F, result, show_progress = false)

        # Double solutions
        S = solutions(result)
        cert = certify(System(F), [S; S])
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 4
    end

    @testset "Parameters: Input = Vector{Expression}" begin
        @var x y
        f = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        # define new variables u₁, u₂ and λ₁
        @var λ[1:1] u[1:2]
        # define the jacobian of F
        J = differentiate([f], [x, y])
        # J' defines the transpose of J
        F = [[x, y] - u - J' * λ; f]
        u₀ = [-0.32, -0.1]
        res = solve(F, parameters = u, target_parameters = u₀)

        # Input: System(F)
        C = System(F, parameters = u)
        cert = certify(C, res; target_parameters = u₀)
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 36
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 8

        # Input: F
        cert = certify(F, res, u₀; parameters = u)
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 36
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 8

        # Invalid solutions
        cert = certify(C, [randn(ComplexF64, 3) for i in 1:10], u₀)
        @test ncertified(cert) == 0
    end

    @testset "Parameters: Input = Vector{MP.Polynomial}" begin
        @polyvar x y
        f = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        # define new variables u₁, u₂ and λ₁
        @polyvar λ[1:1] u[1:2]
        # define the jacobian of F
        J = differentiate([f], [x, y])
        # J' defines the transpose of J
        F = [[x, y] - u - J' * λ; f]
        u₀ = [-0.32, -0.1]
        res = solve(F, parameters = u, target_parameters = u₀)

        # Input: System(F)
        C = System(F, parameters = u)
        cert = certify(C, res; target_parameters = u₀)
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 36
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 8

        # Input: F
        cert = certify(F, res; parameters = u, target_parameters = u₀)
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 36
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 8
    end

    @testset "positive" begin
        @var x y
        f = System([x^2 + y^2 - 1, x - y])
        res = solve(f; compile = false, start_system=:total_degree)
        cert = certify(f, res, compile = false)
        @test count(is_positive, certificates(cert)) == 1
        @test count(is_real, certificates(cert)) == 2
    end
end
