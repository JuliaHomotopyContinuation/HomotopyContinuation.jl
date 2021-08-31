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
        @test ncomplex_certified(cert) == 14
        @test ndistinct_real_certified(cert) == 4
        @test ndistinct_complex_certified(cert) == 14
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
        cert = certify(C, [100 .* randn(ComplexF64, 3) for i = 1:10], u₀)
        @test ncertified(cert) < 10
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
        res = solve(f; compile = false, start_system = :total_degree)
        cert = certify(f, res, compile = false)
        @test count(is_positive, certificates(cert)) == 1
        @test count(s -> is_positive(s, 1), certificates(cert)) == 1
        @test count(is_real, certificates(cert)) == 2
    end

    @testset "3264" begin
        F = steiner()
        real_conics = [
            10124547 // 662488724,
            8554609 // 755781377,
            5860508 // 2798943247,
            -251402893 // 1016797750,
            -25443962 // 277938473,
            1 // 1,
            520811 // 1788018449,
            2183697 // 542440933,
            9030222 // 652429049,
            -12680955 // 370629407,
            -24872323 // 105706890,
            1 // 1,
            6537193 // 241535591,
            -7424602 // 363844915,
            6264373 // 1630169777,
            13097677 // 39806827,
            -29825861 // 240478169,
            1 // 1,
            13173269 // 2284890206,
            4510030 // 483147459,
            2224435 // 588965799,
            33318719 // 219393000,
            92891037 // 755709662,
            1 // 1,
            8275097 // 452566634,
            -19174153 // 408565940,
            5184916 // 172253855,
            -23713234 // 87670601,
            28246737 // 81404569,
            1 // 1,
        ]
        real_sols = read_solutions(joinpath(@__DIR__, "data/3264_real_sols.txt"))
        cert = certify(F, real_sols, real_conics; compile = true)
        @test ndistinct_real_certified(cert) == 3264
        cert = certify(F, real_sols, real_conics; compile = false)
        @test ndistinct_real_certified(cert) == 3264
        @test_throws ArgumentError certify(F, real_sols; compile = false)
        @test_throws ArgumentError certify(F, real_sols; compile = true)
    end

    @testset "certify uses complex inversion" begin
        @var x y
        F = System([y / x + x^2 - 3], parameters = [y])
        monres = monodromy_solve(F)
        cert = certify(F, monres)
        @test ndistinct_certified(cert) >= 1
    end
end
