@testset "NewtonCorrector" begin
    @testset "square" begin
        F = equations(griewank_osborne())

        prob, starts = problem_startsolutions(F; seed=12345)
        S = collect(starts)

        x₀ = S[1] .+ 1e-4
        t = 1.0
        H = HC.HomotopyWithCache(prob.homotopy, x₀, t)
        corr_cache = HC.cache(HC.NewtonCorrector(), H, x₀, t)

        JM = HC.JacobianMonitor(jacobian(H, x₀, t))
        x̄ = similar(x₀)

        R1 = HC.newton!(x̄, H, x₀, t, JM, HC.InfNorm(), corr_cache; tol=1e-3)
        @test R1.iters == 1
        @test R1.accuracy ≤ 1e-3
        @test R1.norm_Δx₀ < 1e-1
        @test isnan(R1.ω)
        @test norm(x̄ - S[1]) < 1e-6


        R2 = HC.newton!(x̄, H, x₀, t, JM, HC.InfNorm(), corr_cache; tol=1e-7)
        @test R2.return_code == HC.NEWT_CONVERGED
        @test R2.iters == 2
        @test R2.accuracy < 1e-7
        @test !isnan(R2.ω + R2.ω₀ + R2.θ + R2.θ₀)
        @test R2.norm_Δx₀ < 1.1e-4
        @test norm(x̄ - S[1]) < 1e-10


        R3 = HC.newton!(x̄, H, [1e-5, 1e-5], 0.0, JM, HC.InfNorm(), corr_cache)
        @test R3.return_code == HC.NEWT_TERMINATED
    end

    @testset "overdetermined" begin
        @polyvar x y
        F = SPSystem([x^2+y^2-1, x-y, (x-y)*(x^2+y^2-1)])
        t = rand()
        s = [√(0.5), √(0.5)]
        x₀ = s .+ 1e-4 .* rand(ComplexF64, 2)
        H = HomotopyWithCache(StraightLineHomotopy(F, F), x₀, t)

        corr_cache = HC.cache(HC.NewtonCorrector(), H, x₀, t)
        JM = HC.JacobianMonitor(jacobian(H, x₀, t))
        x̄ = similar(x₀)

        R1 = HC.newton!(x̄, H, x₀, t, JM, HC.InfNorm(), corr_cache; tol=1e-3)
        @test R1.iters == 1
        @test R1.accuracy ≤ 1e-3
        @test R1.norm_Δx₀ ≤ 1e-3


        R2 = HC.newton!(x̄, H, x₀, t, JM, HC.InfNorm(), corr_cache; tol=1e-7)
        @test R2.return_code == HC.NEWT_CONVERGED
        @test R2.iters == 2
        @test R2.accuracy < 1e-7
        @test !isnan(R2.ω + R2.ω₀ + R2.θ + R2.θ₀)
        @test R2.norm_Δx₀ < 1e-3
        @test norm(x̄ - s) < 1e-10
    end
end
