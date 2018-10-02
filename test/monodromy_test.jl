function toric_ed(A)
    d, n = size(A)
    @polyvar t[1:d] y[1:n] u[1:n]

    φ = map(j -> prod(i -> t[i]^A[i,j], 1:d), 1:n)
    Dφ = [differentiate(φ[i], t[j]) for i=1:n, j=1:d]

    F = [φ + y - u; Dφ' * y]

    # We sample a random t, find a normal vector on this point and then assemble
    # an u where we know the solution to.
    t_rand = randn(Float64, size(A, 1))
    null = nullspace(map(fi -> fi(t => t_rand), Dφ'))
    y_rand = null * randn(Float64, size(null, 2))
    u₀ = map(fi -> fi(t => t_rand), φ) + y_rand

    x₀ = complex.([t_rand; y_rand])

    F, u, u₀, x₀
end

@testset "Monodromy" begin

    @testset "monodromy_solve" begin
        F, p, p₀, x₀ = toric_ed([3 2 1 0; 0 1 2 3])

        solutions = monodromy_solve(F, p₀, x₀, parameters=p, target_solutions_count=21)
        @test length(solutions) == 21

        # test that timeout works
        Random.seed!(123)
        solutions = monodromy_solve(F, p₀, x₀, parameters=p, target_solutions_count=21, timeout=1e-9)
        @test length(solutions) < 21

        solutions = monodromy_solve(F, p₀, x₀, parameters=p, target_solutions_count=21,
                done_callback=(s -> true))
        @test length(solutions) == 2

        solutions = monodromy_solve(F, p₀, x₀, parameters=p, target_solutions_count=21,
            group_action=(s -> begin
                t = cis(π*2/3)
                t² = t * t
                (vcat(t * s[1], t * s[2], s[3:end]),
                 vcat(t² * s[1], t² * s[2], s[3:end]))
            end))
        @test length(solutions) == 21
    end
end
