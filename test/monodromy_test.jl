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

    @testset "Group actions" begin
        f1 = s -> (s * s,);
        f2 = s-> (2s, -s, 5s);
        f3 = s -> (s + 1,);
        action1 = Monodromy.GroupActions(f1)
        action2 = Monodromy.GroupActions(f1, f2)
        action3 = Monodromy.GroupActions(f1, f2, f3)
        @test action1(3) == (9, )
        @test action2(3) == (9, 18, -9, 45)
        @test action3(3) == (9, 18, -9, 45, 10, 19, -8, 46)
    end

    @testset "monodromy_solve" begin
        F, p, p₀, x₀ = toric_ed([3 2 1 0; 0 1 2 3])

        # test that timeout works
        Random.seed!(12345)
        result = monodromy_solve(F, p₀, x₀, parameters=p, target_solutions_count=21, timeout=1e-12)
        @test length(result.solutions) < 21


        result = monodromy_solve(F, p₀, x₀, parameters=p, target_solutions_count=21)
        @test result.returncode == :success
        @test length(result.solutions) == 21
        @test result.statistics.ntrackedpaths ≥ 21
        @test result.statistics.nparametergenerations ≥ 1
        @test length(Utilities.UniquePoints(result.solutions).points) == 21

        # By group_actions=nothing we force that complex conjugation is not used.
        result2 = monodromy_solve(F, p₀, x₀, parameters=p, target_solutions_count=21, group_actions=nothing)
        @test result2.returncode == :success

        result = monodromy_solve(F, p₀, x₀, parameters=p, target_solutions_count=21,
                done_callback=(s -> true))
        @test length(result.solutions) == 2

        result = monodromy_solve(F, p₀, x₀, parameters=p, target_solutions_count=21,
            group_action=(s -> begin
                t = cis(π*2/3)
                t² = t * t
                (vcat(t * s[1], t * s[2], s[3:end]),
                 vcat(t² * s[1], t² * s[2], s[3:end]))
            end))
        @test length(result.solutions) == 21

        result = monodromy_solve(F, p₀, x₀, parameters=p, target_solutions_count=21,
            group_actions=(s -> begin
                t = cis(π*2/3)
                t² = t * t
                (vcat(t * s[1], t * s[2], s[3:end]),
                 vcat(t² * s[1], t² * s[2], s[3:end]))
            end, Monodromy.complex_conjugation))
        @test length(result.solutions) == 21
    end
end
