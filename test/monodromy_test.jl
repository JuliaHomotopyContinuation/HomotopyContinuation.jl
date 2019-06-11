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

        # test that timeout works
        Random.seed!(51232)
        result = monodromy_solve(F, x₀, p₀, parameters=p, target_solutions_count=21, timeout=1e-12)
        @test length(result.solutions) < 21

        result = monodromy_solve(F, x₀, p₀, parameters=p,
                target_solutions_count=21,
                maximal_number_of_iterations_without_progress=20)
        @test result.returncode == :success
        @test length(solutions(result)) == 21
        @test length(solutions(result, onlyreal = true)) >= 1
        @test result.statistics.ntrackedpaths ≥ 21
        @test result.statistics.nparametergenerations ≥ 1
        @test length(HC.UniquePoints(result.solutions).points) == 21
        @test isempty(sprint(show, result)) == false

        # test input of length > 1
        result = monodromy_solve(F, [x₀ for _ in 1:30], p₀, parameters=p)
        @test length(solutions(result)) == 21

        # test with false input
        result = monodromy_solve(F, [rand(6) for _ in 1:10], p₀, parameters=p)
        @test result.returncode == :invalid_startvalue
        result = monodromy_solve(F, [x₀, rand(6)], p₀, parameters=p)
        @test length(solutions(result)) == 21
        result = monodromy_solve(F, [x₀, rand(6)], p₀, parameters=p, check_startsolutions = false)
        @test result.returncode == :invalid_startvalue

        # different distance function
        result = monodromy_solve(F, x₀, p₀, parameters=p, distance = (x,y) -> 0.0)
        @test length(solutions(result)) == 1

        # Test stop heuristic using too high target_solutions_count
        result = monodromy_solve(F, x₀, p₀, parameters=p, target_solutions_count=25)
        @test result.returncode == :heuristic_stop
        # Test stop heuristic with no target solutions count
        result = monodromy_solve(F, x₀, p₀, parameters=p)
        @test result.returncode == :heuristic_stop
        # Test stop heuristic with no target solutions count
        result = monodromy_solve(F, x₀, p₀, parameters=p, strategy=Triangle(useweights=true))
        @test result.returncode == :heuristic_stop


        # By group_actions=nothing we force that complex conjugation is not used.
        result2 = monodromy_solve(F, x₀, p₀, parameters=p,
                        target_solutions_count=21, complex_conjugation=false,
                        maximal_number_of_iterations_without_progress=100)
        @test result2.returncode == :success

        result = monodromy_solve(F, x₀, p₀, parameters=p, target_solutions_count=21,
                done_callback=((_, _) -> true))
        @test length(result.solutions) == 2

        roots_of_unity(s) = begin
            t = cis(π*2/3)
            t² = t * t
            (vcat(t * s[1], t * s[2], s[3:end]),
             vcat(t² * s[1], t² * s[2], s[3:end]))
        end

        result = monodromy_solve(F, x₀, p₀, parameters=p, target_solutions_count=21,
            maximal_number_of_iterations_without_progress=100,
            equivalence_classes=false,
            group_action=roots_of_unity)
        @test length(result.solutions) == 21

        result = monodromy_solve(F, x₀, p₀, parameters=p, target_solutions_count=21,
            maximal_number_of_iterations_without_progress=100,
            complex_conjugation=false, # disable complex conjugation to test it as a group action.
            equivalence_classes=false,
            group_actions=(roots_of_unity, s -> (conj.(s),)))
        @test length(result.solutions) == 21
        @test length(solutions(result)) == 21
        @test length(realsolutions(result)) < 21
        test_treeviews(result)

        # group_actions as a vector
        result = monodromy_solve(F, x₀, p₀, parameters=p, target_solutions_count=21,
            maximal_number_of_iterations_without_progress=100,
            complex_conjugation=false,
            equivalence_classes=false,
            group_actions=[roots_of_unity, s -> (conj.(s),)])
        @test length(result.solutions) == 21


        # equivalence classes
        result = monodromy_solve(F, x₀, p₀, parameters=p,
            equivalence_classes=true,
            target_solutions_count=7,
            maximal_number_of_iterations_without_progress=100,
            group_actions=roots_of_unity)
        @test length(result.solutions) == 7
        # Test that equivalence classes are on by default if we supply a group action
        result = monodromy_solve(F, x₀, p₀, parameters=p,
                            group_action=roots_of_unity,
                            maximal_number_of_iterations_without_progress=20)
        @test length(result.solutions) == 7

        # Test affine tracking
        result = monodromy_solve(F, x₀, p₀, parameters=p, affine_tracking=true,
                        group_action=roots_of_unity,
                        target_solutions_count=7,
                        maximal_number_of_iterations_without_progress=200)
        @test length(result.solutions) == 7

        # AbstractSystem as input
        F_p = SPSystem(F; parameters=p)
        result = monodromy_solve(F_p, x₀, p₀, affine_tracking=true,
                        group_action=roots_of_unity,
                        target_solutions_count=7,
                        maximal_number_of_iterations_without_progress=200)
        @test length(result.solutions) == 7
    end

    @testset "Method of Moments" begin
        Random.seed!(130793)
        @polyvar a[1:3] x[1:3] s[1:3]

        f0 = a[1]+a[2]+a[3];
        f1 = a[1]*x[1]+a[2]*x[2]+a[3]*x[3];
        f2 = a[1]*(x[1]^2+s[1])+a[2]*(x[2]^2+s[2])+a[3]*(x[3]^2+s[3]);
        f3 = a[1]*(x[1]^3+3*s[1]*x[1])+a[2]*(x[2]^3+3*s[2]*x[2])+a[3]*(x[3]^3+3*s[3]*x[3]);
        f4 = a[1]*(x[1]^4+6*s[1]*x[1]^2+3*s[1]^2)+a[2]*(x[2]^4+6*s[2]*x[2]^2+3*s[2]^2)+
        	a[3]*(x[3]^4+6*s[3]*x[3]^2+3*s[3]^2);
        f5= a[1]*(x[1]^5+10*s[1]*x[1]^3+15*x[1]*s[1]^2)+a[2]*(x[2]^5+10*s[2]*x[2]^3+15*x[2]*s[2]^2)+a[3]*(x[3]^5+10*s[3]*x[3]^3+15*x[3]*s[3]^2);
        f6= a[1]*(x[1]^6+15*s[1]*x[1]^4+45*x[1]^2*s[1]^2+15*s[1]^3)+
        	a[2]*(x[2]^6+15*s[2]*x[2]^4+45*x[2]^2*s[2]^2+15*s[2]^3)+
        	a[3]*(x[3]^6+15*s[3]*x[3]^4+45*x[3]^2*s[3]^2+15*s[3]^3);
        f7= a[1]*(x[1]^7+21*s[1]*x[1]^5+105*x[1]^3*s[1]^2+105*x[1]*s[1]^3)+
        	a[2]*(x[2]^7+21*s[2]*x[2]^5+105*x[2]^3*s[2]^2+105*x[2]*s[2]^3)+
        	a[3]*(x[3]^7+21*s[3]*x[3]^5+105*x[3]^3*s[3]^2+105*x[3]*s[3]^3);
        f8= a[1]*(x[1]^8+28*s[1]*x[1]^6+210*x[1]^4*s[1]^2+420*x[1]^2*s[1]^3+105*s[1]^4)+
        	a[2]*(x[2]^8+28*s[2]*x[2]^6+210*x[2]^4*s[2]^2+420*x[2]^2*s[2]^3+105*s[2]^4)+
        	a[3]*(x[3]^8+28*s[3]*x[3]^6+210*x[3]^4*s[3]^2+420*x[3]^2*s[3]^3+105*s[3]^4)

        f = [f0, f1, f2, f3, f4, f5, f6, f7, f8]


        S₃ = SymmetricGroup(3)
        relabeling(v) = map(p -> (v[1:3][p]..., v[4:6][p]..., v[7:9][p]...), S₃)

        @polyvar p[1:9]

        function sample_moments(_=nothing)
        	a₀ = rand(3)
        	a₀ ./= sum(a₀)
        	x₀ = randn(3)
        	s₀ = rand(3) .* 2

        	y₀ = [a₀; x₀; s₀]
        	p₀ = [fᵢ(a => a₀, x => x₀, s => s₀) for fᵢ in f]
        	y₀, p₀
        end

        y₀, p₀ = sample_moments()

        R = monodromy_solve(f - p, y₀, p₀;
        			parameters=p, group_action=relabeling,
        			parameter_sampler=last ∘ sample_moments,
                    show_progress=false,
                    maximal_number_of_iterations_without_progress=5)
        @test length(solutions(R)) ≤ 225

        R = monodromy_solve(f - p, y₀, p₀;
        			parameters=p, group_action=relabeling,
        			target_solutions_count=225,
                    affine_tracking=true,
        			parameter_sampler=last ∘ sample_moments,
                    show_progress=false,
                    maximal_number_of_iterations_without_progress=20)
        @test length(solutions(R)) == 225
    end
end
