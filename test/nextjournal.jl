@testset "nextjournal" begin
    @testset "A first example" begin
        @polyvar x y;

        # define f
        f = [x^2 + 2y, y^2 - 2]

        # solve f
        result = solve(f)
        @test nfinite(result) == 4
    end
	@testset "Start systems" begin
        @polyvar x y

        f = [x^2 + 2y, y^2 - 2]
        g = [x^2 - 1, y^2 - 1]

        start_sols = [
          [1,1], [1,-1], [-1,1], [-1,-1]
        ]

        result = solve(g, f, start_sols)
        @test nfinite(result) == 4
        @test length(real_solutions(result)) == 2
    end
    @testset "Groups of variables" begin
        @polyvar x y z
        J = 3x^3*y+y^2*z^2-2x*y-x*4z^3
        g = x^4+y^4+z^4-1

        # Introduce auxillary variable for Lagrangian
        @polyvar λ

        # define Lagrangian
        L = J - λ * g;
        import DynamicPolynomials: differentiate
        ∇L = differentiate(L, [x, y, z, λ])

        # Now we solve the polynomial system ∇L = 0
        result = solve(∇L; show_progress=false)
        @test nfinite(result) == 108

        result = solve(∇L, variable_groups = [(x,y,z), (λ,)], show_progress = false)
        @test nfinite(result) == 108

        reals = real_solutions(result)
        minval, minindex = findmin(map(s -> J(s[1:3]), reals))
        minarg = reals[minindex][1:3]
        @test minarg ≈ [0.6893448348668392; 0.19072105130305433; 0.9376180557378104] ||
		      minarg ≈ -[0.6893448348668392; 0.19072105130305433; 0.9376180557378104]
    end
    @testset "6R inverse problem" begin
        # initialize the variables
        @polyvar z[1:6,1:3] p[1:3]
        α = randn(5)
        a = randn(9)

        # define the system of polynomials
        f = [z[i,:] ⋅ z[i,:] for i = 2:5]
        g = [z[i,:] ⋅ z[i+1,:] for i = 1:5]
        h = sum(a[i] .* (z[i,:] × z[i+1,:]) for i=1:3) +
            sum(a[i+4] .* z[i,:] for i = 2:5)
        F′ = [f .- 1; g .- cos.(α); h .- p]

        # assign values to z₁ and z₆
        z₁ = normalize!(randn(3))
        z₆ = normalize!(randn(3))
        F = [subs(f, z[1,:]=>z₁, z[6,:]=>z₆, p=>[1, 1, 0]) for f in F′];

        result = solve(F, show_progress = false)
        @test nfinite(result) == 16

        variable_groups=[[z[2,:]; z[4,:]], [z[3,:]; z[5,:]]];
        result = solve(F; variable_groups=variable_groups, show_progress = false)
        @test nfinite(result) == 16

        p_rand = randn(ComplexF64, 3)
        F_rand = [subs(f, z[1,:]=>z₁, z[6,:]=>z₆, p=>p_rand) for f in F′]
        R_rand = solve(F_rand, variable_groups=variable_groups, show_progress=false)

        F̂ = [subs(f, z[1,:]=>z₁, z[6,:]=>z₆) for f in F′]
        q = [2,3,4]
        result = solve(F̂, solutions(R_rand); parameters=p, start_parameters=p_rand, target_parameters=q)
        @test nfinite(result) == 16
    end
    @testset "Monodromy Method" begin
        @polyvar a[1:3] μ[1:3] σ²[1:3]

		m0 = a[1]+a[2]+a[3];
		m1 = a[1]*μ[1]+a[2]*μ[2]+a[3]*μ[3];
		m2 = a[1]*(μ[1]^2+σ²[1])+a[2]*(μ[2]^2+σ²[2])+a[3]*(μ[3]^2+σ²[3]);
		m3 = a[1]*(μ[1]^3+3*σ²[1]*μ[1])+a[2]*(μ[2]^3+3*σ²[2]*μ[2])+a[3]*
					(μ[3]^3+3*σ²[3]*μ[3]);
		m4 = a[1]*(μ[1]^4+6*σ²[1]*μ[1]^2+3*σ²[1]^2)+a[2]*
					(μ[2]^4+6*σ²[2]*μ[2]^2+3*σ²[2]^2)+
		      a[3]*(μ[3]^4+6*σ²[3]*μ[3]^2+3*σ²[3]^2);
		m5 = a[1]*(μ[1]^5+10*σ²[1]*μ[1]^3+15*μ[1]*σ²[1]^2)+
		      a[2]*(μ[2]^5+10*σ²[2]*μ[2]^3+15*μ[2]*σ²[2]^2)+
		      a[3]*(μ[3]^5+10*σ²[3]*μ[3]^3+15*μ[3]*σ²[3]^2);
		m6 = a[1]*(μ[1]^6+15*σ²[1]*μ[1]^4+45*μ[1]^2*σ²[1]^2+15*σ²[1]^3)+
		      a[2]*(μ[2]^6+15*σ²[2]*μ[2]^4+45*μ[2]^2*σ²[2]^2+15*σ²[2]^3)+
		      a[3]*(μ[3]^6+15*σ²[3]*μ[3]^4+45*μ[3]^2*σ²[3]^2+15*σ²[3]^3);
		m7 = a[1]*(μ[1]^7+21*σ²[1]*μ[1]^5+105*μ[1]^3*σ²[1]^2+105*μ[1]*σ²[1]^3)+
		      a[2]*(μ[2]^7+21*σ²[2]*μ[2]^5+105*μ[2]^3*σ²[2]^2+105*μ[2]*σ²[2]^3)+
		      a[3]*(μ[3]^7+21*σ²[3]*μ[3]^5+105*μ[3]^3*σ²[3]^2+105*μ[3]*σ²[3]^3);
		m8 = a[1] * (μ[1]^8+28*σ²[1]*μ[1]^6+210*μ[1]^4*σ²[1]^2+420*μ[1]^2*σ²[1]^3+105*σ²[1]^4)+
		      a[2] * (μ[2]^8+28*σ²[2]*μ[2]^6+210*μ[2]^4*σ²[2]^2+420*μ[2]^2*σ²[2]^3+105*σ²[2]^4)+
		      a[3] * (μ[3]^8+28*σ²[3]*μ[3]^6+210*μ[3]^4*σ²[3]^2+420*μ[3]^2*σ²[3]^3+105*σ²[3]^4)

		f = [m0, m1, m2, m3, m4, m5, m6, m7, m8]

		p₀ =  [1; -1; 3; -5.5; 22.45; -50.75; 243.325; -635.725; 3420.7375]

		@test bezout_number(f - p₀) == 362880

		a₀ = randn(ComplexF64, 3); a₀ = a₀ / sum(a₀)
		μ₀ = rand(ComplexF64, 3)
		σ²₀ = rand(ComplexF64, 3)
		start_sol = [a₀; μ₀; σ²₀]

		q₀ = [m([a; μ; σ²] => start_sol) for m in f]
		@polyvar p[1:9]
		R = monodromy_solve(f - p, start_sol, q₀; parameters=p, target_solutions_count = 1350, show_progress=false)
		R2 = solve(f - p, solutions(R); parameters = p,
  				 start_parameters=q₀, target_parameters = p₀,
				 precision=PRECISION_ADAPTIVE,
				 show_progress = false)

		all_real_sols = real_solutions(R2)
		true_real_solutions  = filter(s -> all(s[7:9] .> -1e-3), all_real_sols);

		S₃ = SymmetricGroup(3)
		relabeling = GroupActions(v -> map(p -> (v[1:3][p]..., v[4:6][p]..., v[7:9][p]...), S₃))
		mults = multiplicities(true_real_solutions, group_action = relabeling)
		@test length(mults) ≥ 1

		R_with_group_action = monodromy_solve(f - p, start_sol, q₀;
					parameters=p, group_action = relabeling,
					show_progress=false,
					target_solutions_count=225)

		@test all(length.([relabeling(s) for s in solutions(R_with_group_action)]) .== 6)
    end
end
