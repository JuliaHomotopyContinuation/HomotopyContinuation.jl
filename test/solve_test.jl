@testset "solve" begin
    @testset "Invalid input" begin
        @polyvar x y z
        @test_throws ArgumentError solve([x-2y+2, 0])
        @test_throws ArgumentError solve([x-2, y-2], [x-2, y-2,y+2], [[2, -3]])
        @test_throws ArgumentError solve([x-2z, y^2+3z, z^3+x^3], homvar=z)
        # overdetermined and abstract system
        @test_throws ArgumentError solve(FPSystem([x-2z, y^2+3z^2, z^3+x^3, z+x]))
        @test_throws ArgumentError solve(FPSystem([x-2z, y^2+3z^2, z^3+x^3, z+x]), homvar=4)

        @test_throws ArgumentError solve([x-2z, y^2+3z^2, z^3+x^3])

        @test_throws ArgumentError solve([x-2z, y^2+3z], homvar=z)
        @test_throws ArgumentError solve([x - 1], patch=OrthogonalPatch())
        # invalid kwargs
        @test_throws ArgumentError solve(equations(cyclic(5)), def=0.4, abc=23)

        # test numerical homogeneous check fails
        @polyvar x y z
        G = FPSystem([x-2z, y^2+3z])
        @test_throws ArgumentError solve(G, homvar=3)
    end

    @testset "solve" begin
        @polyvar x
        @test nfinite(solve([x - 1], threading=false, save_all_paths=true)) == 1
        F = equations(katsura(5))
        @test nfinite(solve(F, threading=false)) == 32
        @test nfinite(solve(F, system=SPSystem, threading=false)) == 32
        @test nfinite(solve(F, system=FPSystem, threading=false)) == 32

        # Simple step size
        F = equations(katsura(5))
        @test nfinite(solve(F, simple_step_size_alg=true, threading=false)) == 32

        # scaling
        F = equations(katsura(5))
        @test nfinite(solve(F, system_scaling=nothing, threading=false)) == 32
        @test nfinite(solve(F, system_scaling=:equations, threading=false)) == 32
        @test nfinite(solve(F, system_scaling=:equations_and_variables, threading=false)) == 32
        @test_throws ArgumentError solve(F, system_scaling=:lalala)

        @test nfinite(solve(F, homotopy=StraightLineHomotopy)) == 32
        result = solve(F, predictor=Euler(), homotopy=StraightLineHomotopy)
        @test nresults(result) == 32
        @test nfinite(solve(F, accuracy=1e-5)) == 32

        result = solve(F)
        @test nfinite(result) == 32
        @test string(result) isa String

        F_hom = homogenize(F)
        @test nfinite(solve(F_hom, threading=false, patch=RandomPatch())) == 32
        @test nfinite(solve(F_hom, threading=false, patch=EmbeddingPatch())) == 32
        @test nfinite(solve(F_hom, threading=false, patch=OrthogonalPatch())) == 32

        @polyvar w
        F = equations(cyclic(5))
        result = solve(homogenize(F, w), threading=false, homvar=w)
        @test nfinite(result) == 70
        @test result isa Result

        @polyvar x y z
        G = [x-2, y+3]
        F = [x+2, y-3]
        @test nfinite(solve(G, F, [[2, -3]])) == 1
        @test nfinite(solve(G, F, [[2+0.0im, -3.0+0im]])) == 1
        @test nfinite(solve(G, F, [@SVector [2, -3]])) == 1

        F = FPSystem(homogenize(equations(cyclic(5))))
        result = solve(F, homvar=6, save_all_paths=true)
        @test nfinite(result) == 70

        F = equations(katsura(5))
        result = solve(homogenize(F), threading=false)
        @test result isa Result{<:ProjectiveVectors.PVector}
        @test isprojective(first(result))
        @test nnonsingular(result) == 32


        prob, startsolutions = problem_startsolutions(TargetSystemInput(F))
        result = solve(prob.homotopy, map(s -> HC.embed(prob, s), startsolutions))
        @test result isa Result{<:Vector}
        @test nnonsingular(result) == 32

        prob, startsolutions = problem_startsolutions(TargetSystemInput(F); affine_tracking=false)
        proj_starts = HC.embed.(prob, startsolutions)
        result = solve(prob.homotopy, proj_starts, homvar=homvars(prob)[1])
        @test result isa Result{<:Vector}
        @test nnonsingular(result) == 32
        @test nfinite(result) == 32

        # Composition
        @polyvar a b c x y z u v
        e = [u + 1, v - 2]
        f = [a * b - 2, a*c- 1]
        g = [x+y, y + 3, x + 2]
        res = solve(e ∘ f ∘ g)
        @test nnonsingular(res) == 2
        @test nnonsingular(solve(e ∘ f ∘ g, system_scaling=nothing)) == 2
        @test nnonsingular(solve(e ∘ f ∘ g, system_scaling=:equations_and_variables)) == 2

        res = solve(e ∘ f ∘ g, system=SPSystem)
        @test nnonsingular(res) == 2
    end


    @testset "Random seed" begin
        R = solve(equations(katsura(5)), seed=1234)
        @test seed(R) == 1234

        @polyvar x y

        G = [x-2, y+3]
        F = [x+2, y-3]
        @test seed(solve(G, F, [[2, -3]], seed=222)) == 222
    end

    @testset "Singular solutions" begin
        @polyvar x y z
        z = 1
        F = [x^2 + 2*y^2 + 2*im*y*z, (18 + 3*im)*x*y + 7*im*y^2 - (3 - 18*im)*x*z - 14*y*z - 7*im*z^2]
        result = solve(F)
        @test nsingular(result) == 3
        @test all(r -> r.winding_number == 3, singular(result))
    end

    @testset "Path jumping" begin
        result = solve(equations(katsura(5)); accuracy=5e-4, refinement_accuracy=1e-8,
            seed=39813, threading=false, affine_tracking=true)
        @test nreal(result) == 16
    end

    @testset "Affine vs projective" begin
        @polyvar x y z
        F = solve([x-2y, y-2z])
        @test length(F[1].solution) == 3
        @test F[1].solution isa ProjectiveVectors.PVector

        G = solve([x-2, y-2])
        @test length(G[1].solution) == 2
        @test G[1].solution isa Vector
    end

    @testset "Parameter Homotopies" begin
        @polyvar x a y b
        F = [x^2-a, x*y-a+b]
        p = [a, b]
        S = solve(F, [[1.0, 1.0 + 0.0*im]], parameters=p, p₁=[1, 0], p₀=[2, 4])

        @test S[1].solution ≈ [complex(√2), -complex(√2)]
        @test nfinite(S) == 1

        @polyvar x a y b z
        F = [x^2-a*z^2, x*y-(a-b)*z^2]
        S = solve(F, [[1.0, 1.0 + 0.0*im, 1.0]];
                parameters=[a, b],
                start_parameters=[1, 0], target_parameters=[2, 4])
        @test S isa Result{<:ProjectiveVectors.PVector}
        @test ProjectiveVectors.affine_chart(solution(S[1])) ≈ [√2, -√2]
        @test nnonsingular(S) == 1

        S2 = solve(F, [[1.0, 1.0 + 0.0*im, 1.0]];
                parameters=[a, b], p₁=[1, 0], p₀=[2, 4],
                homvar=z, threading=false)
        @test S2 isa Result{<:Vector}
        @test solution(S2[1]) ≈ [√2, -√2]
        @test nfinite(S2) == 1

        γ₁, γ₀ =randn(ComplexF64, 2)
        S2 = solve(F, [[1.0, 1.0 + 0.0*im, 1.0]];
                parameters=[a, b], p₁=[1, 0], p₀=[2, 4],
                γ₁=γ₁, γ₀=γ₀, homvar=z)
        @test solution(S2[1]) ≈ [√2, -√2]
        @test nfinite(S2) == 1

        γ₁, γ₀ = randn(ComplexF64, 2)
        # test SVector input
        S2 = solve(F, [@SVector [1.0, 1.0 + 0.0*im, 1.0]];
                    parameters=[a, b], p₁=[1, 0], p₀=[2, 4],
                    γ₁=γ₁, γ₀=γ₀, homvar=z)
        @test solution(S2[1]) ≈ [√2, -√2]
        @test nfinite(S2) == 1
    end

    @testset "Keep invalid start values" begin
        @polyvar x[1:2] a[1:2]
        F = [x[1]^2-a[1], x[1]*x[2]-a[1]+a[2]]
        startsolutions = [[100, 100]]
        p₁ = [1, 1]
        p₀ = [3im, 0.5+2im]
        res = solve(F, startsolutions; parameters=a, start_parameters=p₁, target_parameters=p₀)
        @test nfailed(res) == 1
        @test res[1].return_code == :terminated_invalid_startvalue
    end

    @testset "ParameterHomotopy with Composition" begin
        @polyvar p q a b c x y z u v
        f = [a * b - 2, a*c- 1]
        g = [x + y, y + 3, x + 2]
        res = solve(f ∘ g; system=SPSystem, threading=false)
        @test nnonsingular(res) == 2
        # parameters at the end
        f2 = [a * b - q, a*c- p]
        g = [x + y, y + 3, x + 2]
        r = solve(f2 ∘ g, solutions(res);
                parameters=[p, q], p₁=[1, 2], p₀=[2, 3],
                threading=false)
        @test nnonsingular(r) == 2

        # parameters at the beginning
        f = [a * b - 2, a*c- 1]
        g2 = [x + y, y + u, x + v]
        r = solve(f ∘ g2, solutions(res);
                    parameters=[u, v], p₁=[3, 2],
                    p₀=[-2, 3], threading=false)
        @test nnonsingular(r) == 2

        # parameter in the middle
        e = [u + 1, v - 2]
        res2 = solve(e ∘ f ∘ g; system=SPSystem)
        f2 = [a * b - q, a * c - p]
        r = solve(e ∘ f2 ∘ g, solutions(res2);
                    parameters=[p, q], p₁=[1, 2],
                    p₀=[2, 3], threading=false)
        @test nnonsingular(r) == 2
    end

    @testset "Coefficient Homotopy" begin
        @polyvar x a y b
        E = [[2 1 0; 0 0 0], [1 0; 1 0]]
        start = [[1.0+0im, -3.0, 2.0], [2.0+0im, -2.0]]
        target = [randn(ComplexF64, 3), randn(ComplexF64, 2)]
        H = CoefficientHomotopy(E, start, target)
        result = solve(H, [[1, 1]])
        @test issuccess(result[1])
    end

    @testset "Overdetermined" begin
        Random.seed!(1234567)
        @polyvar x y z w
        a = [0.713865+0.345317im, 0.705182+0.455495im, 0.9815+0.922608im, 0.337617+0.508932im]

        f = [x*z-y^2, y*w-z^2, x*w-y*z]
        L₁ = [1, -1, 1, -1] ⋅ [x, y, z, w]
        L₂ = rand(ComplexF64, 4) ⋅ [x, y, z, w]
        S = solve([f; L₁], [f; L₂], [[1, 1, 1, 1]])
        @test nnonsingular(S) == 1

        f = [x*z-y^2, y-z^2, x-y*z]
        L₁ = [1, -1, 1, -1] ⋅ [x, y, z, 1]
        L₂ = rand(ComplexF64, 4) ⋅ [x, y, z, 1]
        S = solve([f; L₁], [f; L₂], [[1, 1, 1]])
        @test nfinite(S) == 1


        @polyvar x y z
        p₁ = (y - x^2) * (x^2+y^2+z^2 - 1) * (x - 0.5)
        p₂ = (z - x^3) * (x^2+y^2+z^2 - 1) * (y - 0.5)
        p₃ = (y - x^2) * (z - x^3) * (x^2+y^2+z^2 - 1) * (z - 0.5)
        L₁ = [rand(ComplexF64, 4) ⋅ [x, y, z, 1], rand(ComplexF64, 4) ⋅ [x, y, z, 1]]
        L₂ = [rand(ComplexF64, 4) ⋅ [x, y, z, 1], rand(ComplexF64, 4) ⋅ [x, y, z, 1]]

        s = first(solutions(solve([x^2+y^2+z^2 - 1; L₁])))
        tracker, starts = pathtracker_startsolutions(
                [[p₁, p₂, p₃]; L₁], [[p₁, p₂, p₃]; L₂], s;
                affine_tracking=false)
        @test issuccess(track(tracker, s))
    end

    @testset "MultiHomogeneous" begin
        @polyvar x y u v

        f = [x*y - 6, x^2 - 5]
        @test bezout_number(f, variable_groups=[(x,), (y,)]) == 2
        @test bezout_number(f, variable_groups=((x,), (y,))) == 2
        @test bezout_number(f, variable_groups=[[x], [y]]) == 2
        S = solve(f, variable_groups=[(x,), (y,)])
        @test nnonsingular(S) == 2
        S = solve(f, variable_groups=[[x], [y]])
        @test nnonsingular(S) == 2
        g = [x*y - 6u*v, x^2 - u^2]
        S = solve(g, variable_groups=[(x,u), (y,v)], homvars=(u,v))
        @test nnonsingular(S) == 2

        # The 6-R inverse problem

        ## initialize the variables
        @polyvar z[1:6,1:3]
        F = let
            p = [1, 1, 0]
            α = randn(5)
            a = randn(9)
            ## define the system of polynomials
            f = [z[i,:] ⋅ z[i,:] for i=2:5]
            g = [z[i,:] ⋅ z[i+1,:] for i=1:5]
            h = sum(a[i] .* (z[i,:] × z[i+1,:]) for i=1:3) +
                sum(a[i+4] .* z[i,:] for i=2:5)
            F′ = [f .- 1; g .- cos.(α); h .- p]
            ## assign values to z₁ and z₆
            [subs(f, z[1,:] => [1, 0, 0], z[6,:] => [1,0,0]) for f in F′]
        end
        group1 = [[z[2,:]; z[4,:]], [z[3,:]; z[5,:]]]
        group2 = [z[2,:], z[4,:], z[3,:], z[5,:]]

        @test bezout_number(F) == 1024
        @test bezout_number(F; variable_groups=group1) == 320
        @test bezout_number(F; variable_groups=group2) == 576

        @test nnonsingular(solve(F; seed=1234, show_progress=false)) == 16
        @test nnonsingular(solve(F; variable_groups=group1, seed=1234, show_progress=false)) == 16
    end

    @testset "Singular1" begin
        @polyvar x y z
        # This has two roots of multiplicity 6 at the hyperplane z=0
        F = [0.75x^4 + 1.5x^2*y^2-2.5x^2*z^2+0.75*y^4 - 2.5y^2*z^2 + 0.75z^4; 10x^2*z + 10*y^2*z-6z^3]
        res = solve(F; threading=false)
        @test nsingular(res) == 12
        @test length.(multiplicities(solutions(res))) == [6,6]

        F_affine = subs.(F, Ref(y => 1))
        res_affine = solve(F_affine; threading=false)
        @test nsingular(res_affine) == 12
        @test length.(multiplicities(solutions(res_affine))) == [6,6]
        @test_nowarn sprint(show, res_affine)

        res_affine2 = solve(F_affine; affine_tracking=false, threading=false)
        @test nsingular(res_affine2) == 12
        @test length.(multiplicities(solutions(res_affine2))) == [6,6]
    end
end
