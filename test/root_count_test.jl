@testset "Root counts" begin

    @testset "Benchmark systems" begin
        kat10 = solve(equations(katsura(10)); show_progress = false)
        @test nfinite(kat10) == 1024
        @test nsingular(kat10) == 0
        @test nreal(kat10) == 216

        cyclic5 = solve(equations(cyclic(5)); show_progress = false)
        @test nfinite(cyclic5) == 70
        @test nreal(cyclic5) == 10
        @test nsingular(cyclic5) == 0

        cyclic7 = solve(equations(cyclic(7)); show_progress = false)
        @test nfinite(cyclic7) == 924
        @test nreal(cyclic7) == 56
        @test nsingular(cyclic7) == 0
        @test ntracked(cyclic7) == 5040
    end

    @testset "Multi-homogenous" begin
        ## initialize the variables
        @polyvar z[1:6, 1:3]
        F = let
            p = [1, 1, 0]
            α = [
                0.9606655201575953,
                -0.12634585123792536,
                -0.821994270849398,
                -0.4612529209681106,
                0.1914597581898848,
            ]
            a = [
                -0.6508165773155992,
                -0.6252070694344862,
                1.2753203087262615,
                -0.30371050526688304,
                0.45398822665915634,
                -0.4042256983795828,
                -1.1485430048515646,
                -0.22145594990414993,
                0.9590486265536351,
            ]
            ## define the system of polynomials
            f = [z[i, :] ⋅ z[i, :] for i = 2:5]
            g = [z[i, :] ⋅ z[i+1, :] for i = 1:5]
            h = sum(a[i] .* (z[i, :] × z[i+1, :]) for i = 1:3) +
                sum(a[i+4] .* z[i, :] for i = 2:5)
            F′ = [f .- 1; g .- cos.(α); h .- p]
            ## assign values to z₁ and z₆
            [subs(f, z[1, :] => [1, 0, 0], z[6, :] => [1, 0, 0]) for f in F′]
        end
        group1 = [[z[2, :]; z[4, :]], [z[3, :]; z[5, :]]]
        @test nnonsingular(solve(
            F;
            system = FPSystem,
            seed = 141232,
            show_progress = false,
        )) == 16
        @test nnonsingular(solve(
            F;
            variable_groups = group1,
            system = FPSystem,
            seed = 12312,
            show_progress = false,
        )) == 16
    end

    @testset "Bio-chemical reaction networks" begin
        @testset "bacillus_subtilis" begin
            f = equations(PolynomialTestSystems.bacillus_subtilis())
            res = solve(f; show_progress = false, seed = 173221)
            @test nsolutions(res) == 44
            @test nsingular(res) == 0
            res = solve(f; seed = 231231, start_system = :polyhedral, show_progress = false)
            @test nsolutions(res) == 44
        end

        # These test got prepared by Torkel Loman to test the BioDiffEq integration
        @polyvar x[1:3] p[1:10]
        results1_direct = fill(0, 12)
        results1_template = fill(0, 12)
        results2_direct = fill(0, 12)
        results2_template = fill(0, 12)
        results3_direct = fill(0, 12)
        results3_template = fill(0, 12)
        results4_direct = fill(0, 12)
        results4_template = fill(0, 12)

        for i = 1:100
            pol = [
                -x[1] * x[3] * p[3] - x[1] * p[2] + p[1],
                x[1] * x[2] * x[3] * p[8] * p[9] + x[1] * x[3] * p[7] * p[8] * p[9] -
                x[2] * x[3] * p[5] * p[6] - x[2] * p[5] * p[6] * p[10] -
                x[2] * x[3] * p[4] - x[2] * p[4] * p[10],
                -x[1] * x[2] * x[3] * p[8] * p[9] - x[1] * x[3] * p[7] * p[8] * p[9] +
                x[2] * x[3] * p[5] * p[6] + x[2] * p[5] * p[6] * p[10] +
                x[2] * x[3] * p[4] + x[2] * p[4] * p[10],
                x[2] + x[3] - 1.0,
            ]
            p_vals = [0.04, 0.04, 1.0, 1.0, 10.0, 0.0, 0.04, 35.0, .1, .04]
            pol_pars = DynamicPolynomials.subs.(pol, Ref(p => p_vals))
            sol = solutions(HomotopyContinuation.solve(pol_pars, show_progress = false))
            results1_direct[length(sol)+1] += 1

            p_template = randn(ComplexF64, 10)
            f_template = DynamicPolynomials.subs.(pol, Ref(p => p_template))
            result_template = solutions(HomotopyContinuation.solve(
                f_template,
                show_progress = false,
            ))
            (length(result_template) > 0) && (sol_again = solutions(HomotopyContinuation.solve(
                pol,
                result_template,
                parameters = p,
                p₁ = p_template,
                p₀ = ComplexF64.(p_vals),
                show_progress = false,
            )))
            (length(result_template) > 0) ? (results1_template[length(sol_again)+1] += 1) :
            (results1_template[1] += 1)

            pol2 = [
                -x[1]^5 * x[2] * p[5] + x[1]^4 * x[3] * p[6] * p[8] -
                x[1] * x[2] * p[3]^4 * p[5] + x[3] * p[3]^4 * p[6] * p[8] - x[1]^5 * p[7] +
                x[1]^4 * x[3] * p[4] - x[1] * p[3]^4 * p[7] +
                x[3] * p[3]^4 * p[4] + x[1]^4 * p[1] + x[1]^4 * p[2] + p[1] * p[3]^4,
                -x[1]^5 * x[2] * p[5] - x[1] * x[2] * p[3]^4 * p[5] - x[1]^4 * x[2] * p[7] +
                x[1]^4 * x[3] * p[4] - x[2] * p[3]^4 * p[7] +
                x[3] * p[3]^4 * p[4] + x[1]^4 * p[1] + x[1]^4 * p[2] + p[1] * p[3]^4,
                x[1] * x[2] * p[5] - x[3] * p[6] * p[8] - x[3] * p[4] - x[3] * p[7],
            ]
            p2_vals = [0.005, 0.1, 2.8, 10, 100, 0.1, 0.01, 0.0]
            pol2_pars = DynamicPolynomials.subs.(pol2, Ref(p => p2_vals))
            sol2 = solutions(HomotopyContinuation.solve(pol2_pars, show_progress = false))
            results2_direct[length(sol2)+1] += 1

            p2_template = randn(ComplexF64, 8)
            f2_template = DynamicPolynomials.subs.(pol2, Ref(p => p2_template))
            solve_templ2 = HomotopyContinuation.solve(f2_template, show_progress = false)
            result2_template = solutions(solve_templ2)
            sol2_again = solutions(HomotopyContinuation.solve(
                pol2,
                result2_template,
                parameters = p[1:8],
                p₁ = p2_template,
                p₀ = ComplexF64.(p2_vals),
                precision_strategy = :adaptive_always,
                max_steps = 10_000,
                show_progress = false,
            ))
            results2_template[length(sol2_again)+1] += 1

            pol3 = pol2
            p3_vals = [0.005, 0.1, 2.8, 10, 100, 0.1, 0.01, 1.0] #Last parameter is 1 and not 0. Here there should be 3 real solutions instead of 1.
            pol3_pars = DynamicPolynomials.subs.(pol3, Ref(p => p3_vals))
            sol3 = solutions(HomotopyContinuation.solve(pol3_pars, show_progress = false))
            results3_direct[length(sol3)+1] += 1

            p3_template = randn(ComplexF64, 8)
            f3_template = DynamicPolynomials.subs.(pol3, Ref(p => p3_template))
            result3_template = solutions(HomotopyContinuation.solve(
                f3_template,
                show_progress = false,
            ))
            sol3_again = solutions(HomotopyContinuation.solve(
                pol3,
                result3_template,
                parameters = p[1:8],
                p₁ = p3_template,
                p₀ = ComplexF64.(p3_vals),
                precision_strategy = :adaptive_always,
                max_steps = 10_000,
                show_progress = false,
            ))
            results3_template[length(sol3_again)+1] += 1


            pol4 = [-x[1]^3 * p[3] - x[1] * p[2]^2 * p[3] + x[1]^2 * p[1]]
            p4_vals = [1.0, 0.2, 1.0]
            pol4_pars = DynamicPolynomials.subs.(pol4, Ref(p => p4_vals))
            sol4 = solutions(HomotopyContinuation.solve(pol4_pars, show_progress = false))
            results4_direct[length(sol4)+1] += 1

            p4_template = randn(ComplexF64, 3)
            f4_template = DynamicPolynomials.subs.(pol4, Ref(p => p4_template))
            result4_template = solutions(HomotopyContinuation.solve(
                f4_template,
                show_progress = false,
            ))
            sol4_again = solutions(HomotopyContinuation.solve(
                pol4,
                result4_template,
                parameters = p[1:3],
                p₁ = p4_template,
                p₀ = ComplexF64.(p4_vals),
                precision_strategy = :adaptive_always,
                max_steps = 10_000,
                show_progress = false,
            ))
            results4_template[length(sol4_again)+1] += 1
        end

        @test results1_template[4] == 100
        @test 98 <= results1_direct[4] <= 100
        @test 98 <= results2_direct[7] <= 100
        @test 98 ≤ results2_template[7] ≤ 100
        @test 98 ≤ results3_direct[7] ≤ 100
        @test 98 ≤ results3_template[7] ≤ 100
        @test 98 ≤ results4_direct[4] ≤ 100
        @test 98 ≤ results4_template[4] ≤ 100
    end
end
