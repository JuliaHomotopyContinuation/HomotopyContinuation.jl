using HomotopyContinuation, DynamicPolynomials, LinearAlgebra

# These tests I got send by Torkel Loman to test the BioDiffEq integration.
# These compare a direct solve of a specific instance with the approach of first solving
# a generic instance and then doing a paramter homotopy.
@polyvar x[1:3] p[1:10]

N = 1000
results1_direct = fill(0,8)
results1_template = fill(0,8)
results2_direct = fill(0,8)
results2_template = fill(0,8)
results3_direct = fill(0,8)
results3_template = fill(0,8)
results4_direct = fill(0,8)
results4_template = fill(0,8)

@time for i = 1:N
    pol = [ -x[1]*x[3]*p[3] - x[1]*p[2] + p[1],
            x[1]*x[2]*x[3]*p[8]*p[9] + x[1]*x[3]*p[7]*p[8]*p[9] - x[2]*x[3]*p[5]*p[6] - x[2]*p[5]*p[6]*p[10] - x[2]*x[3]*p[4] - x[2]*p[4]*p[10],
            x[2] + x[3] - 1.0 ]
    p_vals = [0.04,0.04,1.,1.,10.0,0.,0.04,35.,.1,.04]
    pol_pars = DynamicPolynomials.subs.(pol, Ref(p => p_vals))
    sol = solutions(solve(pol_pars, show_progress=false))
    results1_direct[length(sol)+1] +=1

    p_template = randn(ComplexF64, 10)
    f_template = DynamicPolynomials.subs.(pol, Ref(p => p_template))
    result_template = solutions(solve(f_template, show_progress=false))
    sol_again = solutions(solve(pol, result_template, parameters=p, p₁=p_template,
            p₀=ComplexF64.(p_vals), show_progress=false))
    results1_template[length(sol_again)+1] += 1

    pol2 = [ -x[1]^5*x[2]*p[5] + x[1]^4*x[3]*p[6]*p[8] - x[1]*x[2]*p[3]^4*p[5] + x[3]*p[3]^4*p[6]*p[8] - x[1]^5*p[7] + x[1]^4*x[3]*p[4] - x[1]*p[3]^4*p[7] + x[3]*p[3]^4*p[4] + x[1]^4*p[1] + x[1]^4*p[2] + p[1]*p[3]^4,
             -x[1]^5*x[2]*p[5] - x[1]*x[2]*p[3]^4*p[5] - x[1]^4*x[2]*p[7] + x[1]^4*x[3]*p[4] - x[2]*p[3]^4*p[7] + x[3]*p[3]^4*p[4] + x[1]^4*p[1] + x[1]^4*p[2] + p[1]*p[3]^4,
             x[1]*x[2]*p[5] - x[3]*p[6]*p[8] - x[3]*p[4] - x[3]*p[7] ]
    p2_vals = [0.005, 0.1, 2.8,  10, 100, 0.1, 0.01, 0.]
    pol2_pars = DynamicPolynomials.subs.(pol2, Ref(p => p2_vals))
    sol2 = solutions(solve(pol2_pars, show_progress=false))
    results2_direct[length(sol2)+1] +=1

    p2_template = randn(ComplexF64, 8)
    f2_template = DynamicPolynomials.subs.(pol2, Ref(p => p2_template))
    solve_templ2 = solve(f2_template, show_progress=false)
    result2_template = solutions(solve_templ2)
    sol2_again = solutions(solve(pol2, result2_template, precision=PRECISION_ADAPTIVE, parameters=p[1:8], p₁=p2_template, p₀=ComplexF64.(p2_vals), show_progress=false))

    results2_template[length(sol2_again)+1] +=1

    pol3 = pol2
    p3_vals = [0.005, 0.1, 2.8,  10, 100, 0.1, 0.01, 1.] #Last parameter is 1 and not 0. Here there should be 3 real solutions instead of 1.
    pol3_pars = DynamicPolynomials.subs.(pol3, Ref(p => p3_vals))
    sol3 = solutions(solve(pol3_pars, precision=PRECISION_ADAPTIVE, show_progress=false))
    results3_direct[length(sol3)+1] +=1

    p3_template = randn(ComplexF64, 8)
    f3_template = DynamicPolynomials.subs.(pol3, Ref(p => p3_template))
    result3_template = solutions(solve(f3_template, show_progress=false))
    sol3_again = solutions(solve(pol3, result3_template, parameters=p[1:8], p₁=p3_template, p₀=ComplexF64.(p3_vals), show_progress=false))
    results3_template[length(sol3_again)+1] +=1

    pol4 = [-x[1]^3*p[3] - x[1]*p[2]^2*p[3] + x[1]^2*p[1]]
    p4_vals = [1., 0.2, 1.];
    pol4_pars = DynamicPolynomials.subs.(pol4, Ref(p => p4_vals))
    sol4 = solutions(solve(pol4_pars, show_progress=false))
    results4_direct[length(sol4)+1] +=1

    p4_template = randn(ComplexF64, 3)
    f4_template = DynamicPolynomials.subs.(pol4, Ref(p => p4_template))
    result4_template = solutions(solve(f4_template, show_progress=false))
    sol4_again = solutions(solve(pol4, result4_template, parameters=p[1:3], p₁=p4_template, p₀=ComplexF64.(p4_vals), show_progress=false))
    results4_template[length(sol4_again)+1] +=1
end

# These numbers should all coincide
println("Results for the first system")
println(results1_direct)
println(results1_template,"\n")
println("Results for the second system")
println(results2_direct)
println(results2_template,"\n")
println("Results for the third system")
println(results3_direct)
println(results3_template,"\n")
println("Results for the fourth system")
println(results4_direct)
println(results4_template)
