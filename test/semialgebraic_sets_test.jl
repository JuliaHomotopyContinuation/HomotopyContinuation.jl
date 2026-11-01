using Test, HomotopyContinuation
import SemialgebraicSets, MixedSubdivisions
import MultivariatePolynomials as MP
import MultivariateBases as MB

@testset "MultivariatePolynomials conversion" begin
    @polyvar x y a
    for T in (Expression, Variable, HomotopyContinuation.ModelKit.ExpressionRef)
        @test float(T) === Expression
    end
    F = System([x^2 - y, y^2 - a]; variables = [x, y], parameters = [a])
    @test F([2, 4], [16]) == [0, 0]

    polynomials = [zero(x + y), x - y]
    F = System(polynomials; variables = [y, x])
    @test F([4, 2]) == [0, -2]
    G, c = HomotopyContinuation.ModelKit.system_with_coefficents_as_params(
        polynomials;
        variables = [y, x],
    )
    @test G([4, 2], c) == [0, -2]
    @test length(c) == 2

    F = System([x + y]; variable_groups = [[x], [y]])
    G, c = HomotopyContinuation.ModelKit.system_with_coefficents_as_params(
        [x + y];
        variable_groups = [[x], [y]],
    )
    for H in (F, G)
        @test string.(only.(HomotopyContinuation.ModelKit.variable_groups(H))) == ["x", "y"]
    end
    @test F([2, 3]) == [5]
    @test G([2, 3], c) == [5]

    # MP 0.6 also represents polynomials in bases other than monomials.
    if isdefined(MP, :BasisElement)
        @testset "$B" for (B, p, values, expected) in (
            (MB.Chebyshev, 2.0x^2 - 1, [2], 7),
            (MB.ScaledMonomial, 1.0x * y, [2, 3], 6),
        )
            q = MB.convert_basis(MB.FullBasis{B}(MP.variables(p)), p)
            @test System([q])(values) ≈ [expected]
            G, c = HomotopyContinuation.ModelKit.system_with_coefficents_as_params([q])
            @test G(values, c) ≈ [expected]
            @test c == MP.coefficients(q)
        end
    end

    polynomials = [x^2 - 1, y^2 - 1]
    @test MixedSubdivisions.mixed_volume(polynomials; show_progress = false) == 4
    result = solve(polynomials; compile = false, show_progress = false, seed = UInt32(1234))
    @test length(real_solutions(result)) == 4
end

@testset "SemialgebraicSets" begin
    solver = SemialgebraicSetsHCSolver(; compile = false)
    @polyvar x y
    V = SemialgebraicSets.@set x^2 == 1 && y^2 == 2 solver
    S = collect(V)
    S = sort!(map(s -> round.(s, digits = 2), S))
    @test S == [[-1.0, -1.41], [-1.0, 1.41], [1.0, -1.41], [1.0, 1.41]]

    # Inspired from https://jump.dev/SumOfSquares.jl/v0.4.6/generated/Polynomial%20Optimization/
    f =
        (ε, atol) -> begin
            @testset "ε=$ε atol=$atol" begin
                solver = SemialgebraicSetsHCSolver(;
                    excess_residual_tol = atol,
                    real_atol = atol,
                    compile = false,
                )
                @test sprint(show, solver) ==
                      "SemialgebraicSetsHCSolver(; excess_residual_tol = $atol, real_atol = $atol, real_rtol = 0.0, compile = false)"
                o = 1 + ε
                V = SemialgebraicSets.algebraicset(
                    [-x - y + o, -o * x * y + o * y^2 - y, -o * x^2 + y^2 - 2y + o],
                    solver,
                )
                S = collect(V)
                S = sort!(map(s -> round.(s, digits = 2), S))
                @test length(S) == 2
                @test S[1] ≈ [0.0, 1.0] atol = atol
                @test S[2] ≈ [1.0, 0.0] atol = atol
            end
        end
    f(1e-5, 1e-4)
    f(1e-4, 1e-3)
    f(1e-3, 1e-2)
    f(1e-2, 1e-1)
end
