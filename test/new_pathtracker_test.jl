using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));

using HomotopyContinuation, HomotopyContinuationGym, Test, ProjectiveVectors
using LinearAlgebra
const HC = HomotopyContinuation
import PolynomialTestSystems: equations, katsura, griewank_osborne, cyclic
import PolynomialTestSystems

@testset "PathTracker" begin
    # Simple
    @testset "Simple" begin
        @polyvar x y
        f = [x^2 - 2, x + y - 1]
        @test pathtracker(f; system = FPSystem) isa PathTracker{Vector{ComplexF64}}
        @test pathtracker(
            f;
            projective_tracking = true,
            system = FPSystem,
        ) isa PathTracker{PVector{ComplexF64,1}}

        tracker, starts = pathtracker_startsolutions(f; system = FPSystem)
        S = collect(starts)
        for i = 1:2
            @test is_success(track!(tracker, S[i]))
            s = solution(tracker)
            @test abs(s[1]) ≈ sqrt(2) atol = 1e-8
            @test sum(s) ≈ 1 atol = 1e-8
            @test isnothing(winding_number(tracker))
        end

        g = katsura(5)
        tracker, starts = pathtracker_startsolutions(f; system = FPSystem)
        @test all(s -> is_success(track!(tracker, s)), S)
    end

    # Univariate singular
    @testset "Univariate singular" begin
        @polyvar x
        for n = 2:12
            f = [(x - 3)^n]
            g = [x^n - 1]
            S = [[cis(i * 2π / n)] for i = 0:(n-1)]

            tracker = pathtracker(g, f, S; seed = 842121, system = FPSystem)
            @test all(S) do s
                is_success(track!(tracker, s)) &&
                winding_number(tracker) == n &&
                isapprox(solution(tracker)[1], 3.0; atol = 1e-6)
            end
        end
    end

    # Wilkinson
    @testset "Wilkinson" begin
        @polyvar x
        for n = 2:18
            g = [x^n - 1]
            f = [prod(x - i for i = 1:n)]
            S = [[cis(i * 2π / n)] for i = 0:(n-1)]
            tracker = pathtracker(g, f, S; seed = 512321, system = FPSystem)
            @test all(s -> is_success(track!(tracker, s)), S)
        end
    end

    # Bivariate singular
    @testset "Bivariate singular" begin
        f = equations(griewank_osborne())
        tracker, starts = pathtracker_startsolutions(f; seed = 78373, system = FPSystem)
        @test all(starts) do s
            retcode = track!(tracker, s)
            (is_success(retcode) && winding_number(tracker) == 3) || is_at_infinity(retcode)
        end

        # Bivariate projetive singular
        @polyvar x z y
        # This has two roots of multiplicity 6 at the hyperplane z=0.
        # But the winding numbers are only 3 at each singularity
        F = [
            0.75 * x^4 + 1.5 * x^2 * y^2 - 2.5 * x^2 * z^2 + 0.75 * y^4 - 2.5 * y^2 * z^2 +
            0.75 * z^4,
            10 * x^2 * z + 10 * y^2 * z - 6 * z^3,
        ]
        @test pathtracker(F) isa PathTracker{PVector{ComplexF64,1}}
        tracker, starts = pathtracker_startsolutions(F; seed = 29831, system = FPSystem)
        @test all(s -> is_success(track!(tracker, s)), starts)
        @test all(starts) do s
            retcode = track!(tracker, s)
            (is_success(retcode) && winding_number(tracker) == 3)
        end
        # affine
        F̂ = subs.(F, Ref(y => 1))
        tracker, starts = pathtracker_startsolutions(F̂; seed = 29831, system = FPSystem)
        @test all(starts) do s
            retcode = track!(tracker, s)
            (is_success(retcode) && winding_number(tracker) == 3)
        end

        # Fix z == 1 --> all solutions at infinity
        G = subs.(F, Ref(z => 1))
        tracker, starts = pathtracker_startsolutions(G; seed = 29831, system = FPSystem)
        @test all(starts) do s
            is_at_infinity(track!(tracker, s))
        end
    end

    @testset "Non-Singular bad-conditioned" begin
        f = equations(PolynomialTestSystems.bacillus_subtilis())
        tracker, starts = pathtracker_startsolutions(f; seed = 78373, system = FPSystem)
        @test count(s -> is_success(track!(tracker, s)), starts) == 44
    end
end
