@testset "Numerical Irreducible Decomposition" begin
    Random.seed!(0x8b868a92)

    @testset "Union of 1 2-dim, 2 1-dim and 8 points" begin
        @var x, y, z
        p = (x * y - x^2) + 1 - z
        q = x^4 + x^2 - y - 1
        F = [
            p * q * (x - 3) * (x - 5)
            p * q * (y - 3) * (y - 5)
            p * (z - 3) * (z - 5)
        ]

        W = regeneration(F)
        @test degree.(W) == [2, 8, 8]

        W = regeneration(F; sorted = false)
        @test degree.(W) == [2, 8, 8]

        N = decompose(W)
        @test isa(N, NumericalIrreducibleDecomposition)

        # seed
        s = 0x42c9d504
        N = nid(F; seed = s)
        @test seed(N) == s

        N = nid(F; seed = nothing)
        @test isnothing(seed(N))

        # progress
        N = nid(F; show_progress = false)
        @test isa(N, NumericalIrreducibleDecomposition)

        N = nid(F; show_monodromy_progress = true)
        @test isa(N, NumericalIrreducibleDecomposition)

        N = nid(F; warning = false)
        @test isa(N, NumericalIrreducibleDecomposition)

        # options
        N = nid(F; sorted = false)
        @test isa(N, NumericalIrreducibleDecomposition)

        N_fails = nid(F; endgame_options = EndgameOptions(; max_endgame_steps = 1))
        @test isempty(witness_sets(N_fails))

        N2 = nid(F; tracker_options = TrackerOptions(; extended_precision = false))
        @test isa(N2, NumericalIrreducibleDecomposition)

        N3 = nid(F; monodromy_options = MonodromyOptions(; trace_test_tol = 1e-12))
        @test isa(N3, NumericalIrreducibleDecomposition)
    end

    @testset "Hypersurface of degree 5" begin
        @var x[1:4]
        f = rand_poly(ComplexF64, x, 5)
        Hyp = System([f], variables = x)

        N_Hyp = numerical_irreducible_decomposition(Hyp)
        @test degrees(N_Hyp) == Dict(3 => [5])
        @test n_components(N_Hyp) == 1
    end

    @testset "Curve of degree 6" begin
        @var x[1:3]
        f = rand_poly(ComplexF64, x, 2)
        g = rand_poly(ComplexF64, x, 3)
        Curve = System([f; g], variables = x)

        N_Curve = nid(Curve)
        @test degrees(N_Curve) == Dict(1 => [6])
        @test n_components(N_Curve) == 1
    end

    @testset "Overdetermined Test" begin
        @var x y z
        TwistedCubicSphere =
            [x * z - y^2; y - z^2; x - y * z; rand_poly(ComplexF64, [x; y; z], 1)]

        N_TwistedCubicSphere = nid(TwistedCubicSphere)
        @test degrees(N_TwistedCubicSphere) == Dict(0 => [3])
        @test n_components(N_TwistedCubicSphere) == 1
    end

    @testset "Three Lines" begin
        ###  Example from 
        ### https://bertini.nd.edu/BertiniExamples/Chap8PositiveDimensional/ThreeLines.input
        @var x, y, z
        f = x * z + y
        g = y * z + x
        ThreeLines = System([f; g], variables = [x; y; z])

        N_ThreeLines = nid(ThreeLines)
        @test degrees(N_ThreeLines) == Dict(1 => [1; 1; 1])
        @test n_components(N_ThreeLines) == 3
    end

    @testset "Bricard6R" begin
        ### Example from
        ### https://bertini.nd.edu/BertiniExamples/Chap8PositiveDimensional/BricardSixR.input

        z1x = 1
        z1y = 0
        z1z = 0
        z6x = 0
        z6y = 0
        z6z = 1
        a1 = 1
        a2 = 1
        a3 = 1
        a4 = 1
        a5 = 1
        d2 = 0
        d3 = 0
        d4 = 0
        d5 = 0
        c1 = 0
        c2 = 0
        c3 = 0
        c4 = 0
        c5 = 0
        px = 0
        py = 1
        pz = 0

        @var z2x, z2y, z2z, z3x, z3y, z3z, z4x, z4y, z4z, z5x, z5y, z5z

        unit2 = z2x^2 + z2y^2 + z2z^2 - 1
        unit3 = z3x^2 + z3y^2 + z3z^2 - 1
        unit4 = z4x^2 + z4y^2 + z4z^2 - 1
        unit5 = z5x^2 + z5y^2 + z5z^2 - 1
        twist1 = z1x * z2x + z1y * z2y + z1z * z2z - c1
        twist2 = z2x * z3x + z2y * z3y + z2z * z3z - c2
        twist3 = z3x * z4x + z3y * z4y + z3z * z4z - c3
        twist4 = z4x * z5x + z4y * z5y + z4z * z5z - c4
        twist5 = z5x * z6x + z5y * z6y + z5z * z6z - c5
        # form cross products
        x1x = z1y * z2z - z1z * z2y
        x2x = z2y * z3z - z2z * z3y
        x3x = z3y * z4z - z3z * z4y
        x4x = z4y * z5z - z4z * z5y
        x5x = z5y * z6z - z5z * z6y
        x1y = z1z * z2x - z1x * z2z
        x2y = z2z * z3x - z2x * z3z
        x3y = z3z * z4x - z3x * z4z
        x4y = z4z * z5x - z4x * z5z
        x5y = z5z * z6x - z5x * z6z
        x1z = z1x * z2y - z1y * z2x
        x2z = z2x * z3y - z2y * z3x
        x3z = z3x * z4y - z3y * z4x
        x4z = z4x * z5y - z4y * z5x
        x5z = z5x * z6y - z5y * z6x
        # position
        X =
            a1 * x1x +
            d2 * z2x +
            a2 * x2x +
            d3 * z3x +
            a3 * x3x +
            d4 * z4x +
            a4 * x4x +
            d5 * z5x +
            a5 * x5x - px
        Y =
            a1 * x1y +
            d2 * z2y +
            a2 * x2y +
            d3 * z3y +
            a3 * x3y +
            d4 * z4y +
            a4 * x4y +
            d5 * z5y +
            a5 * x5y - py
        Z =
            a1 * x1z +
            d2 * z2z +
            a2 * x2z +
            d3 * z3z +
            a3 * x3z +
            d4 * z4z +
            a4 * x4z +
            d5 * z5z +
            a5 * x5z - pz

        Bricard6R = System(
            [unit2, unit3, unit4, unit5, twist1, twist2, twist3, twist4, twist5, X, Y, Z],
            variables = [z2x, z2y, z2z, z3x, z3y, z3z, z4x, z4y, z4z, z5x, z5y, z5z],
        )

        N_Bricard6R = nid(Bricard6R)
        @test degrees(N_Bricard6R) == Dict(1 => [8])
        @test n_components(N_Bricard6R) == 1
    end

    #     @testset "426" begin
    #         ### Example thanks to Julian Vill
    #         @var a, b, c, d, e


    #         f1 =
    #             -18 * a^4 * c^5 + 27 * a^3 * b * c^5 - 9 * a^2 * b^2 * c^5 - 18 * a^4 * c^4 +
    #             81 * a^3 * b * c^4 - 90 * a^2 * b^2 * c^4 + 27 * a * b^3 * c^4 - 9 * a^2 * c^6 +
    #             c^8 - 18 * a^4 * c^3 + 81 * a^3 * b * c^3 - 126 * a^2 * b^2 * c^3 +
    #             81 * a * b^3 * c^3 - 18 * b^4 * c^3 + 36 * a * b * c^5 - 9 * b^2 * c^5 +
    #             c^7 +
    #             27 * a^3 * b * c^2 - 90 * a^2 * b^2 * c^2 + 81 * a * b^3 * c^2 -
    #             18 * b^4 * c^2 + 18 * a^2 * c^4 - 36 * a * b * c^4 - 8 * c^6 -
    #             9 * a^2 * b^2 * c + 27 * a * b^3 * c - 18 * b^4 * c - 36 * a * b * c^3 +
    #             18 * b^2 * c^3 +
    #             7 * c^5 - 9 * a^2 * c^2 + 36 * a * b * c^2 - 2 * c^4 - 9 * b^2 * c + 7 * c^3 -
    #             8 * c^2 +
    #             c +
    #             1

    #         f2 =
    #             -18 * a^4 * e^5 + 27 * a^3 * d * e^5 - 9 * a^2 * d^2 * e^5 - 18 * a^4 * e^4 +
    #             81 * a^3 * d * e^4 - 90 * a^2 * d^2 * e^4 + 27 * a * d^3 * e^4 - 9 * a^2 * e^6 +
    #             e^8 - 18 * a^4 * e^3 + 81 * a^3 * d * e^3 - 126 * a^2 * d^2 * e^3 +
    #             81 * a * d^3 * e^3 - 18 * d^4 * e^3 + 36 * a * d * e^5 - 9 * d^2 * e^5 +
    #             e^7 +
    #             27 * a^3 * d * e^2 - 90 * a^2 * d^2 * e^2 + 81 * a * d^3 * e^2 -
    #             18 * d^4 * e^2 + 18 * a^2 * e^4 - 36 * a * d * e^4 - 8 * e^6 -
    #             9 * a^2 * d^2 * e + 27 * a * d^3 * e - 18 * d^4 * e - 36 * a * d * e^3 +
    #             18 * d^2 * e^3 +
    #             7 * e^5 - 9 * a^2 * e^2 + 36 * a * d * e^2 - 2 * e^4 - 9 * d^2 * e + 7 * e^3 -
    #             8 * e^2 +
    #             e +
    #             1

    #         f3 =
    #             -9 * b^2 * c^5 * d^2 * e + 27 * b * c^5 * d^3 * e - 18 * c^5 * d^4 * e +
    #             27 * b^3 * c^4 * d * e^2 - 90 * b^2 * c^4 * d^2 * e^2 +
    #             81 * b * c^4 * d^3 * e^2 - 18 * c^4 * d^4 * e^2 - 18 * b^4 * c^3 * e^3 +
    #             81 * b^3 * c^3 * d * e^3 - 126 * b^2 * c^3 * d^2 * e^3 +
    #             81 * b * c^3 * d^3 * e^3 - 18 * c^3 * d^4 * e^3 - 18 * b^4 * c^2 * e^4 +
    #             81 * b^3 * c^2 * d * e^4 - 90 * b^2 * c^2 * d^2 * e^4 +
    #             27 * b * c^2 * d^3 * e^4 - 18 * b^4 * c * e^5 + 27 * b^3 * c * d * e^5 -
    #             9 * b^2 * c * d^2 * e^5 - 9 * c^6 * d^2 * e - 9 * b^2 * c^5 * e^2 +
    #             36 * b * c^5 * d * e^2 - 36 * b * c^4 * d * e^3 +
    #             18 * c^4 * d^2 * e^3 +
    #             18 * b^2 * c^3 * e^4 - 36 * b * c^3 * d * e^4 + 36 * b * c^2 * d * e^5 -
    #             9 * c^2 * d^2 * e^5 - 9 * b^2 * c * e^6 +
    #             c^8 +
    #             c^7 * e - 8 * c^6 * e^2 + 7 * c^5 * e^3 - 2 * c^4 * e^4 + 7 * c^3 * e^5 -
    #             8 * c^2 * e^6 +
    #             c * e^7 +
    #             e^8

    #         F = System([f1; f2; f3], variables = [a, b, c, d, e])

    #         N = nid(
    #             F;
    #             seed = 0xc7cca254,
    #             endgame_options = EndgameOptions(; sing_accuracy = 1e-10),
    #         )

    #         @test n_components(N) == 1
    #         @test degrees(N) == Dict(2 => [426])
    #     end
end
