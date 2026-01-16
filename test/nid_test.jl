@testset "Numerical Irreducible Decomposition" begin
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

        W = regeneration(F; max_codim = 2)
        @test degree.(W) == [2, 8]

        dec = decompose(W)
        N = NumericalIrreducibleDecomposition(dec)
        @test isa(N, NumericalIrreducibleDecomposition)

        # no threading 
        N = nid(F; threading = false)
        @test isa(N, NumericalIrreducibleDecomposition)

        # seed
        s = 0x42c9d504
        N = nid(F; seed = s)
        @test seed(N) == s

        N = nid(F; seed = nothing)
        @test isnothing(seed(N))

        # bad seed
        N = nid(F; seed = 0xc770fa47)
        degs = degrees(N)
        @test degs[2] == [2]
        @test degs[1] == [4, 4]

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

        N_fails = nid(F; endgame_options = EndgameOptions(; max_endgame_steps = 0))
        @test isempty(witness_sets(N_fails))

        N2 = nid(F; tracker_options = TrackerOptions(; extended_precision = false))
        @test isa(N2, NumericalIrreducibleDecomposition)

        N3 = nid(F; monodromy_options = MonodromyOptions(; trace_test_tol = 1e-12))
        @test isa(N3, NumericalIrreducibleDecomposition)

        # number of components
        @test ncomponents(N3) == 11
        @test ncomponents(N3, dims = [1, 2]) == 3
        @test ncomponents(N3, 1) == 2
        @test n_components(N3) == 11
        @test n_components(N3, dims = [1, 2]) == 3
        @test n_components(N3, 1) == 2

        # max_codim = 1
        N4 = nid(F; max_codim = 1)
        @test isa(N4, NumericalIrreducibleDecomposition)
    end

    @testset "Non complete superwitness set after regeneration" begin
        @var p1, p2, p3, p4, p5, p6, p7, p8, u1, u2, u3, u4, u5, u6, u7, u8, λ1, λ2, λ3
        p = [p1, p2, p3, p4, p5, p6, p7, p8]
        u = [u1, u2, u3, u4, u5, u6, u7, u8]
        λ = [λ1, λ2, λ3]
        E = System(
            [
                u1 - p1 * (λ3 - p6 * λ1 - p7 * λ2),
                u2 - p2 * (λ3 + p5 * λ1 + p8 * λ2),
                u3 - p3 * (λ3 + p5 * λ2 + p8 * λ1),
                u4 - p4 * (λ3 - p6 * λ2 - p7 * λ1),
                u5 - p5 * (λ3 + p2 * λ1 + p3 * λ2),
                u6 - p6 * (λ3 - p1 * λ1 - p4 * λ2),
                u7 - p7 * (λ3 - p1 * λ2 - p4 * λ1),
                u8 - p8 * (λ3 + p2 * λ2 + p3 * λ1),
                p2 * p5 + p3 * p8 - (p1 * p6 + p4 * p7),
                p2 * p8 + p5 * p3 - (p1 * p7 + p6 * p4),
                -1 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8,
                u2 * u5 + u8 * u3 - (u1 * u6 + u7 * u4),
                u3 * u5 + u8 * u2 - (u1 * u7 + u6 * u4),
                -1 + u1 + u2 + u3 + u4 + u5 + u6 + u7 + u8,
            ],
            variables = [p; u; λ],
        )
        s = 0x7eec3900
        R = regeneration(E; seed = s)
        D = decompose(R; seed = s)
        N = NumericalIrreducibleDecomposition(D, s)
        @test ncomponents(N) == 6
    end

    @testset "Hypersurface of degree 5" begin
        @var x[1:4]
        f = rand_poly(ComplexF64, x, 5)
        Hyp = System([f], variables = x)

        N_Hyp = numerical_irreducible_decomposition(Hyp)
        @test degrees(N_Hyp) == Dict(3 => [5])
        @test ncomponents(N_Hyp) == 1
    end

    @testset "Curve of degree 6" begin
        @var x[1:3]
        f = rand_poly(ComplexF64, x, 2)
        g = rand_poly(ComplexF64, x, 3)
        Curve = System([f; g], variables = x)

        N_Curve = nid(Curve)
        @test degrees(N_Curve) == Dict(1 => [6])
        @test ncomponents(N_Curve) == 1
    end

    @testset "Overdetermined Test" begin
        @var x y z
        TwistedCubicSphere =
            [x * z - y^2; y - z^2; x - y * z; rand_poly(ComplexF64, [x; y; z], 1)]

        N_TwistedCubicSphere = nid(TwistedCubicSphere)
        @test degrees(N_TwistedCubicSphere) == Dict(0 => [1, 1, 1])
        @test ncomponents(N_TwistedCubicSphere) == 3
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
        @test ncomponents(N_ThreeLines) == 3
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
        @test ncomponents(N_Bricard6R) == 1

        N_Bricard6R_c4 = nid(Bricard6R, max_codim = 4)
        @test ncomponents(N_Bricard6R_c4) == 0
    end

    @testset "ACR" begin
        @var xx_Di, xx_Da, xx_Ya, xx_Yi
        @var xx_CXY, xx_CXYp, xx_G, xx_CNA, xx_X
        @var xx_A, xx_N, xx_CXA, xx_T, xx_CXT, xx_Xp


        F_ACR = [
            -(5 / 3) * xx_Di + (2 / 3) * xx_Da,
            (5 / 3) * xx_Di - (2 / 3) * xx_Da,
            -(1 / 3) * xx_Ya * xx_X - (2 / 3) * xx_Ya +
            (5 / 6) * xx_Yi +
            (6 / 5) * xx_CXY +
            4 * xx_CXYp,
            -8 * xx_Da * xx_Yi + (5 / 7) * xx_G * xx_CNA + (2 / 3) * xx_Ya -
            (49 / 30) * xx_Yi,
            8 * xx_Da * xx_Yi - (5 / 7) * xx_G * xx_CNA + (4 / 5) * xx_Yi,
            8 * xx_Da * xx_Yi - (5 / 7) * xx_G * xx_CNA +
            (2 / 3) * xx_A * xx_N +
            (4 / 5) * xx_Yi - (8 / 3) * xx_CNA,
            -(2 / 3) * xx_A * xx_N - (7 / 8) * xx_A * xx_X +
            (8 / 3) * xx_CNA +
            (5 / 7) * xx_CXA,
            -(2 / 3) * xx_A * xx_N + (8 / 3) * xx_CNA - (5 / 2) * xx_N + 1,
            -(1 / 3) * xx_Ya * xx_X - (7 / 8) * xx_A * xx_X - xx_X * xx_T - (4 / 3) * xx_X +
            (6 / 5) * xx_CXY +
            3 * xx_CXT +
            (5 / 7) * xx_CXA +
            1,
            (1 / 3) * xx_Ya * xx_X - (26 / 5) * xx_CXY,
            4 * xx_CXY - 4 * xx_CXYp,
            4 * xx_CXYp - 5 * xx_Xp,
            -xx_X * xx_T + 3 * xx_CXT,
            xx_X * xx_T - 3 * xx_CXT,
            (7 / 8) * xx_A * xx_X - (5 / 7) * xx_CXA,
        ]

        N_ACR = nid(F_ACR)
        degrees(N_ACR) == Dict(4 => [7])
    end


    @testset "Union of a sphere, a line, and a point" begin
        @var x, y, z

        S = [x^2 + y^2 + z^2 - 1] #The sphere
        L = [2 * x - z, 2 * y - z] #The line (t,t,2t)
        P = [x + y + 2 * z - 4, y - z, x - z] #The point (1,1,1)

        F = System([s * l * p for s in S for l in L for p in P])

        NID = numerical_irreducible_decomposition(F; seed = 0x7a4845b9)

        @test ncomponents(NID, 0) == 1
        @test degrees(NID)[0] == [1]
    end
end
