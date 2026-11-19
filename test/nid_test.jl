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

        dec = decompose(W)
        @test all(is_irreducible, dec)

        # sorting
        W = regeneration(F; sorted = false, show_progress = false)
        @test degree.(W) == [2, 8, 8]

        # limited codimension
        W = regeneration(F; max_codim = 2, show_progress = false)
        @test degree.(W) == [2, 8]

        # no threading 
        N = nid(F; threading = false, show_progress = false)
        @test isa(N, NumericalIrreducibleDecomposition)

        # bad seed
        s = UInt32(62)
        N = nid(F; seed = s, show_progress = false)
        degs = degrees(N)
        @test seed(N) == s
        @test degs[2] == [2]
        @test degs[1] == [4, 4]

        N = nid(F; seed = nothing, show_progress = false)
        @test isnothing(seed(N))

        # options
        N_fails = nid(
            F;
            endgame_options = EndgameOptions(; max_endgame_steps = 0),
            show_progress = false,
        )
        @test isempty(witness_sets(N_fails))

        N2 = nid(
            F;
            monodromy_options = MonodromyOptions(; trace_test_tol = 1e-5),
            show_progress = false,
        )
        @test isa(N2, NumericalIrreducibleDecomposition)

        N = nid(F; show_monodromy_progress = true, show_progress = false)
        @test isa(N, NumericalIrreducibleDecomposition)

        N = nid(F; warning = false, show_progress = false)
        @test isa(N, NumericalIrreducibleDecomposition)

        # number of components
        @test ncomponents(N2) == 11
        @test ncomponents(N2, dims = [1, 2]) == 3
        @test ncomponents(N2, 1) == 2
        @test n_components(N2) == 11
        @test n_components(N2, dims = [1, 2]) == 3
        @test n_components(N2, 1) == 2
    end

    @testset "randomization" begin
        @var x y
        f = System([(x - 1)^2; x - 1], variables = [x, y])
        N = nid(f; sorted = :randomized, show_progress = false)
        degs = degrees(N)
        @test degs[1] == [1]
    end

    @testset "rational systems" begin
        @var x y z
        g = System([x^2 + y^2 - z; x / (y - 1) + y + z - 1], variables = [x, y, z])
        N = nid(g)
        W = first(witness_sets(N)[1])
        @test ncomponents(N) == 1
        @test degree(W) == 4
        p = randn(3)
        q = solutions(W)[1]
        @test !membership(p, W)
        @test membership(q, W)
        U = intersect(W, x / (y - 1) + y + z - 1)
        @test degree(U) == 4
        V = intersect(W, x * y^3 - z^4 + x^3 - 8)
        @test degree(V) == 16
    end

    @testset "Hypersurface of degree 5" begin
        @var x[1:4]
        f = rand_poly(ComplexF64, x, 5)
        Hyp = System([f], variables = x)

        N_Hyp = numerical_irreducible_decomposition(Hyp, show_progress = false)
        @test degrees(N_Hyp) == Dict(3 => [5])
        @test ncomponents(N_Hyp) == 1
    end

    @testset "Curve of degree 6" begin
        @var x[1:3]
        f = rand_poly(ComplexF64, x, 2)
        g = rand_poly(ComplexF64, x, 3)
        Curve = System([f; g], variables = x)

        N_Curve = nid(Curve, show_progress = false)
        @test degrees(N_Curve) == Dict(1 => [6])
        @test ncomponents(N_Curve) == 1
    end

    @testset "Homogeneous systems" begin
        @var x[1:4]

        a = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
        b = x[1]^3 + x[2]^3 + 2x[3]^3 + 3x[4]^3
        c = x[1]^4 + 2x[2]^4 + 4x[3]^4 - x[4]^4
        G = System([a * c; b * c]; variables = x)
        N = nid(G; show_progress = false)
        @test degrees(N) == Dict(2 => [4], 1 => [6])
        @test ncomponents(N) == 2

        Ws = witness_sets(N)
        W1, W2 = Ws[1][1], Ws[2][1]
        @test all(W -> W.projective, [W1; W2])
        @test all(W -> is_linear(linear_subspace(W)), [W1; W2])
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

        N_ThreeLines = nid(ThreeLines, show_progress = false)
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

        N_Bricard6R = nid(Bricard6R, show_progress = false)
        @test degrees(N_Bricard6R) == Dict(1 => [8])
        @test ncomponents(N_Bricard6R) == 1

        N_Bricard6R_c4 = nid(Bricard6R, max_codim = 4, show_progress = false)
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

        N_ACR = nid(F_ACR, show_progress = false)
        degrees(N_ACR) == Dict(4 => [7])
    end


    @testset "Union of a sphere, a line, and a point" begin
        @var x, y, z

        S = [x^2 + y^2 + z^2 - 1] #The sphere
        L = [2 * x - z, 2 * y - z] #The line (t,t,2t)
        P = [x + y + 2 * z - 4, y - z, x - z] #The point (1,1,1)

        F = System([s * l * p for s in S for l in L for p in P])

        NID =
            numerical_irreducible_decomposition(F; seed = 0x7a4845b9, show_progress = false)

        @test ncomponents(NID, 0) == 1
        @test degrees(NID)[0] == [1]
    end

    @testset "trace test with cached point traces" begin
        @var x y
        circle, cubic = x^2 + y^2 - 1, y - x^3 + 2x - 1
        line1, line2 = x - 2y, x + y - 3
        function setup(f, seed)
            W = only(regeneration(System([f]); show_progress = false, seed = seed))
            opts = HC.decompose_with_monodromy_options(MonodromyOptions())
            L = linear_subspace(W)
            MS = HC.MonodromySolver(system(W), L; compile = false, options = opts)
            PT = HC.PointTraces(L, opts) # one fixed displacement for all points
            HC.update_traces!(PT, solutions(W), MS, false)
            W, solutions(W), PT, opts.trace_test_tol
        end
        on(g, P) = findall(p -> abs(System([g], variables = [x, y])(p)[1]) < 1e-8, P)
        function certify(W, P, PT, tol, orbits)
            certified = Vector{Int}[]
            up = HC.UniquePoints(first(P), 1; distance = HC.InfNorm())
            dec = WitnessSet[]
            HC.certify_orbit_unions!(
                Set{Int}[],
                dec,
                up,
                orbits,
                Set{Int}(),
                P,
                PT,
                system(W),
                linear_subspace(W),
                tol;
                atol = 1e-14,
                rtol = 1e-8,
            ) do idx, _
                push!(certified, idx)
            end
            sort(certified), dec
        end

        # two components, each split into several orbits: exactly the components are certified
        W, P, PT, tol = setup(circle * cubic, 0x2d4c1a73)
        A, B = on(circle, P), on(cubic, P)
        trace(S) = HC.cached_trace(PT, P[S])
        @test trace(A) < tol && trace(B) < tol && trace(B[1:2]) > tol
        orbits = [Set([A[1]]), Set([A[2]]), Set(B[1:2]), Set([B[3]])]
        certified, dec = certify(W, P, PT, tol, orbits)
        @test certified == sort([A, B]) && sort(degree.(dec)) == [2, 3]

        # the union of the two components has a vanishing trace but is not minimal
        AB = sort([A; B])
        @test HC.is_certifiable(PT, P[A], trace(A), tol)
        @test !HC.is_certifiable(PT, P[AB], trace(AB), tol)
        @test Set(HC.split_nonminimal_orbits([Set(AB)], P, PT, Set{Int}(), tol)) ==
              Set([Set(A), Set(B)])
        @test isempty(
            first(certify(W, P, PT, tol, [Set([A[1], B[1]]), Set([A[2]; B[2:3]])])),
        )

        # points on lines are components on their own: they are not combined, and a single point needs the degree 1 rule; links of such points need a second loop
        W, P, PT, tol = setup(line1 * line2 * circle, 0x5be0f3a1)
        l1, l2, C = only(on(line1, P)), only(on(line2, P)), on(circle, P)
        @test HC.small_trace_points(PT, P, tol) == Set([l1, l2])
        certified, _ =
            certify(W, P, PT, tol, [Set([l1]), Set([l2]), Set([C[1]]), Set([C[2]])])
        @test certified == [sort(C)]
        @test !HC.is_certifiable(PT, P[[l1]], trace([l1]), tol)
        @test HC.is_certifiable(PT, P[[l1]], trace([l1]), tol; allow_single_point = true)
        jump = reshape([i == l1 ? l2 : i == l2 ? l1 : i for i = 1:length(P)], :, 1)
        links = HC.SmallTraceLinks()
        orbits, _ = HC.get_orbits_from_monodromy_permutations(
            jump,
            Set{Int}[],
            P,
            Set([l1, l2]),
            links,
        )
        @test only(filter(o -> l1 in o, orbits)) == Set([l1])
        orbits, _ =
            HC.get_orbits_from_monodromy_permutations(jump, orbits, P, Set([l1, l2]), links)
        @test only(filter(o -> l1 in o, orbits)) == Set([l1, l2])
    end

    @testset "trace_precheck!" begin
        @var x y
        W = only(
            regeneration(
                System([(x^2 + y^2 - 1) * (y - x^3 + 2x - 1)]);
                show_progress = false,
                seed = 0x2d4c1a73,
            ),
        )
        L = linear_subspace(W)
        opts = HC.decompose_with_monodromy_options(MonodromyOptions())
        MS = HC.MonodromySolver(system(W), L; compile = false, options = opts)
        # complete witness points are certified, incomplete ones are not
        @test HC.trace_precheck!(MS, solutions(W), L, false)
        @test !HC.trace_precheck!(MS, solutions(W)[1:3], L, false)
        @test is_success(
            monodromy_solve(MS, solutions(W)[1:3], L, rand(UInt32); show_progress = false),
        )
    end
end
