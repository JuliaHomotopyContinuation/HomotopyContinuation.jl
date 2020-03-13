function cyclic(n)
    @var z[1:n]
    eqs = [sum(prod(z[(k-1)%n+1] for k = j:j+m) for j = 1:n) for m = 0:(n-2)]
    push!(eqs, prod(z) - 1)
    System(eqs, z)
end

function bacillus()
    @var w w2 w2v v w2v2 vP sigmaB w2sigmaB vPp phos
    poly = [
        (-1 * 0.7 * w + -2 * 3600.0 * (w^2 / 2) + 2 * 18.0 * w2) *
        (0.2 + sigmaB) + 4.0 * 0.4 * (1 + 30.0sigmaB),
        -1 * 0.7 * w2 +
        3600.0 * (w^2 / 2) +
        -1 * 18.0 * w2 +
        -1 * 3600.0 * w2 * v +
        18.0w2v +
        36.0w2v +
        -1 * 3600.0 * w2 * sigmaB +
        18.0w2sigmaB,
        -1 * 0.7 * w2v +
        3600.0 * w2 * v +
        -1 * 18.0 * w2v +
        -1 * 3600.0 * w2v * v +
        18.0w2v2 +
        -1 * 36.0 * w2v +
        36.0w2v2 +
        1800.0 * w2sigmaB * v +
        -1 * 1800.0 * w2v * sigmaB,
        (
            -1 * 0.7 * v +
            -1 * 3600.0 * w2 * v +
            18.0w2v +
            -1 * 3600.0 * w2v * v +
            18.0w2v2 +
            -1 * 1800.0 * w2sigmaB * v +
            1800.0 * w2v * sigmaB +
            180.0vPp
        ) * (0.2 + sigmaB) + 4.5 * 0.4 * (1 + 30.0sigmaB),
        -1 * 0.7 * w2v2 +
        3600.0 * w2v * v +
        -1 * 18.0 * w2v2 +
        -1 * 36.0 * w2v2,
        -1 * 0.7 * vP + 36.0w2v + 36.0w2v2 + -1 * 3600.0 * vP * phos + 18.0vPp,
        (
            -1 * 0.7 * sigmaB +
            -1 * 3600.0 * w2 * sigmaB +
            18.0w2sigmaB +
            1800.0 * w2sigmaB * v +
            -1 * 1800.0 * w2v * sigmaB
        ) * (0.2 + sigmaB) + 0.4 * (1 + 30.0sigmaB),
        -1 * 0.7 * w2sigmaB +
        3600.0 * w2 * sigmaB +
        -1 * 18.0 * w2sigmaB +
        -1 * 1800.0 * w2sigmaB * v +
        1800.0 * w2v * sigmaB,
        -1 * 0.7 * vPp +
        3600.0 * vP * phos +
        -1 * 18.0 * vPp +
        -1 * 180.0 * vPp,
        (phos + vPp) - 2.0,
    ]
    System(poly, [w, w2, w2v, v, w2v2, vP, sigmaB, w2sigmaB, vPp, phos])
end

@testset "Endgame AD4: $AD4" for AD4 in [true, false]
    @testset "Cyclic 7" begin
        f = cyclic(7)
        H, starts = total_degree_homotopy(f; gamma = 0.2 + 0.4im)
        S = collect(starts)
        tracker = HC2.EndgameTracker(Tracker(
            H,
            automatic_differentiation = (true, true, true, AD4),
        ))
        res = track.(tracker, S)

        @test count(is_success, res) == 924
        @test count(is_at_infinity, res) == 4116
    end

    @testset "Hyperbolic - 6,6" begin
        # 2 solutions with multiplicity 6, projective
        @var x z
        y = 1
        # This has two roots of multiplicity 6 at the hyperplane z=0
        # each root has winding number 3
        F = [
            0.75 * x^4 + 1.5 * x^2 * y^2 - 2.5 * x^2 * z^2 + 0.75 * y^4 -
            2.5 * y^2 * z^2 + 0.75 * z^4
            10 * x^2 * z + 10 * y^2 * z - 6 * z^3
        ]
        H, starts = total_degree_homotopy(F, [x, z])
        S = collect(starts)
        tracker = HC2.EndgameTracker(Tracker(
            H,
            automatic_differentiation = (true, true, true, AD4),
        ))
        res = track.(tracker, S)
        @test count(r -> r.winding_number == 3, res) == 12
    end

    @testset "Wilkinson 19" begin
        @var x
        d = 19
        f = expand(prod(x - i for i = 1:d))
        H, starts = total_degree_homotopy([f], [x], gamma = 0.4 + 1.3im)
        S = collect(starts)
        tracker = HC2.EndgameTracker(Tracker(
            H,
            automatic_differentiation = (true, true, true, AD4),
        ))
        res = track.(tracker, S)
        @test all(is_success, res)
        @test round.(Int, real.(sort(first.(solution.(res)); by = abs))) == 1:19
        @test maximum(abs.(imag.(first.(solution.(res))))) < 1e-8
        @test count(r -> isnothing(r.winding_number), res) == 19
    end

    @testset "(x-10)^$d" for d in [2, 8, 12, 16]
        @var x
        f = [(x - 10)^d]
        H, starts = total_degree_homotopy(f, [x])
        S = collect(starts)
        tracker = HC2.EndgameTracker(Tracker(
            H,
            automatic_differentiation = (true, true, true, AD4),
        ))
        res = track.(tracker, S)
        @test count(r -> r.winding_number == d, res) == d
    end

    @testset "Beyond Polyhedral Homotopy Example" begin
        @var x y
        f = [
            2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3,
            2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5,
        ]
        H, starts = total_degree_homotopy(f, [x, y]; gamma = 1.3im + 0.4)
        S = collect(starts)
        tracker = HC2.EndgameTracker(Tracker(
            H,
            automatic_differentiation = (true, true, true, AD4),
        ))
        res = track.(tracker, S)
        @test count(is_success, res) == 2
        @test count(is_at_infinity, res) == 2
    end

    @testset "Winding Number Family d=$d" for d = 2:2:6
        @var x y
        a = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32]
        f1 = (a[1] * x^d + a[2] * y) * (a[3] * x + a[4] * y) + 1
        f2 = (a[1] * x^d + a[2] * y) * (a[5] * x + a[6] * y) + 1
        H, starts =
            total_degree_homotopy([f1, f2], [x, y]; gamma = 1.3im + 0.4)
        S = collect(starts)
        tracker = HC2.EndgameTracker(Tracker(
            H,
            automatic_differentiation = (true, true, true, AD4),
        ))
        res = track.(tracker, S)
        @test count(is_success, res) == d + 1
        @test count(is_at_infinity, res) == (d + 1)^2 - d - 1
    end

    @testset "Bacillus Subtilis" begin
        H, starts = total_degree_homotopy(bacillus())
        tracker = HC2.EndgameTracker(Tracker(H))
        tracker = HC2.EndgameTracker(Tracker(
            H,
            automatic_differentiation = (true, true, true, AD4),
        ))
        res = track.(tracker, starts)
        @test count(is_success, res) == 44
    end

    @testset "Mohab" begin
        # Communicated by Mohab Safey El Din
        @var x y z

        F =
            ModelKit.horner.([
                -9091098778555951517 * x^3 * y^4 * z^2 +
                5958442613080401626 * y^2 * z^7 +
                17596733865548170996 * x^2 * z^6 -
                17979170986378486474 * x * y * z^6 -
                2382961149475678300 * x^4 * y^3 -
                15412758154771986214 * x * y^3 * z^3 + 133,
                -10798198881812549632 * x^6 * y^3 * z -
                11318272225454111450 * x * y^9 -
                14291416869306766841 * y^9 * z -
                5851790090514210599 * y^2 * z^8 +
                15067068695242799727 * x^2 * y^3 * z^4 +
                7716112995720175148 * x^3 * y * z^3 +
                171,
                13005416239846485183 * x^7 * y^3 +
                4144861898662531651 * x^5 * z^4 - 8026818640767362673 * x^6 -
                6882178109031199747 * x^2 * y^4 +
                7240929562177127812 * x^2 * y^3 * z +
                5384944853425480296 * x * y * z^4 +
                88,
            ])

        H, starts =
            total_degree_homotopy(System(F, [x, z, y]); gamma = 0.2 + 0.4im)
        egtracker = HC2.EndgameTracker(Tracker(
            H,
            automatic_differentiation = (true, true, true, AD4),
        ))
        egres = HC2.track.(egtracker, starts)
        @test count(is_success, egres) == 693
    end
end
