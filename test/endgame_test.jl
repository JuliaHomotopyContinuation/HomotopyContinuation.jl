function cyclic(n)
    @var z[1:n]
    eqs = [sum(prod(z[(k-1)%n+1] for k = j:j+m) for j = 1:n) for m = 0:(n-2)]
    push!(eqs, prod(z) - 1)
    System(eqs, z)
end

@testset "Endgame" begin

    @testset "Cyclic 7" begin
        f = cyclic(7)
        H, starts = total_degree_homotopy(f; gamma = 0.2 + 0.4im)
        S = collect(starts)
        tracker = HC2.EndgameTracker(Tracker(H))
        @time res = track.(tracker, S)

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
        tracker = HC2.EndgameTracker(Tracker(H))
        res = track.(tracker, S)
        @test count(r -> r.winding_number == 3, res) == 12

    end

    @testset "Wilkinson 19" begin

        @var x
        d = 19
        f = expand(prod(x - i for i = 1:d))
        H, starts = total_degree_homotopy([f], [x], gamma = 0.4 + 1.3im)
        S = collect(starts)
        tracker = HC2.EndgameTracker(Tracker(H))
        res = track.(tracker, S)

        @test round.(Int, real.(sort(first.(solution.(res)); by = abs))) == 1:19
        @test maximum(abs.(imag.(first.(solution.(res))))) < 1e-8
        @test count(r -> isnothing(r.winding_number), res) == 19
    end

    @testset "(x-10)^20" begin
        @var x
        f = [(x - 10)^20]
        H, starts = total_degree_homotopy(f, [x])
        S = collect(starts)
        tracker = HC2.EndgameTracker(Tracker(H))
        res = track.(tracker, S)
        @test count(r -> r.winding_number == 20, res) == 20
    end

    @testset "Beyond Polyhedral Homotopy Example" begin
        @var x y
        f = [
            2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3,
            2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5,
        ]
        H, starts = total_degree_homotopy(f, [x, y]; gamma = 1.3im + 0.4)
        S = collect(starts)
        tracker = HC2.EndgameTracker(Tracker(H))
        res = track.(tracker, S)
        @test count(is_success, res) == 2
        @test count(is_at_infinity, res) == 2
    end

    @testset "Winding Number Family" begin
        for d = 2:6
            d = 6
            @var x y
            a = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32]
            f1 = (a[1] * x^d + a[2] * y) * (a[3] * x + a[4] * y) + 1
            f2 = (a[1] * x^d + a[2] * y) * (a[5] * x + a[6] * y) + 1
            H, starts =
                total_degree_homotopy([f1, f2], [x, y]; gamma = 1.3im + 0.4)
            S = collect(starts)
            tracker = HC2.EndgameTracker(Tracker(H))
            res = track.(tracker, S)
            @test count(is_success, res) == d + 1
            @test count(is_at_infinity, res) == (d + 1)^2 - d - 1
        end
    end
end
