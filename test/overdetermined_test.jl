@testset "Overdetermined" begin
    @testset "square up" begin
        @var x y
        f = System([(x^2 + y^2 + x * y - 3) * (x + 3), (x^2 + y^2 + x * y  - 3) * (y - x + 2), 2x + 5y - 3])

        F = square_up(f)
        G = ModelKitSystem(System(randn(ComplexF64) .* [(x^3 - 1), y^3 - 1]))
        H = StraightLineHomotopy(G, F)
        S = collect(HC2.TotalDegreeStarts([3,3]))
        tracker = PathTracker(H)
        res = track.(tracker, S)

        @test count(is_at_infinity, res) == 4
        @test count(is_success, res) == 5

        count(is_success, map(filter(is_success, res)) do r
            newton(
                f,
                solution(r);
                extended_precision = true,
                max_rel_norm_first_update = 100r.accuracy
            )
        end)
    end
end
