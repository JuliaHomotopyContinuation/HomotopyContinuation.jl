@testset "Overdetermined" begin
    @testset "Simple system" begin
        @var x y
        f = System([
            (x^2 + y^2 + x * y - 3) * (x + 3),
            (x^2 + y^2 + x * y - 3) * (y - x + 2),
            2x + 5y - 3,
        ])
        F = square_up(f)
        H, S = total_degree_homotopy(F)
        tracker = PathTracker(H)
        res = track.(tracker, S)

        @test count(is_at_infinity, res) == 4
        @test count(is_success, res) == 5
        count(is_success, excess_solution_check(F).(res)) == 2
    end

    @testset "3 by 5 minors" begin
        F = square_up(minors())
        H, S = total_degree_homotopy(F)
        tracker = PathTracker(H, automatic_differentiation = 1)
        res = excess_solution_check(F).(track.(tracker, S))
        @test count(is_success, res) == 80
    end
end
