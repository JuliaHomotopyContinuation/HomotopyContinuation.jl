@testset "Overdetermined systems" begin
        @polyvar x y z
        f = [x * z - y^2, y - z^2, x - y * z, x + y + z + 1]

        tracker, starts = pathtracker_startsolutions(f; seed = 213412)
        @test count(is_success, track.(tracker, starts)) == 3
        @test count(r -> r.return_code == :excess_solution, track.(tracker, starts)) == 2

        # check that it also works with polyhedral homotopy
        tracker, starts = pathtracker_startsolutions(
                f;
                seed = 213412,
                start_system = :polyhedral,
        )
        HC.prepare!(tracker, starts)
        @test count(is_success, track.(tracker, starts)) == 3

        minors = include(joinpath(@__DIR__, "examples", "3_by_5_minors.jl"))
        tracker, starts = pathtracker_startsolutions(minors)
        @test count(is_success, track.(tracker, starts)) == 80

        # chemical reaction network examples
        g = let
                @polyvar x[1:3]
                p = [0.04, 0.04, 1.0, 1.0, 10.0, 0.0, 0.04, 35.0, .1, .04]
                [
                 -x[1] * x[3] * p[3] - x[1] * p[2] + p[1],
                 x[1] * x[2] * x[3] * p[8] * p[9] + x[1] * x[3] * p[7] * p[8] * p[9] -
                 x[2] * x[3] * p[5] * p[6] - x[2] * p[5] * p[6] * p[10] -
                 x[2] * x[3] * p[4] - x[2] * p[4] * p[10],
                 -x[1] * x[2] * x[3] * p[8] * p[9] - x[1] * x[3] * p[7] * p[8] * p[9] +
                 x[2] * x[3] * p[5] * p[6] + x[2] * p[5] * p[6] * p[10] +
                 x[2] * x[3] * p[4] + x[2] * p[4] * p[10],
                 x[2] + x[3] - 1.0,
                ]

        end

        tracker, starts = pathtracker_startsolutions(g; seed = 213412)
        @test count(is_success, track.(tracker, starts)) == 3
end
