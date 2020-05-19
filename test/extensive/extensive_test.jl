@testset "Lines on a Quintic surface in 3-space" begin
    @testset "Lines on a quintic surface in 3-space" begin
        sys, q₀ = fano_quintic()
        @time res =
            solve(sys; gamma = gamma, target_parameters = q₀, start_system = :total_degree)
        @test nsolutions(res) == 2875

        @time poly_res = solve(sys; target_parameters = q₀, start_system = :polyhedral)
        @test nsolutions(poly_res) = 2875
    end

end


@time poly_res = solve(sys; target_parameters = q₀, start_system = :polyhedral,
    # tracker_options = (parameters = :conservative,),
    # path_tracker_options = PathTrackerOptions(val_at_infinity_tol = 1e-8),
        seed = 0x56db29a3)


# succ = path_number.(results(poly_res))

succ2 = path_number.(results(poly_res))


setdiff(succ, succ2)

solver, starts =solver_startsolutions(sys; target_parameters = q₀, start_system = :polyhedral,
   # tracker_options = (parameters = :conservative,),
   # path_tracker_options = PathTrackerOptions(val_at_infinity_tol = 1e-8),
       seed = 0x56db29a3)
S = collect(starts)


solver.tracker.generic_tracker.options.min_coord_growth = 1000
track(solver.tracker, S[1785]; debug = true)



F = six_revolute()

res = solve(F; start_system = :total_degree)
a x
