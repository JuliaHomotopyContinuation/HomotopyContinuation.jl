@testset "path_info" begin
    @polyvar x y
    f = [x^2 + y^2 + 1, x + y - 3]
    tracker, starts = coretracker_startsolutions(f)
    ctinfo = path_info(tracker, first(starts))
    @test !isempty(sprint(show, ctinfo))

    tracker = pathtracker(f)
    ptinfo = path_info(tracker, first(starts))
    @test !isempty(sprint(show, ptinfo))
end
