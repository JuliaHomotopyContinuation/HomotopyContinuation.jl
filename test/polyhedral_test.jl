@testset "Polyhedral" begin
    f = cyclic(5)
    tracker, starts = HC2.polyhedral(f)
    res = track.(tracker, starts)
    @test count(is_success, res) == 70

    f = cyclic(7)
    tracker, starts = HC2.polyhedral(f)
    @time res = track.(tracker, starts);
    @test count(is_success, res) == 924
end
