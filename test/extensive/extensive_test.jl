@testset "Lines on a Quintic surface in 3-space" begin
    n = 4
    @var x[1:n]

    ν = monomials(x, 5)
    N = length(ν)

    q₀ = randn(ComplexF64, N-1)
    @var q[1:N-1]
    F = sum(q[i] * ν[i] for i in 1:N-1) + 1

    @var a[1:n-1] b[1:n-1] t
    L = [a;1] .* t + [b;0]
    FcapL = last(ModelKit.exponents_coefficients(subs(F, x=>L), [t]))
    sys = System(ModelKit.horner.(FcapL), [a;b], q)
    H, starts = total_degree_homotopy(sys; gamma = 0.42 - 0.52im, target_parameters = q₀)
    tracker = HC2.EndgameTracker(Tracker(H))
    @time res = track.(tracker, starts)
    @test count(is_success, res) == 2875
end
