module TestSystems

    import TypedPolynomials
    using ..HomConBase

    export cyclic5, cyclic5Solutions, cyclic7, cyclic7Solutions

    """
        cyclical(iter, n)

    Cycles through an `iter` with n consecutive elements.

    ###Example:
        cyclical([1, 2, 3, 4], 2) == [[1, 2], [2, 3], [3, 4], [4, 1]]
    """
    function cyclical(iter, n)
        m = length(iter)
        values = collect(iter);
        res = Vector{Vector{eltype(values)}}()
        for k=0:m-1
            push!(res, [values[((k + i) % m) + 1] for i=0:n-1])
        end
        res
    end

    function cyclical_polys(vars)
        n = length(vars)
        F = map(k -> sum(map(prod, cyclical(vars, k))), 1:n-1)
        push!(F, prod(vars) - 1)
        PolySystem(F)
    end

    """
        cyclic5()

    Cyclic5 example problem folowing [^1]

    [^1]: A faster way to count the solutions of inhomogeneous systems of algebraic equations,
    with applications to cyclic n-roots
    """
    function cyclic5()
        TypedPolynomials.@polyvar z[1:5]
        cyclical_polys(z)
    end

    """
        cyclic5Solutions()

    Solutions to the Cyclic5 example.
    """
    function cyclic5Solutions()
        # (1,w,w^{2k},w^{3k},w^{4k}) with w = exp(2πi/5) permuted cyclical, 1 ≦ k ≦ 4
        w = exp(2.0*π*im/5.0)
        classic_sols = vec([ map(i -> w^(i*k), perm) for k in 1:4, perm in cyclical(0:4, 5)])

        #solution of ɛ^2+3ɛ+1=0
        ɛ = -1.5 - 0.5 * √5.0 + 0*im
        extended_sols = vec([ rev ? reverse!(w^k * perm) : w^k * perm for perm in cyclical([ɛ, 1/ɛ, 1, 1, 1], 5),
                                                                            k in 0:4,
                                                                            rev in [false, true]])
        vcat(classic_sols, extended_sols)
    end
    """
        cyclic7()

    Cyclic7 problem folowing [^1]

    [^1]: A faster way to count the solutions of inhomogeneous systems of algebraic equations,
    with applications to cyclic n-roots
    """
    function cyclic7()
        TypedPolynomials.@polyvar z[1:7]
        cyclical_polys(z)
    end

    """
        cyclic7Solutions()

    Solutions to the Cyclic7 example.
    """
    function cyclic7Solutions()
        # (1,w,w^{2k},w^{3k},w^{4k},w^{5k},w^{6k}) with w = exp(2πi/7) permuted cyclical, 1 ≦ k ≦ 6
        w = exp(2.0*π*im/7.0)
        classic_sols = vec([ map(i -> w^(i*k), perm) for k in 1:6, perm in cyclical(0:6, 7)])

        #solution of ɛ^2+5ɛ+1=0
        ɛ = -2.5 - 0.5 * √21.0 + 0*im
        extended_sols = vec([ rev ? reverse!(w^k * perm) : w^k * perm for perm in cyclical([ɛ, 1/ɛ, 1, 1, 1, 1, 1], 7),
                                                                            k in 0:6,
                                                                            rev in [false, true]])

        # index two solutions
        d = 0.25*im*(√7 + 3*im)
        index_two = vec([ rev ? reverse!(w^k * perm) : w^k * perm for perm in cyclical([d, conj(d), d, 1, conj(d), 1, 1], 7),
                                                                            k in 0:6,
                                                                            rev in [false, true]])
        index_two_conj = vec([ rev ? reverse!(w^k * perm) : w^k * perm for perm in cyclical([conj(d), d, conj(d), 1, d, 1, 1], 7),
                                                                            k in 0:6,
                                                                            rev in [false, true]])

        # index three solutions
        a = 2.738895317095 + 0im
        b = 3.436680125767 + 0im
        c = -0.1298393513967 + 0im
        sol_1 = [a, b, c, 1, 1/c, 1/b, 1/a]
        sol_2 = [a*b,c,1/(b*c),1,b*c,1/c, 1/(a*b)]
        sol_3 = [a*b*c,1/(b*c),b,1,1/b,b*c,1/(a*b*c)]

        perms = [cyclical(sol_1, 7); cyclical(sol_2, 7); cyclical(sol_3, 7)]
        index_three_first = vec([ rev ? reverse!(w^k * perm) : w^k * perm for perm in perms, k in 0:6, rev in [false, true]])

        a = exp(4.3128389787245im)
        b = exp(1.356227956787im)
        c = exp(1.900668281165im)
        sol_1 = [a, b, c, 1, 1/c, 1/b, 1/a]
        sol_2 = [a*b,c,1/(b*c),1,b*c,1/c, 1/(a*b)]
        sol_3 = [a*b*c,1/(b*c),b,1,1/b,b*c,1/(a*b*c)]
        perms = [cyclical(sol_1, 7); cyclical(sol_2, 7); cyclical(sol_3, 7)]
        index_three_second = vec([ rev ? reverse!(w^k * perm) : w^k * perm for perm in perms, k in 0:6, rev in [false, true]])

        vcat(classic_sols, extended_sols, index_two, index_two_conj, index_three_first, index_three_second)
    end
end
