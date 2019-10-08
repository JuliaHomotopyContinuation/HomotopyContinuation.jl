@testset "PathTracker" begin
    @testset "Constructor + Options" begin
        @polyvar x y z
        f = [x^2 - 2, x + y - 1]
        tracker, starts = pathtracker_startsolutions(f; system = FPSystem)
        @test tracker isa PathTracker{Vector{ComplexF64}}
        test_show_juno(tracker)
        @test !isempty(sprint(show, tracker))
        test_show_juno(tracker.options)
        @test !isempty(sprint(show, tracker.options))
        test_show_juno(tracker.state)
        @test !isempty(sprint(show, tracker.state))
        @test Base.broadcastable(tracker) isa Ref{typeof(tracker)}
        tracker = pathtracker(f; system = FPSystem, precision_strategy = :adaptive_never)
        @test tracker.options.precision_strategy == HC.PREC_STRATEGY_NEVER
        tracker = pathtracker(f; system = FPSystem, precision_strategy = :adaptive_always)
        @test tracker.options.precision_strategy == HC.PREC_STRATEGY_ALWAYS
        @test_throws ArgumentError pathtracker(
            f;
            system = FPSystem,
            precision_strategy = :NOT_AN_OPTION,
        )

        @test_throws ArgumentError pathtracker(
            f;
            not_an_option = true,
            this_option = "nope",
        )

        @test_throws ArgumentError pathtracker(f, patch = OrthogonalPatch())
        @test_throws ArgumentError pathtracker(f, system_scaling = :lalala)
    end

    # Simple
    @testset "Simple" begin
        @polyvar x y z
        f = [x^2 - 2, x + y - 1]
        tracker, starts = pathtracker_startsolutions(f; system = FPSystem)

        for xᵢ in starts
            @test is_success(track!(tracker, xᵢ))
            s = solution(tracker)
            @test abs(s[1]) ≈ sqrt(2) atol = 1e-8
            @test sum(s) ≈ 1 atol = 1e-8
            @test isnothing(winding_number(tracker))
            @test !isnan(tracker.state.solution_cond)
            @test !isnan(tracker.state.solution_accuracy)
            @test tracker.state.solution_accuracy < 1e-12
            @test tracker.state.solution_residual < 1e-12
        end

        # test iterator
        init!(tracker, first(starts))
        for _ in tracker
        end
        @test is_success(status(tracker))

        tracker, starts = pathtracker_startsolutions(
            f;
            projective_tracking = true,
            system = FPSystem,
        )
        @test tracker isa PathTracker{PVector{ComplexF64,1}}
        for x in starts
            @test is_success(track!(tracker, x))
            s = solution(tracker)
            @test s isa Vector{ComplexF64}
            @test abs(s[1]) ≈ sqrt(2) atol = 1e-8
            @test sum(s) ≈ 1 atol = 1e-8
            @test isnothing(winding_number(tracker))
            @test !isnan(tracker.state.solution_cond)
            @test !isnan(tracker.state.solution_accuracy)
            @test !isnan(tracker.state.solution_residual)
        end

        tracker, starts = pathtracker_startsolutions(homogenize(f, z); system = FPSystem)
        @test tracker isa PathTracker{PVector{ComplexF64,1}}
        for x in starts
            @test is_success(track!(tracker, x))
            s = solution(tracker)
            @test s isa PVector{ComplexF64,1}
        end

        g = equations(katsura(5))
        tracker, starts = pathtracker_startsolutions(g; system = FPSystem)
        @test all(s -> is_success(track!(tracker, s)), starts)

        # start solution is also target solution
        @polyvar x
        tracker = pathtracker([x - 1]; system = FPSystem)
        @test is_success(track(tracker, [1]))
    end

    # Univariate singular
    @testset "Univariate singular" begin
        @polyvar x
        print("Univariate singular: ")
        for n = 2:12
            print(n, ",")
            f = [(x - 3)^n]
            g = [x^n - 1]
            S = [[cis(i * 2π / n)] for i = 0:(n-1)]

            tracker = pathtracker(g, f, S; seed = 842121, system = FPSystem)
            @test all(S) do s
                is_success(track!(tracker, s)) &&
                winding_number(tracker) == n &&
                isapprox(solution(tracker)[1], 3.0; atol = 1e-6) &&
                !isnan(tracker.state.solution_cond) &&
                tracker.state.solution_accuracy < 1e-5
            end
        end
        println("")
    end

    # Wilkinson
    @testset "Wilkinson" begin
        @polyvar x
        print("Wilkinson: ")
        for n = 1:18
            print(n, ",")
            g = [x^n - 1]
            f = [prod(x - i for i = 1:n)]
            S = [[cis(i * 2π / n)] for i = 0:(n-1)]
            tracker = pathtracker(g, f, S; seed = 512321, system = FPSystem)
            @test all(s -> is_success(track!(tracker, s)), S)
            # check that we actually get the correct results
            @test sort(round.(
                Int,
                real.(first.(solution.(track.(tracker, S)))),
            )) == collect(1:n)
        end
        println("")
    end

    # Bivariate singular
    @testset "Bivariate singular" begin
        f = equations(griewank_osborne())
        tracker, starts = pathtracker_startsolutions(f; seed = 78373, system = FPSystem)
        @test all(starts) do s
            retcode = track!(tracker, s)
            (is_success(retcode) && winding_number(tracker) == 3) || is_at_infinity(retcode)
        end

        # Bivariate projetive singular
        @polyvar x z y
        # This has two roots of multiplicity 6 at the hyperplane z=0.
        # But the winding numbers are only 3 at each singularity
        F = [
            0.75 * x^4 + 1.5 * x^2 * y^2 - 2.5 * x^2 * z^2 + 0.75 * y^4 - 2.5 * y^2 * z^2 +
            0.75 * z^4,
            10 * x^2 * z + 10 * y^2 * z - 6 * z^3,
        ]
        @test pathtracker(F) isa PathTracker{PVector{ComplexF64,1}}
        tracker, starts = pathtracker_startsolutions(F; seed = 29831, system = FPSystem)
        @test all(s -> is_success(track!(tracker, s)), starts)
        @test all(starts) do s
            is_success(track!(tracker, s)) &&
            winding_number(tracker) == 3 && tracker.state.solution_cond > 1e10
        end
        # affine
        F̂ = subs.(F, Ref(y => 1))
        tracker, starts = pathtracker_startsolutions(F̂; seed = 29831, system = FPSystem)
        @test all(starts) do s
            is_success(track!(tracker, s)) &&
            winding_number(tracker) == 3 && tracker.state.solution_cond > 1e10
        end

        # Fix z == 1 --> all solutions at infinity
        G = subs.(F, Ref(z => 1))
        tracker, starts = pathtracker_startsolutions(G; seed = 29831, system = FPSystem)
        @test all(starts) do s
            is_at_infinity(track!(tracker, s))
        end
    end

    @testset "Non-Singular bad-conditioned" begin
        f = equations(PolynomialTestSystems.bacillus_subtilis())
        tracker, starts = pathtracker_startsolutions(f; seed = 78373, system = FPSystem)
        @test count(s -> is_success(track!(tracker, s)), starts) == 44
    end

    @testset "Diverging path without ill-conditioning" begin
        f = equations(cyclic(7))
        tracker, starts = pathtracker_startsolutions(f; seed = 78373, system = FPSystem)
        S = collect(starts)
        @test is_at_infinity(track!(tracker, S[7]))
        @test is_at_infinity(track!(tracker, S[9]))
        @test is_at_infinity(track!(tracker, S[5016]))
    end

    @testset "Cyclic 7 correct root count" begin
        f = equations(cyclic(7))
        tracker, starts = pathtracker_startsolutions(f; seed = 78373, system = FPSystem)
        results = map(s -> track!(tracker, s), starts)
        @test count(is_success, results) == 924

        successfull_paths = [1 4 6 10 13 14 16 18 43 44 46 59 63 66 69 70 73 78 85 93 94 98 104 106 111 119 121 122 134 141 148 152 155 157 160 163 171 173 182 185 189 200 204 206 208 211 215 218 225 227 230 231 236 244 245 246 256 258 259 261 264 272 277 280 284 287 298 305 318 319 326 336 340 342 347 358 367 380 384 397 400 405 410 411 416 420 423 424 430 431 439 442 447 470 476 483 486 497 510 511 550 552 558 566 569 575 583 585 587 589 591 599 601 602 605 610 615 622 627 630 634 638 639 645 650 651 653 656 658 669 673 679 681 683 689 691 693 696 707 709 715 719 722 723 724 733 746 755 757 765 767 770 773 786 796 805 807 808 811 816 818 832 837 840 846 850 857 858 864 873 874 875 880 882 900 905 909 910 911 918 919 920 923 925 929 935 936 942 943 947 955 960 972 979 982 990 1008 1010 1022 1031 1035 1037 1047 1059 1060 1062 1070 1075 1076 1082 1092 1096 1100 1102 1103 1109 1111 1115 1121 1124 1129 1147 1148 1152 1166 1168 1178 1189 1197 1198 1200 1201 1215 1216 1217 1223 1229 1235 1237 1242 1247 1251 1258 1261 1262 1272 1273 1277 1282 1293 1313 1318 1333 1336 1352 1362 1376 1377 1406 1408 1409 1419 1421 1422 1423 1424 1425 1442 1444 1446 1453 1454 1460 1462 1465 1473 1490 1503 1509 1520 1522 1523 1526 1527 1530 1531 1534 1549 1552 1558 1559 1561 1565 1579 1585 1594 1600 1621 1629 1631 1637 1640 1642 1644 1645 1649 1655 1656 1662 1665 1669 1671 1679 1705 1706 1712 1731 1734 1736 1737 1749 1750 1751 1753 1757 1759 1769 1784 1790 1795 1799 1800 1801 1802 1803 1807 1809 1817 1819 1825 1832 1833 1842 1847 1848 1852 1853 1856 1861 1865 1867 1871 1876 1885 1890 1896 1900 1903 1914 1925 1926 1927 1930 1933 1937 1947 1951 1954 1956 1962 1968 1971 1974 1980 1982 1988 1997 2000 2001 2014 2015 2017 2020 2030 2033 2045 2047 2054 2058 2062 2064 2066 2070 2080 2085 2099 2102 2103 2104 2111 2113 2122 2133 2136 2137 2139 2140 2145 2148 2149 2156 2159 2161 2163 2166 2171 2174 2182 2185 2190 2193 2195 2201 2202 2224 2225 2226 2229 2230 2233 2243 2260 2262 2267 2268 2272 2273 2292 2293 2295 2300 2301 2302 2306 2310 2314 2321 2326 2334 2338 2345 2356 2357 2358 2362 2363 2368 2369 2388 2391 2400 2405 2407 2416 2426 2430 2434 2435 2451 2454 2456 2463 2465 2473 2476 2477 2493 2495 2506 2509 2521 2527 2528 2551 2560 2562 2564 2566 2572 2574 2576 2603 2608 2613 2615 2624 2633 2635 2645 2650 2653 2658 2660 2668 2675 2682 2688 2692 2693 2696 2697 2699 2740 2741 2742 2767 2768 2777 2778 2784 2785 2799 2806 2808 2818 2822 2825 2836 2840 2852 2865 2875 2882 2894 2900 2902 2909 2910 2915 2925 2927 2935 2937 2938 2953 2956 2958 2974 2982 2985 2987 2991 2996 3000 3012 3013 3021 3024 3025 3026 3030 3032 3033 3042 3046 3064 3068 3077 3085 3091 3092 3093 3096 3102 3107 3110 3114 3119 3128 3135 3148 3150 3151 3154 3161 3164 3167 3169 3183 3187 3188 3190 3194 3201 3204 3206 3207 3212 3224 3230 3237 3242 3245 3250 3256 3260 3262 3263 3268 3271 3275 3277 3281 3283 3285 3288 3294 3305 3311 3319 3321 3324 3327 3335 3341 3344 3357 3367 3368 3371 3384 3388 3389 3397 3409 3412 3429 3434 3437 3444 3458 3464 3468 3473 3479 3492 3501 3505 3511 3513 3515 3517 3527 3539 3540 3549 3550 3553 3560 3565 3567 3577 3578 3584 3588 3600 3607 3608 3612 3614 3632 3659 3667 3670 3671 3674 3679 3680 3684 3686 3693 3700 3709 3712 3713 3729 3731 3735 3736 3737 3738 3749 3751 3752 3783 3798 3799 3804 3808 3815 3819 3826 3846 3847 3851 3855 3860 3864 3865 3871 3874 3880 3881 3886 3905 3915 3916 3917 3924 3941 3950 3956 3958 3963 3976 3978 3982 3986 3987 3989 3990 4002 4014 4015 4016 4018 4020 4023 4025 4032 4036 4045 4055 4064 4065 4068 4073 4078 4084 4087 4088 4097 4100 4107 4110 4116 4117 4135 4138 4143 4148 4150 4153 4159 4164 4168 4189 4193 4194 4199 4200 4201 4206 4208 4213 4219 4228 4238 4241 4243 4252 4254 4258 4259 4263 4264 4294 4303 4305 4314 4315 4324 4329 4359 4360 4363 4370 4379 4381 4383 4387 4390 4400 4403 4415 4436 4441 4452 4453 4461 4465 4466 4471 4478 4479 4480 4485 4504 4510 4526 4534 4550 4552 4553 4556 4566 4570 4576 4580 4584 4597 4600 4603 4606 4608 4613 4615 4616 4617 4618 4632 4639 4650 4663 4679 4683 4689 4691 4698 4700 4709 4718 4719 4720 4721 4724 4725 4728 4729 4734 4738 4739 4741 4744 4750 4753 4754 4759 4762 4768 4775 4776 4780 4781 4785 4788 4804 4807 4808 4814 4818 4821 4823 4826 4828 4830 4838 4843 4844 4864 4870 4872 4873 4889 4895 4909 4916 4918 4920 4922 4934 4936 4939 4949 4954 4962 4972 4973 4975 4976 4983 4989 4994 4996 4997 5003 5006 5008 5013 5015 5023 5025 5029 5035]

        @test all(successfull_paths' .== findall(is_success, results))
    end

    @testset "PathResult" begin
        # simple affine
        @polyvar x y
        f = [x^2 - 2, x + y - 1]
        tracker, starts = pathtracker_startsolutions(f; system = FPSystem)
        for (i, sᵢ) in enumerate(starts)
            result = track(tracker, sᵢ; path_number = i, details = :extensive)
            @test is_success(result)
            s = solution(tracker)
            @test abs(s[1]) ≈ sqrt(2) atol = 1e-8
            @test sum(s) ≈ 1 atol = 1e-8

            @test accuracy(result) < 1e-8
            @test residual(result) < 1e-8

            @test isnothing(winding_number(result))
            @test isnothing(multiplicity(result))
            @test condition_jacobian(result) < 1e3
            @test cond(result) == condition_jacobian(result)
            @test path_number(result) == i
            @test is_affine(result)
            @test !is_projective(result)
            @test is_finite(result)
            @test !is_singular(result)
            @test is_nonsingular(result)
            @test start_solution(result) == sᵢ
            @test is_nonsingular(result, 1e10)
            @test is_real(result)
            @test isreal(result)
            @test !is_failed(result)

            test_show_juno(result)
            @test !isempty(sprint(show, result))
        end

        # simple projective
        tracker, starts = pathtracker_startsolutions(homogenize(f); system = FPSystem)
        @test tracker isa PathTracker{PVector{ComplexF64,1}}
        for (i, sᵢ) in enumerate(starts)
            result = track(tracker, sᵢ; path_number = i, details = :extensive)
            @test is_success(result)
            @test is_projective(result)
            test_show_juno(result)
            @test !isempty(sprint(show, result))
        end


        # 1 singular solution (multiplicity 3) and 3 paths going to infinity. Checked earlier.
        f = equations(griewank_osborne())
        tracker, starts = pathtracker_startsolutions(f; seed = 78373, system = FPSystem)
        for sᵢ in starts
            result = track(tracker, sᵢ; details = :extensive)

            test_show_juno(result)
            @test !isempty(sprint(show, result))

            if is_at_infinity(result)
                @test result.t > 0
                @test !is_finite(result)
            elseif is_success(result)
                @test solution(result) ≈ [0.0, 0.0] atol = 1e-6
                @test is_real(result)
                @test winding_number(result) == 3
                @test accuracy(result) < 1e-6
                @test condition_jacobian(result) > 1e10
                @test is_singular(result)
                @test is_finite(result)
                @test isnothing(path_number(result))
            else
                @test false
            end
        end
    end

    @testset "Status codes" begin
        @test HC.path_tracker_status(HC.CoreTrackerStatus.success) == HC.PathTrackerStatus.success
        @test HC.path_tracker_status(HC.CoreTrackerStatus.terminated_invalid_startvalue) == HC.PathTrackerStatus.terminated_invalid_startvalue
        @test HC.path_tracker_status(HC.CoreTrackerStatus.terminated_maximal_iterations) == HC.PathTrackerStatus.terminated_max_iters
        @test HC.path_tracker_status(HC.CoreTrackerStatus.terminated_step_size_too_small) == HC.PathTrackerStatus.terminated_step_size_too_small
        @test HC.path_tracker_status(HC.CoreTrackerStatus.terminated_ill_conditioned) == HC.PathTrackerStatus.terminated_ill_conditioned
        @test HC.path_tracker_status(HC.CoreTrackerStatus.tracking) == HC.PathTrackerStatus.tracking
    end

    @testset "track with parameters change" begin
        @polyvar x a y b
        F = [x^2 - a, x * y - a + b]
        p = [a, b]

        tracker, starts = pathtracker_startsolutions(
            F,
            [1.0, 1.0 + 0.0 * im],
            parameters = p,
            p₁ = [5, 5], # wrong params
            p₀ = [10, 10],
            affine_tracking = true,
        )
        p₁ = [2, 0]
        p₀ = [2, 4]
        res = track(
            tracker,
            starts[1];
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )

        @test is_success(res)
        @test isa(solution(res), Vector{ComplexF64})
        @test length(solution(res)) == 2
    end

    @testset "Overdetermined" begin
        @polyvar x y z

        p₁ = (x^2 + y^2 + z^2 - 1) * (x - 0.5)
        p₂ = (x^2 + y^2 + z^2 - 1) * (y - 0.5)
        p₃ = (z - x^2 - 2) * (x^2 + y^2 + z^2 - 1) * (z - 0.5)
        L₁ = randn(ComplexF64, 2, 4) * [x, y, z, 1]
        L₂ = randn(ComplexF64, 2, 4) * [x, y, z, 1]

        s = let
            tracker, starts = pathtracker_startsolutions(
                [x^2 + y^2 + z^2 - 1; L₁];
                system = FPSystem,
            )
            solution(track(tracker, first(starts)))
        end

        tracker = pathtracker([[p₁, p₂, p₃]; L₁], [[p₁, p₂, p₃]; L₂], s; system = FPSystem)
        @test is_success(track(tracker, s))

        tracker, starts = pathtracker_startsolutions(
            [[p₁, p₂, p₃]; L₁],
            [[p₁, p₂, p₃]; L₂],
            s;
            system = FPSystem,
            projective_tracking = true,
        )
        @test is_success(track(tracker, s))
    end

    @testset "Homotopy input" begin
        @polyvar x a y b
        E = [[2 1 0; 0 0 0], [1 0; 1 0]]
        start = [[1.0 + 0im, -3.0, 2.0], [2.0 + 0im, -2.0]]
        target = [randn(ComplexF64, 3), randn(ComplexF64, 2)]
        H = CoefficientHomotopy(E, start, target)
        tracker = pathtracker(H, [[1, 1]])
        @test is_success(track!(tracker, [1, 1]))
    end
end
