@testset "Certification" begin
    @testset "Simple 1: Input = Vector{Expression}" begin
        @var x y
        # define the polynomials
        f₁ = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        f₂ = x^2 + 2x * y^2 - 2 * y^2 - 1 / 2
        F = [f₁, f₂]
        result = solve(F)

        # Input: System(F)
        cert = certify(System(F), result)
        @test nnonsingular(result) == 18
        @test ncertified(cert) == 18
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 4
        @test ncomplex_certified(cert) == 14
        @test ndistinct_real_certified(cert) == 4
        @test ndistinct_complex_certified(cert) == 14
        save("tmp_cert.txt", cert)
        @test !isempty(read("tmp_cert.txt", String))

        # Input: F
        cert = certify(F, result)
        @test nnonsingular(result) == 18
        @test ncertified(cert) == 18
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 4
        @test ndistinct_real_certified(cert) == 4

        # Control Display
        cert = certify(F, result, show_progress = false)

        # Double solutions
        S = solutions(result)
        cert = certify(F, [S; S]; extended_certificate = true)
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 4
    end

    @testset "Simple 2: Input = Vector{MP.Polynomial}" begin
        @polyvar x y
        # define the polynomials
        f₁ = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        f₂ = x^2 + 2x * y^2 - 2 * y^2 - 1 / 2
        F = [f₁, f₂]
        result = solve(F)

        # Input: System(F)
        cert = certify(System(F), result; extended_certificate = true)
        @test nnonsingular(result) == 18
        @test ncertified(cert) == 18
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 4
        @test ndistinct_real_certified(cert) == 4

        # Input: F
        cert = certify(F, result)
        @test nnonsingular(result) == 18
        @test ncertified(cert) == 18
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 4
        @test ndistinct_real_certified(cert) == 4

        # Control Display
        cert = certify(F, result, show_progress = false)

        # Double solutions
        S = solutions(result)
        cert = certify(System(F), [S; S])
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 18
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 4
    end

    @testset "Parameters: Input = Vector{Expression}" begin
        @var x y
        f = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        # define new variables u₁, u₂ and λ₁
        @var λ[1:1] u[1:2]
        # define the jacobian of F
        J = differentiate([f], [x, y])
        # J' defines the transpose of J
        F = [[x, y] - u - J' * λ; f]
        u₀ = [-0.32, -0.1]
        res = solve(F, parameters = u, target_parameters = u₀)

        # Input: System(F)
        C = System(F, parameters = u)
        cert = certify(C, res; target_parameters = u₀)
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 36
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 8

        # Input: F
        cert = certify(F, res, u₀; parameters = u)
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 36
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 8

        # Invalid solutions
        cert = certify(C, [100 .* randn(ComplexF64, 3) for i = 1:10], u₀)
        @test ncertified(cert) < 10
    end

    @testset "Parameters: Input = Vector{MP.Polynomial}" begin
        @polyvar x y
        f = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        # define new variables u₁, u₂ and λ₁
        @polyvar λ[1:1] u[1:2]
        # define the jacobian of F
        J = differentiate([f], [x, y])
        # J' defines the transpose of J
        F = [[x, y] - u - J' * λ; f]
        u₀ = [-0.32, -0.1]
        res = solve(F, parameters = u, target_parameters = u₀)

        # Input: System(F)
        C = System(F, parameters = u)
        cert = certify(C, res; target_parameters = u₀)
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 36
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 8

        # Input: F
        cert = certify(F, res; parameters = u, target_parameters = u₀)
        @test ncertified(cert) == 36
        @test ndistinct_certified(cert) == 36
        @test nreal_certified(cert) == 8
        @test ndistinct_real_certified(cert) == 8
    end

    @testset "positive" begin
        @var x y
        f = System([x^2 + y^2 - 1, x - y])
        res = solve(f; compile = false, start_system = :total_degree)
        cert = certify(f, res, compile = false)
        @test count(is_positive, certificates(cert)) == 1
        @test count(s -> is_positive(s, 1), certificates(cert)) == 1
        @test count(is_real, certificates(cert)) == 2
    end

    @testset "3264" begin
        F = steiner()
        real_conics = [
            10124547 // 662488724,
            8554609 // 755781377,
            5860508 // 2798943247,
            -251402893 // 1016797750,
            -25443962 // 277938473,
            1 // 1,
            520811 // 1788018449,
            2183697 // 542440933,
            9030222 // 652429049,
            -12680955 // 370629407,
            -24872323 // 105706890,
            1 // 1,
            6537193 // 241535591,
            -7424602 // 363844915,
            6264373 // 1630169777,
            13097677 // 39806827,
            -29825861 // 240478169,
            1 // 1,
            13173269 // 2284890206,
            4510030 // 483147459,
            2224435 // 588965799,
            33318719 // 219393000,
            92891037 // 755709662,
            1 // 1,
            8275097 // 452566634,
            -19174153 // 408565940,
            5184916 // 172253855,
            -23713234 // 87670601,
            28246737 // 81404569,
            1 // 1,
        ]
        real_sols = read_solutions(joinpath(@__DIR__, "data/3264_real_sols.txt"))
        cert = certify(F, real_sols, real_conics; compile = true)
        @test ndistinct_real_certified(cert) == 3264
        cert = certify(F, real_sols, real_conics; compile = false)
        @test ndistinct_real_certified(cert) == 3264
        @test_throws ArgumentError certify(F, real_sols; compile = false)
        @test_throws ArgumentError certify(F, real_sols; compile = true)

        dcs = DistinctCertifiedSolutions(F, real_conics)
        for s in real_sols
            add_solution!(dcs, s, 1)
        end
        @test length(solutions(dcs)) == 3264

        dcs2 = distinct_certified_solutions(
            F,
            [real_sols; real_sols[1:100]],
            real_conics,
            threading = true,
            show_progress = false,
        )
        @test length(solutions(dcs2)) == 3264
    end

    @testset "certify uses complex inversion" begin
        @var x y
        F = System([y / x + x^2 - 3], parameters = [y])
        monres = monodromy_solve(F)
        cert = certify(F, monres)
        @test ndistinct_certified(cert) >= 1
    end

    # This solution is fairly bad conditioned and requires that the approx. inverse of the jacobian 
    # becomes more accurate over time. Added this test as a result of the bug where we just used inv!
    # to compute C but this had too large of an error bound to certify the solution.
    @testset "certify uses approximate inverse of jacobian" begin
        @var x[1:8]

        M = [1 0 0 1 1 1 1 1; 0 1 0 1 x[1] x[2] x[3] x[4]; 0 0 1 1 x[5] x[6] x[7] x[8]]

        n = size(M, 2)

        S = collect(Combinatorics.combinations(collect(1:n), 3))
        nonconstant_minors = filter(s -> !ModelKit.is_number(s), [det(M[:, s]) for s in S])

        V = zeros(Expression, 6, 8)
        for i = 1:8
            c = 1
            for j = 1:3
                for k = j:3
                    V[c, i] = M[j, i] * M[k, i]
                    c = c + 1
                end
            end
        end
        conics =
            [det(V[:, s]) for s in Combinatorics.combinations(collect(1:size(V, 2)), 6)]


        @var s[1:length(nonconstant_minors)+length(conics)]

        expressions = HomotopyContinuation.differentiate(
            sum(s .* log.([nonconstant_minors; conics])),
            x,
        )
        F = System(expressions; parameters = s)

        p = ComplexF64[
            0.04694929498353116-0.9795943487338107im,
            0.03814936256080075+0.6036805123359278im,
            0.19811294615170094-0.04655749375473129im,
            -0.568612454485359-0.30035127049952426im,
            -0.15638090220760134-1.2654555749207392im,
            -0.41345309978217826-0.1784389293895055im,
            -1.1113178721109995-0.22606998166881093im,
            0.9982708200513642+0.8528947638503744im,
            0.20166791921434227+0.9988764609171762im,
            -0.8100994025378792+0.3969273031583972im,
            0.6071086706963922-1.4311940916205672im,
            -0.01611865898465837-0.4502118183049486im,
            -0.110443619721666-0.7045834618227158im,
            0.7234274346282669+0.6283366831106044im,
            0.36809026072491596-0.8203266151063194im,
            1.1813346349972564+0.3998055846562327im,
            0.7027850166072859-1.2954603858129299im,
            1.0399846079963362-0.08005516566457833im,
            0.9198982114121447+0.9781446557067801im,
            -0.6891065881783133-0.4053802660764455im,
            0.483503689861307-0.3225817056899272im,
            -1.3251128587644545+0.897858125440237im,
            0.0747608828119116-0.5836512862131934im,
            -0.17485527683038962+0.44459749365390716im,
            0.21960852538572917-0.5379498288177373im,
            -0.8014331739697617+0.059816621817923445im,
            -0.43460375391295053-0.2546402146687261im,
            -0.6046592521930908-0.2957968257500265im,
            -1.1053855142695876+1.590614119832034im,
            0.9840486357930686+0.7702769853904382im,
            0.16170642515220965-0.1222743301057992im,
            0.4795973979603369+0.2810297596371718im,
            -1.2000892295608054-0.8923206318916687im,
            -0.6077910779725929+0.26436888705825917im,
            -1.0376531919986127-0.33652306556055805im,
            -0.06733756067158703-0.639783597527494im,
            -0.37371460434888953+0.24901918161963316im,
            0.017508442680565478+0.32123450723164876im,
            1.4128224393512774+0.33400143986123293im,
            -0.0434837081190993-0.8985222271337256im,
            0.06449022830922942+0.04671222198604458im,
            0.24633610502585773+0.5801362854382477im,
            -0.6955336053945378+1.0796258697355416im,
            -0.6329313759525533-0.2354860006860427im,
            0.6109430855309662-0.6729029511956512im,
            -0.260922761581385-0.12434880231119286im,
            0.059517746441893456+0.38136326366327844im,
            0.12550005778394707+0.4949701677886087im,
            0.21111192677412854+0.6810392555780735im,
            0.6320528544930368-0.2833511182719277im,
            -0.6050405544591456-0.9440915895571896im,
            -1.4345871717374217+0.22556237643401553im,
            -0.6887743690027825-0.47027128695679876im,
            -0.6415806120384493-0.5014670765374285im,
            0.08502741985519593-0.6235217666945221im,
            0.21493766396179834+0.8552649062604912im,
            0.3046927960034058-1.590707502314377im,
            0.44767647377604874+0.8188582740588161im,
            1.1192822043579032+0.3234506727151911im,
            -0.008792767666757002-0.21480481868217485im,
            -1.1022743248740579+1.439084779398887im,
            -1.25202306664751-0.4856621088482818im,
            0.5609474151450039-0.06309637623414076im,
            -0.20050911884693357+0.2973185023874754im,
            -0.09243154166047889-0.49984234324259064im,
            -0.28169898290001166+0.9680337642860322im,
            1.8202758158992145-0.7826768803555774im,
            0.46305418751090566+0.19925413382595714im,
            -0.017534482018503578-0.6315561160225208im,
            0.28487129868590444+1.299813316065981im,
            -0.507158318173951-0.7613059977835542im,
            -1.5246336574768562+0.15842771871857747im,
            0.19964798257040645-0.028394263163747824im,
            -1.3150248862774019+0.43865160271978626im,
            -0.1447364027631027+0.04047743522664874im,
            0.18766444754121442+0.8548396241105418im,
        ]
        sol = ComplexF64[
            -267.34389183407836+90.16013119150394im,
            -0.5175728762869246-0.4762083613875728im,
            -0.5275621759154175-0.48081770228311727im,
            -0.5187465849951617-0.49134401394255406im,
            -0.14084183600441835+1.181750659496792im,
            -0.4215223684209179+0.1331890129058157im,
            -0.42154412086687626+0.13324039516624503im,
            -0.42153938328192286+0.13323037611055516im,
        ]
        f = fixed(F)
        c = HomotopyContinuation.CertificationCache(f)
        q = HomotopyContinuation.certification_parameters(p)
        cert = HomotopyContinuation.certify_solution(f, sol, q, c, 1)
        @test is_certified(cert)
    end
end
