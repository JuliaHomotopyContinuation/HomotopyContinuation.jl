@testset "Utilities" begin

    @testset "Group actions" begin
        f1 = s -> (s * s,);
        f2 = s-> (2s, -s, 5s);
        f3 = s -> (s + 1,);
        action1 = GroupActions(f1)
        action2 = GroupActions(f1, f2)
        action3 = GroupActions(f1, f2, f3)
        r1 = Int[]
        HC.apply_actions(action1, 3) do s
            push!(r1, s)
            false
        end
        @test r1 == [9]

        r2 = Int[]
        HC.apply_actions(action2, 3) do s
            push!(r2, s)
            false
        end
        r2
        @test r2 == [9, 18, -9, 45, 6, -3, 15]

        r3 = Int[]
        HC.apply_actions(action3, 3) do s
            push!(r3, s)
            false
        end
        r3
        @test r3 == [9, 18, 19, -9, -8, 45, 46, 10, 6, 7, -3, -2, 15, 16, 4]
        @test action3(3) == [3, 9, 18, 19, -9, -8, 45, 46, 10, 6, 7, -3, -2, 15, 16, 4]

        # also test with arrays
        g1 = s -> [s * s];
        g2 = s-> [2s, -s, 5s];
        g3 = s -> [s + 1];
        action1 = GroupActions(g1)
        action2 = GroupActions(g1, g2)
        action3 = GroupActions(g1, g2, g3)
        r1 = Int[]
        HC.apply_actions(action1, 3) do s
            push!(r1, s)
            false
        end
        @test r1 == [9]

        r2 = Int[]
        HC.apply_actions(action2, 3) do s
            push!(r2, s)
            false
        end
        @test r2 == [9, 18, -9, 45, 6, -3, 15]

        r3 = Int[]
        HC.apply_actions(action3, 3) do s
            push!(r3, s)
            false
        end
        r3
        @test r3 == [9, 18, 19, -9, -8, 45, 46, 10, 6, 7, -3, -2, 15, 16, 4]
    end

    @testset "UniquePoints" begin
        Random.seed!(1234)
        X = [randn(ComplexF64, 10) for _ = 1:2_000]
        indices = HC.unique!(rand(1:2000, 20))
        data = HC.UniquePoints(X)

        test_show_juno(data)
        @test string(data) == "UniquePoints{Array{Complex{Float64},1},Float64,typeof(euclidean_distance),Nothing} with 2000 points"

        @test length(data) == 2_000

        for i ∈ indices
            @test HC.iscontained(data, X[i])
            @test HC.iscontained(data, X[i], Val(true)) == i
            @test HC.iscontained(data, X[i] .+ 1e-4) == false
            @test HC.iscontained(data, X[i] .+ 1e-9, Val(true)) == i
            @test HC.iscontained(data, X[i] .+ 1e-9) == true
            @test HC.add!(data, X[i]) == false
            @test HC.add!(data, X[i], Val(true)) == i
            @test data[i] == X[i]
            @test HC.add!(data, X[i] .+ 1e-4) == true
            @test HC.add!(data, X[i] .- 1e-4, Val(true)) == -1
        end

        empty!(data)
        @test length(points(data)) == 0

        # Test many points with nearly indentical distance to the inserted point
        p = shuffle!([[cis(k/100*2π)] for k=0:99])
        data = HC.UniquePoints(p)
        @test HC.iscontained(data, [0.0im]) == false

        # Test with group action
        x = randn(ComplexF64, 4)
        permutation1(x) = ([x[2]; x[1]; x[3]; x[4]],)
        permutation2(x) = ([x[1]; x[2]; x[4]; x[3]],)
        X = [v for v in GroupActions(permutation1, permutation2)(x)]

        # One group action
        data = UniquePoints(X, group_action = permutation1)
        @test length(data) == 2

        # Two group actions
        data = HC.UniquePoints(X, group_actions = GroupActions(permutation1, permutation2))
        @test length(data) == 1

        # Group action and reality check
        x = randn(4)
        X = [im.*x, x, randn(ComplexF64, 4)]

        data = HC.UniquePoints(X, group_action = x -> (im.*x, (-1).*x, (-im).*x), check_real=false)
        @test HC.isrealvector(points(data)[1]) == false
        @test length(points(data)) == 2

        data = HC.UniquePoints(X, group_action = x -> (im.*x, (-1).*x, (-im).*x), check_real = true)
        @test HC.isrealvector(points(data)[1]) == true
        @test length(points(data)) == 2

        data = HC.UniquePoints(randn(ComplexF64, 4), check_real = true)
        k = add!(data, randn(4), Val{true}())
        @test k == -2
        k = add!(data, randn(ComplexF64, 4), Val{true}())
        @test k == -1
    end

    @testset "Multiplicities" begin
        V = [randn(3) + im.* randn(3) for i in 1:10]
        W = [map(v -> v + [1e-7 * v[1] ;1.0;1.0], V); V; map(v -> v + [1e-6 * v[1]; 0; 0], V)]
        U = [[cis(rand()) for _=1:3] for _=1:20]
        X = [[0.499064-0.707804im, -0.164223+0.0933025im, -0.00695332-0.139516im, 0.000716617-0.135055im, 0.009846-0.12984im, 0.020204-0.123609im, 0.34326+0.161402im], [-0.795428-0.0975931im, -0.20942-0.382162im, -0.119855-0.0414809im, -0.123105-0.0293171im, -0.128161-0.0148951im, -0.134641+0.00272884im, 0.258239+0.185188im], [0.605937-0.289401im, 0.42308-0.188036im, -0.057496+0.200812im,0.284444+0.0477843im, 0.191189-0.225824im, -0.233732-0.181827im, -0.0939606+0.154403im], [-0.493372+0.512561im, -0.385166+0.204525im, -0.120853-0.142237im, 0.0742939+0.163197im, -0.298456+0.202635im, 0.100108-0.219369im, 0.193075+0.105514im], [0.483174-0.581023im, 0.23432-0.0528133im, -0.111264-0.067633im, -0.292867-0.014453im, 0.24805-0.267034im, 0.141232-0.13631im, 0.168603+0.261794im], [-0.455116+0.4482im, -0.366004+0.264082im, -0.0578486+0.110708im, -0.296085-0.0619238im, 0.216907-0.327015im, 0.207934+0.242906im, -0.121251+0.101268im], [0.526416-0.465825im, 0.39214-0.0908531im, -0.260222+0.00391701im, -0.110348-0.1151im, 0.101057-0.190265im, 0.308379-0.137148im, 0.0488701+0.29665im], [0.147886-0.750864im, -0.235168-0.257487im, 0.0394847+0.194512im, -0.0325228-0.227443im, -0.111526-0.214717im, 0.0245065+0.162303im, 0.304189-0.150105im], [-0.372175-0.538886im, -0.117563-0.341994im, 0.0752294-0.0141672im, 0.143671-0.249302im, -0.374165-0.228281im, -0.0971581+0.346633im, 0.168499+0.0394346im], [-0.503759+0.48im, -0.447123+0.262057im, -0.214529+0.109917im, 0.148493-0.0330798im, -0.0510215-0.262778im, -0.110165+0.198317im, 0.200704+0.0508037im], [0.623685-0.513102im, 0.0522679+0.0223332im, -0.188478-0.10101im, 0.229926-0.157889im, -0.146424-0.0562049im, 0.209112-0.185147im, 0.205229+0.276164im], [0.127017-0.655061im, -0.147053-0.569588im, -0.134863+0.0901295im, 0.152888+0.049255im, -0.124224+0.0624774im, 0.152924+0.0778815im, 0.105575-0.310906im], [0.575711+0.384346im, 0.336129+0.0531313im,0.155698-0.152359im, -0.104555+0.197167im, -0.254208-0.215823im, 0.327345+0.0728999im, -0.0274605+0.288691im], [0.479139-0.593553im, -0.0960412-0.315975im, 0.0438803-0.245038im, 0.000516985+0.193173im, -0.0828034+0.172297im, -0.00479366-0.21881im, 0.353615-0.0168729im], [-0.37103+0.609109im, -0.201984+0.259807im, -0.158043+0.252162im, -0.211405-0.106454im, 0.0274292-0.262351im, 0.333102+0.0455649im, -0.0847849+0.220536im], [-0.200697+0.682536im, -0.247073+0.359358im, -0.172902-0.0702972im, 0.141297+0.110399im, -0.171049+0.317616im, -0.0125869-0.240803im, 0.219982+0.00435721im], [-0.0896611-0.805671im, -0.0393637-0.183565im, 0.224508-0.167715im, 0.0918438-0.180778im, -0.0879493-0.157071im, -0.250607-0.0700599im, 0.0521971+0.291749im], [-0.226054+0.605127im, -0.210133+0.55141im,0.218162-0.0576634im, 0.08772-0.0936im, -0.0820928-0.101495im, -0.22419-0.0507194im, -0.0146449+0.311527im], [0.498516+0.360158im, 0.219394+0.0657585im, -0.0836376+0.0106629im, 0.00702255-0.322935im, 0.481634+0.13367im, -0.133322+0.397338im, -0.171875-0.0532251im], [0.440673-0.62568im, -0.110019-0.330907im, -0.043096+0.193744im, 0.0628576-0.220991im, -0.0144861-0.24152im, -0.0436694+0.158228im, 0.338948-0.0133119im], [0.622749+0.219955im, 0.628838+0.216896im, 0.0331851+0.145177im, -0.123675-0.0636893im, 0.0348129+0.117408im, -0.133898-0.0420959im, 0.155277-0.144857im], [0.652186-0.0569317im, 0.552174-0.169926im, -0.0977967-0.0521929im, -0.161552-0.169677im, 0.00437006+0.237406im, -0.0164646+0.116656im, 0.310812-0.0599935im], [0.448293+0.574542im, 0.247779+0.106134im, -0.0955677-0.215494im, 0.0935008+0.252581im, 0.336093-0.0655203im, -0.242972-0.0256941im, -0.0177608+0.301502im], [0.445422+0.671948im, -0.0726023+0.313409im, -0.044775+0.21726im,0.154346+0.0261763im, -0.0305679+0.188435im, 0.20006+0.0391233im, -0.0810809-0.297159im], [-0.00979382-0.81029im, 0.106391-0.195197im, 0.108622-0.139201im, 0.234362-0.1806im, -0.20136-0.131455im, -0.11823-0.118552im, -0.0308846+0.297292im], [0.318611+0.584183im, 0.105547+0.221636im, 0.248174+0.0583734im, -0.00256072-0.290164im, -0.302039-0.0824955im, -0.114841+0.409322im, 0.25213+0.0736154im], [0.696901+0.0684572im, 0.505426+0.266373im, -0.152493+0.0711982im, 0.197532-0.0412531im, 0.214441+0.0165324im, -0.096866+0.0847148im, -0.103987-0.201388im], [0.115188-0.657244im, -0.157297-0.566844im, -0.133216+0.092546im, 0.153751+0.0464909im, -0.123077+0.0647066im, 0.154303+0.0751121im, 0.0999536-0.312759im], [0.653199-0.277909im, 0.556407-0.00960311im, 0.204776-0.112518im, -0.112199+0.0771487im, -0.0695128+0.163961im, 0.167741-0.0723464im, -0.16148-0.148729im], [0.411538+0.550657im, 0.413916+0.420472im, -0.0555052+0.089203im, -0.0609812+0.227821im, -0.0913546-0.187807im, -0.0607857-0.0922333im, 0.237526+0.0198403im], [0.602049+0.184366im, 0.695878+0.270759im, -0.0420156+0.0359236im, -0.0431941+0.0331572im, -0.0446403+0.0302374im, -0.046326+0.027106im, 0.0979775-0.156785im], [-0.270799+0.567571im, -0.113714+0.377218im, 0.0592799-0.219174im, 0.226587+0.091218im, 0.0269123+0.326853im, -0.402766+0.100383im, -0.018186-0.240642im]]

        M = HC.multiplicities(V, distance=HC.infinity_norm)
        @test length(M) == 0

        N = HC.multiplicities(W, distance=HC.infinity_norm)
        sort!(N, by=first)
        @test length(N) == 10
        @test unique([length(m) for m in N]) == [2]
        @test N == [[i, i+10] for i=11:20]

        O = HC.multiplicities([U;U], distance=HC.infinity_norm)
        @test length(O) == 20

        P = HC.multiplicities(X, distance=(x,y) -> 1-abs(LinearAlgebra.dot(x,y)))
        @test length(P) == 3

        # Test with group action
        x = randn(ComplexF64, 4)
        permutation1(x) = ([x[2]; x[1]; x[3]; x[4]],)
        permutation2(x) = ([x[1]; x[2]; x[4]; x[3]],)
        X = GroupActions(permutation1, permutation2)(x)

        # One group action
        m = multiplicities(X, group_action = permutation1) |> sort
        @test m == [[1, 2], [3, 4]]

        # Two group actions
        m = multiplicities(X, group_actions = GroupActions(permutation1, permutation2))
        @test m == [[1, 2, 3, 4]]
    end

    @testset "Polynomials" begin
        @polyvar p q a b c x y z
        e = [p + 1]
        f = [a * b * c ]
        g = [x+y, y + z, x + z]
        @test expand(f ∘ g) == [x^2*y + x^2*z + x*y^2 + 2*x*y*z + x*z^2 + y^2*z + y*z^2]
        @test expand(e ∘ f ∘ g) == [x^2*y + x^2*z + x*y^2 + 2*x*y*z + x*z^2 + y^2*z + y*z^2 + 1]

        @test validate(f ∘ g) == true
        @test validate(f ∘ g, parameters=[z]) == true
        @test validate(f ∘ g, parameters=[c]) == false
        @test validate(g ∘ f) == false

        @test HC.nvariables(f ∘ g) == 3
        @test HC.nvariables(f ∘ g, parameters=[z]) == 2

        @polyvar x y z
        @test ishomogeneous([x^2+y^2+x*y, x^5])
        @test ishomogeneous([x^2+y^2+x*y, x^4+1]) == false

        #test weighted degree
        @test ishomogeneous(x^3+x*y, [(x, 2), (y, 4)])
        @test homogenize(x+x*y, [(x, 2), (y, 4)], z) == x*z^4+x*y

        @test ishomogeneous(homogenize([x^2+y^2+x*y, x^4+1]))
        @test ishomogeneous([x^2+z^2 + y, x^4+z^4], [x,z]) == false
        @test ishomogeneous(homogenize([x^2+z^2*y, x^4+z^4*y], [x,z]), [x,z]) == true

        @test ishomogeneous(f ∘ g) == true
        h = [a * b * c  + p * c^3] ∘ [x+y, y + z, x + z]
        @test ishomogeneous(h, parameters=[p]) == true
        h2 = [a * b * c  + p] ∘ [x+y, y + z, x + z]
        @test ishomogeneous(h2, parameters=[p]) == false


        #homogenize
        @test homogenize(x^2+y+1, z) == x^2+y*z + z^2
        # This needs to be an array due to a compiler bug
        @test homogenize([x^2+y+p], z, parameters=[p]) == [x^2+y*z + z^2*p]

        @test homogenize([x^3+p, 1+y], z, parameters=[p]) == [x^3+p*z^3, z+y]

        h2 = [a * b * c  + p] ∘ [x+y, y + z, x + q]
        @test validate(homogenize(h2, parameters=[p, q]), parameters=[p, q])
        h3 = [a * b * c  + 1] ∘ [x+y, y + z, x + 1]
        @test validate(homogenize(h3))

        # Weylnorm
        @polyvar x y z
        f = 3.0x^2 + 2x*y - y^2
        g = (-2.5+2im) * x^2 - 3.0*x*y + 4y^2
        @test HC.weyldot(f, g) == 3.0 * conj(-2.5 + 2im) + 2.0 * (-3.0) / 2 + (-1.0) * 4.0
        @test HC.weyldot(f, f) == 9.0 + 4.0  / 2 + 1.0
        @test HC.weylnorm(f)^2 ≈ HC.weyldot(f,f)
        @test HC.weyldot([f, f], [g, g]) == 2 * HC.weyldot(f, g)
        @test HC.weylnorm([f, f]) == √HC.weyldot([f, f], [f, f])
    end

    @testset "Multihomogeneous" begin
        @polyvar x y
        f = [x*y-2, x^2-4]

        affine_hominfo = HC.HomogenizationInformation(variable_groups = ((x,), (y,)))
        @test HC.ishomogeneous(f, affine_hominfo) == false

        @polyvar v w
        g = [x * y - 2v * w, x^2-4v^2]
        hominfo = HC.HomogenizationInformation(variable_groups = ((x,v), (y,w)))
        @test HC.ishomogeneous(g, hominfo) == true

        hominfo = HC.HomogenizationInformation(homvars=(v,w), variable_groups = ((x,v), (y,w)))
        @test HC.ishomogeneous(g, hominfo) == true

        @test HC.homogenize(f, hominfo) == g
        HC.multidegrees(f, ([x], [y])) == [1 2; 1 0]
    end

    @testset "Misc" begin
        A = rand(Complex{Float64}, 12, 12)
        b = rand(Complex{Float64}, 12)
        C, d = copy(A), copy(b)
        @test norm(HC.solve!(C, d) - A \ b) < 1e-10

        A = rand(Complex{Float64}, 15, 12)
        b = rand(Complex{Float64}, 15)
        C, d = copy(A), copy(b)
        HC.solve!(C, d)
        @test norm(d[1:12] - A \ b) < 1e-10

        @test HC.unpack(5, 2) == 5
        @test HC.unpack(nothing, 2) == 2

        x = rand()
        @test HC.nthroot(x, 1) == x
        @test HC.nthroot(x, 0) == one(x)

        x = rand(6)
        @test HC.isrealvector(x)

        y = x + [0, 0, 2, 0, 0, 0]
        @test infinity_distance(x, y) ≈ 2
        @test infinity_distance(complex.(x), complex.(y)) ≈ 2


        segment = HC.ComplexSegment(2, 4)
        @test length(segment) ≈ 2
        test_show_juno(segment)
        @test string(segment) == "ComplexSegment(2.0 + 0.0im, 4.0 + 0.0im)"
    end
end
