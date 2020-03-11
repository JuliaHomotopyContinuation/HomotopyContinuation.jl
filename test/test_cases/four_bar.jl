@testset "Iterative Refinement for predictor (FourBar)" begin
    @var x a y b x_hat a_hat y_hat b_hat
    gamma = [Variable(Symbol(:gamma, i)) for i = 1:8]
    delta = [Variable(Symbol(:delta, i)) for i = 1:8]
    gamma_hat = [Variable(Symbol(:gamma_hat, i)) for i = 1:8]
    delta_hat = [Variable(Symbol(:delta_hat, i)) for i = 1:8]
    #system of polynomials
    D1 = map(1:8) do i
        (a_hat * x - delta_hat[i] * x) * gamma[i] +
        (a * x_hat - delta[i] * x_hat) * gamma_hat[i] +
        (a_hat - x_hat) * delta[i] +
        (a - x) * delta_hat[i] - delta[i] * delta_hat[i]
    end
    D2 = map(1:8) do i
        (b_hat * y - delta_hat[i] * y) * gamma[i] +
        (b * y_hat - delta[i] * y_hat) * gamma_hat[i] +
        (b_hat - y_hat) * delta[i] +
        (b - y) * delta_hat[i] - delta[i] * delta_hat[i]
    end
    D3 = [gamma[i] * gamma_hat[i] + gamma[i] + gamma_hat[i] for i = 1:8]
    F = System(
        [D1; D2; D3],
        [x; a; y; b; x_hat; a_hat; y_hat; b_hat; gamma; gamma_hat],
        [delta; delta_hat],
    )

    s = Complex{Float64}[
        -2.6370887+0.35329306im,
        -3.7223306+2.5267087im,
        -0.70272787-0.74162847im,
        -0.52112188-1.5249im,
        -0.29444549+0.2698171im,
        0.23877446+1.1593521im,
        -0.2226372+0.24523739im,
        -0.23666974+0.87507397im,
        2.5277031-6.3353057im,
        -0.88010494+1.3133313im,
        -0.73606124+0.15544107im,
        0.80744564-3.0613903im,
        -0.82122921+0.80545901im,
        0.68534751-2.5000172im,
        0.65872999-2.3719828im,
        -0.87767826-0.35435132im,
        -0.93290889+0.12048708im,
        -0.93106365-0.75512927im,
        1.8130785-1.6567022im,
        -0.85699423+0.24221833im,
        -0.73738109-1.1832401im,
        -0.81460307+0.2750148im,
        -0.80200623+0.28313097im,
        -0.12955283+2.5215805im,
    ]
    p = Complex{Float64}[
        7.297148065259444-3.7969299471707454im,
        0.06125202773860465+0.9251492821389162im,
        1.4586512320807385+0.641527879044937im,
        3.3162292880980537-2.979942348702623im,
        -0.6613217458488018+0.8321586519581641im,
        2.6212360706406703-2.6766427578876666im,
        2.4522388848243555-2.5938408784246976im,
        -1.7265529790073868-0.29337325907998707im,
        0.24011427062306998+1.375412535460542im,
        -0.19336068638852033+0.5669754468190339im,
        0.147943349839693-0.2784406507316879im,
        -0.01641513182688631+1.5963735195205522im,
        -0.33880136827427765+0.2772103200265199im,
        -0.19788506286822416+1.6632361622942196im,
        -0.26097545384785903+1.6764055030939693im,
        0.44012684131549085+1.0212717350758722im,
    ]
    q = Complex{Float64}[
        -0.35263257737537096+1.2598825875492277im,
        -0.5353479512983561-0.5963278468670689im,
        1.4096839876427854-1.2227632666172732im,
        -0.2550909822371234-0.641270846350799im,
        0.24180730168363096-0.46611650665985527im,
        0.29953332004980354-0.8469219548644632im,
        0.35386303195558433-0.3683866291866925im,
        1.2021910424897007+1.1148794002014533im,
        -0.35263257737537096-1.2598825875492277im,
        -0.5353479512983561+0.5963278468670689im,
        1.4096839876427854+1.2227632666172732im,
        -0.2550909822371234+0.641270846350799im,
        0.24180730168363096+0.46611650665985527im,
        0.29953332004980354+0.8469219548644632im,
        0.35386303195558433+0.3683866291866925im,
        1.2021910424897007-1.1148794002014533im,
    ]

    tracker = Tracker(ParameterHomotopy(F, p, q))
    res = track(tracker, s, 1, 0)
    @test is_success(res)

    tracker.options.automatic_differentiation = (true, true, true, false)
    res = track(tracker, s, 1, 0)
    @test is_success(res)

    tracker.options.max_steps = 10_000
    tracker.options.automatic_differentiation = (true, true, false, false)
    res = track(tracker, s, 1, 0)
    @test is_success(res)
end
