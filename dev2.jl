using HomotopyContinuation
using HomotopyContinuation: ComplexDF64
import Arblib
import Arblib: Acb, Arb, AcbMatrix, AcbVector
using LinearAlgebra

using BenchmarkTools
# include("src/new_tracker/adaptive_tracker.jl")

include("test/test_systems.jl")
include("test/engame_test_systems.jl")

# f = eg_system_1(2)
# (T, xs) = total_degree(f)


# function estimate_system_distance(H)
#     m, n = size(H)
#     u1 = zeros(ComplexF64, m)
#     u0 = zeros(ComplexF64, m)
#     v = zeros(m)
#     normalizer = zeros(m)

#     nsamples = 5

#     for i = 1:nsamples
#         x = randn(ComplexF64, n)
#         x ./= norm(x)
#         evaluate!(u1, H, x, 1.0)
#         evaluate!(u0, H, x, 0.0)
#         normalizer .+= abs.(u1)
#         v .+= abs.(u1 .- u0)
#     end
#     normalizer ./= nsamples
#     v .= v ./ nsamples
#     v .= v ./ normalizer

#     return v
# end

# estimate_system_distance(T.tracker.homotopyq)


# x = collect(xs)[1]

# path_info(T.tracker, x)















F = four_bar();

q = read_parameters("four_bar_params_start.txt");
p = read_parameters("four_bar_params_target.txt");
xs = read_solutions("four_bar_sols.txt");

H = ParameterHomotopy(F, q, p; compile = false);

T = Tracker(H)
t = x -> track(T, x)
@time track(T, xs[1])

for i = 1:length(xs)
    @show i
    @time track(T, xs[i])
end


path_info(T, xs[2])
R = map(t, xs)
count(is_success, R)
@time R = map(t, xs);

bad_ = findall(!is_success, R)




bad_paths = [
    82,
    285,
    527,
    702,
    1695,
    1837,
    1845,
    1850,
    1935,
    2008,
    2132,
    2491,
    2501,
    2798,
    2799,
    2808,
    3156,
    3221,
    3992,
    4054,
    4266,
    4318,
    4541,
    5052,
    5213,
    5816,
    5846,
    5888,
    6062,
    6528,
    6672,
    6842,
    6946,
    7132,
    7351,
    7681,
    8221,
    8256,
    8455,
    8478,
]

# p = randn(ComplexF64, 16)

# s = Complex{Float64}[
#     -2.6370887+0.35329306im,
#     -3.7223306+2.5267087im,
#     -0.70272787-0.74162847im,
#     -0.52112188-1.5249im,
#     -0.29444549+0.2698171im,
#     0.23877446+1.1593521im,
#     -0.2226372+0.24523739im,
#     -0.23666974+0.87507397im,
#     2.5277031-6.3353057im,
#     -0.88010494+1.3133313im,
#     -0.73606124+0.15544107im,
#     0.80744564-3.0613903im,
#     -0.82122921+0.80545901im,
#     0.68534751-2.5000172im,
#     0.65872999-2.3719828im,
#     -0.87767826-0.35435132im,
#     -0.93290889+0.12048708im,
#     -0.93106365-0.75512927im,
#     1.8130785-1.6567022im,
#     -0.85699423+0.24221833im,
#     -0.73738109-1.1832401im,
#     -0.81460307+0.2750148im,
#     -0.80200623+0.28313097im,
#     -0.12955283+2.5215805im,
# ]

# s_q = Complex{Float64}[
#     -0.11052551340465448+0.1689556312781767im,
#     0.919958012855403+0.24366762989216786im,
#     0.901696352903974-1.6064926647211917im,
#     -8.338280112808492-8.919160579975916im,
#     0.006804245219856744+0.6666626834147511im,
#     -0.010616472136481038+0.6647920903228585im,
#     0.0054175505582474415+0.6594332738427415im,
#     -0.06615618428899468+0.7160098329660299im,
#     -0.6653465961148477+0.05690565132133756im,
#     -0.8365994540764602+1.2293643774139256im,
#     -1.165781038353456-0.40227703934042947im,
#     -0.636262867891756-7.511530491963335im,
#     -0.1555564468886374-4.36664898084851im,
#     2.1006494673301965-5.4340054725677485im,
#     -0.03216014143666483-1.1684119692446884im,
#     -0.7479475298436062-0.1751684440725321im,
#     1.9041920748392966-0.4938391173146806im,
#     -0.8937602945251837-0.7993076684015589im,
#     -1.8757112160341456+2.124962654367882im,
#     -0.9935684803917256+0.13281722261129017im,
#     -0.9573097412376701+0.22075291383283657im,
#     -0.9207854885450023+0.1388264275883713im,
#     -0.5795474872001224+0.5075857788948817im,
#     1.675306384751507+1.859252783928112im,
# ]


# p = Complex{Float64}[
#     7.297148065259444-3.7969299471707454im,
#     0.06125202773860465+0.9251492821389162im,
#     1.4586512320807385+0.641527879044937im,
#     3.3162292880980537-2.979942348702623im,
#     -0.6613217458488018+0.8321586519581641im,
#     2.6212360706406703-2.6766427578876666im,
#     2.4522388848243555-2.5938408784246976im,
#     -1.7265529790073868-0.29337325907998707im,
#     0.24011427062306998+1.375412535460542im,
#     -0.19336068638852033+0.5669754468190339im,
#     0.147943349839693-0.2784406507316879im,
#     -0.01641513182688631+1.5963735195205522im,
#     -0.33880136827427765+0.2772103200265199im,
#     -0.19788506286822416+1.6632361622942196im,
#     -0.26097545384785903+1.6764055030939693im,
#     0.44012684131549085+1.0212717350758722im,
# ]
# q = Complex{Float64}[
#     -0.35263257737537096+1.2598825875492277im,
#     -0.5353479512983561-0.5963278468670689im,
#     1.4096839876427854-1.2227632666172732im,
#     -0.2550909822371234-0.641270846350799im,
#     0.24180730168363096-0.46611650665985527im,
#     0.29953332004980354-0.8469219548644632im,
#     0.35386303195558433-0.3683866291866925im,
#     1.2021910424897007+1.1148794002014533im,
#     -0.35263257737537096-1.2598825875492277im,
#     -0.5353479512983561+0.5963278468670689im,
#     1.4096839876427854+1.2227632666172732im,
#     -0.2550909822371234+0.641270846350799im,
#     0.24180730168363096+0.46611650665985527im,
#     0.29953332004980354+0.8469219548644632im,
#     0.35386303195558433+0.3683866291866925im,
#     1.2021910424897007-1.1148794002014533im,
# ]


H = ParameterHomotopy(F, q, p; compile = false);

T = Tracker(H)
HomotopyContinuation.path_info(T, s_q)





x₀ = xs[1];
# tracker = AdaptivePathTracker(H; predictor_order = 3);
# tracker3 = AdaptivePathTracker(H; predictor_order = 3);
# tracker4 = AdaptivePathTracker(H; predictor_order = 4);
@time S = adaptive_track(tracker, s_q);
@time tracker_path_info(tracker, s_q)


function threads_track(tracker, xs)
    ys = Vector{TrackerResult}(undef, length(xs))
    trackers = [deepcopy(tracker) for _ = 1:Threads.nthreads()]

    Threads.@threads for i = 1:length(xs)
        ys[i] = track(trackers[Threads.threadid()], xs[i])
    end

    ys
end


ys = []

const t = x -> track(T, x)
@time ys = threads_track(T, xs[1:2000]);
@time R2 = map(t, xs[1:2000]);

D = xs
N = length(D)
count(s -> is_success(s), R2)



R1 = map(s -> adaptive_track(tracker4, s), D);
sum(r -> r.iter, R1) / N
R1_3 = map(s -> adaptive_track(tracker3, s), D);
sum(r -> r.iter, R1_3) / N
sum(r -> steps(r), R2) / N
# count(s -> s.code === :success, R1_3)





M = tracker.state.prec_ComplexF64.M


M.A .= randn.(ComplexF64)


HomotopyContinuation.updated!(M)

x = randn(ComplexF64, 24)
v = randn(ComplexF64, 24)

@benchmark (HomotopyContinuation.updated!($M); HomotopyContinuation.factorize!($M))

@benchmark LA.ldiv!($v, $M, $x)

@benchmark HomotopyContinuation.mixed_precision_iterative_refinement!($v, $M, $x)
@benchmark HomotopyContinuation.fixed_precision_iterative_refinement!($v, $M, $x)


T = Tracker(H)
HomotopyContinuation.path_info(T, s_q)

findall(s -> s.code !== :success, S)

path
@time r = adaptive_track(tracker, xs[1])
@benchmark adaptive_track($tracker, $(xs[1:100]))



@benchmark adaptive_track($H, $x₀, $S, 5000)
t = 1.0



@time rs = map(x -> track(T, x), xs);
count(r -> is_success(r), rs)

@benchmark map(x -> track($T, x), $(xs[1:100]))
@benchmark track($(Tracker(H)), $(xs[1]))


u = zeros(ComplexF64, size(H, 1))
tracker_state = AdaptiveTrackerState{ComplexF64}(H)
predictor = PredictorState(H)

x = copy(x₀)
# predicted
x̂ = similar(x)
# corrected
x̄ = similar(x)

# Check if really zero
evaluate!(u, H, xs[1], 1.0)
if norm(u) > 1e-12
    @warn ("Norm of given value: ", norm(u))
end

t = 1.0
successes = 0
Δt = -1e-2

iter = 0
max_iters = 1000
while abs(Δt) > 16 * eps()
    iter += 1

    evaluate_and_jacobian!(u, tracker_state.prec_ComplexF64.M.A, H, x, t)
    HomotopyContinuation.updated!(tracker_state.prec_ComplexF64.M)

    x̂ = predict(H, x, t, Δt, predictor.cache_ComplexF64, tracker_state.prec_ComplexF64)
    t̂ = t + Δt

    (x̄, code) = correct(H, x̂, t̂)
    if code == :success
        x = x̄
        t = t̂
        successes += 1
        if (successes >= 3)
            Δt *= 2
            successes = 0
        end
        Δt = -min(t, abs(Δt))
    else
        successes = 0
        Δt *= 0.5
    end

    (iter > max_iters) && break

end


return x, t, iter






Arblib.Acb(x::ComplexDF64; prec = 256) = Acb(
    Arb(real(x).hi, prec = prec) + Arb(real(x).lo, prec = prec),
    Arb(imag(x).hi, prec = prec) + Arb(imag(x).lo, prec = prec),
)

F = four_bar()
# S = monodromy_solve(F)
# write_solutions("four_bar_sols.txt", solutions(S))
# write_parameters("four_bar_params_start.txt", parameters(S))
# write_parameters("four_bar_params_target.txt", randn(ComplexF64, length(parameters(S))))

q = read_parameters("four_bar_params_start.txt");
p = read_parameters("four_bar_params_target.txt");
xs = read_solutions("four_bar_sols.txt");

H = ParameterHomotopy(F, q, p; compile = false);
x = ComplexDF64.(xs[1])
t = 1.0
t_start = 1.0
t_target = 0.0
u = zeros(ComplexDF64, size(H, 1))

evaluate!(u, H, xs[1], 1.0)

Base.@kwdef struct PredictorPrecisionState
    cache_ComplexF64::PredictorCache{ComplexF64}
end

Base.@kwdef mutable struct PredictorCache{T}
    tx⁰::TaylorVector{1,T}
    tx¹::TaylorVector{2,T}
    tx²::TaylorVector{3,T}
    tx³::TaylorVector{4,T}
    t::ComplexF64 = complex(NaN)

    u::Vector{T}
    x::Vector{T}
end

function PredictorCache{T}(H::AbstractHomotopy) where {T}
    m, n = size(H)
    tx³ = TaylorVector{4}(T, n)
    PredictorCache(
        tx⁰ = TaylorVector{1}(tx³),
        tx¹ = TaylorVector{2}(tx³),
        tx² = TaylorVector{3}(tx³),
        tx³ = tx³,
        x = zeros(ComplexF64, n),
        u = zeros(ComplexF64, m),
    )
end


function predict(H, x, t, Δt)
    m, n = size(H)
    u = similar(x, n)
    x̂ = similar(x)
    tx³ = TaylorVector{4}(eltype(x), n)
    tx² = TaylorVector{3}(tx³)
    tx¹ = TaylorVector{2}(tx³)
    tx⁰ = TaylorVector{1}(tx³)
    x⁰, x¹, x², x³ = vectors(tx³)

    H_x = similar(x, m, n)
    evaluate_and_jacobian!(u, H_x, H, x, t)

    x⁰ .= x

    taylor!(u, Val(1), H, x, t)
    x¹ .= -(H_x \ u)

    taylor!(u, Val(2), H, tx¹, t)
    x² .= -(H_x \ u)

    taylor!(u, Val(3), H, tx², t)
    x³ .= -(H_x \ u)

    for (i, (xi, xi¹, xi², xi³)) in enumerate(tx³)
        δᵢ = 1 - Δt * xi³ / xi²
        x̂[i] = xi + Δt * (xi¹ + Δt * xi² / δᵢ)
    end


    x̂
end

function correct(H, x̂, t̂)
    m, n = size(H)
    u = similar(x̂, m)
    H_x = similar(x̂, m, n)

    Δx̄ = similar(x̂, m)
    x̄ = copy(x̂)


    for iter = 1:3
        evaluate_and_jacobian!(u, H_x, H, x̄, t̂)
        Δx̄ = H_x \ u
        x̄ = x̄ - Δx̄
        if norm(ComplexF64.(Δx̄)) < 1e-11 * norm(ComplexF64.(x̄))
            return (x̄, :success)
        end
    end
    return (x̄, :failure)
end

function adaptive_track(H::AbstractHomotopy, x₀; max_iters = 5000)
    u = zeros(ComplexF64, size(H, 1))
    predictor = HomotopyContinuation.Predictor(H)
    init!(predictor)
    x = copy(x₀)
    # predicted
    x̂ = similar(x)
    # corrected
    x̄ = similar(x)

    # Check if really zero
    evaluate!(u, H, xs[1], 1.0)
    if norm(u) > 1e-12
        @warn ("Norm of given value: ", norm(u))
    end

    t = 1.0
    successes = 0
    Δt = -1e-2

    iter = 0
    while abs(Δt) > 16 * eps()
        iter += 1

        x̂ = predict(H, x, t, Δt)
        t̂ = t + Δt

        (x̄, code) = correct(H, x̂, t̂)
        if code == :success
            x = x̄
            t = t̂
            successes += 1
            if (successes >= 3)
                Δt *= 2
                successes = 0
            end
            Δt = -min(t, abs(Δt))
        else
            successes = 0
            Δt *= 0.5
        end

        (iter > max_iters) && break
    end


    return x, t, iter
end

PredictorCache{ComplexF64}(H)

@time (y, s, iter) = adaptive_track(H, ComplexDF64.(xs[4]))



x̂ = Arblib.AcbVector(xs[4], prec = 512)
m, n = size(H)
u = similar(x̂, m)
H_x = similar(x̂, m, n)

x̄ = similar(x̂)
x̄ .= x̂

ComplexF64.(x̄)

Arblib.solve!(Arblib.AcbMatrix(x̄), H_x, Arblib.AcbMatrix(u))



@time (y, s, iter) = adaptive_track(H, xs[4])

# STUCK VALUE 
y = ComplexF64[
    -33.213067473353185+196.16004540833595im,
    -4.372213609293898-1.2382436314255036im,
    0.00016494687429948436-0.23993422730066152im,
    -0.013018800881789907-0.23788388709977465im,
    0.06973265575634985-0.043518328548475406im,
    -0.002595700251195619+0.09674979810227978im,
    0.01500781437870132-0.04354957242962884im,
    1.6656964130969656-1.6424036458209432im,
    -0.9782460633689788-0.007355169388765582im,
    -0.17328989701659664+1.2945066183445961im,
    -0.9799614017516131-0.021594169302720706im,
    -0.9681656336296728-0.009581728449883626im,
    -0.9786799692941063-0.010187900339149333im,
    -1.627541258474613-0.18837586094939im,
    -0.9820291843620675-0.0030835896423094334im,
    -0.978383348836015-0.008242659028210843im,
    40.252807070850196+13.947883958413712im,
    -0.6495798739009003-0.548706458043055im,
    22.089837958491632+24.882272894921403im,
    27.8032192067137+8.669392747155422im,
    37.18487730909027+18.24686511263379im,
    -2.4618000603555408+0.4388043673410284im,
    53.054273471980366+9.275104767584017im,
    39.38826902197336+15.400476593822487im,
]
s = 0.6059543885663016

m, n = size(H)
u = similar(y, m)
H_x = similar(y, m, n)
evaluate_and_jacobian!(u, H_x, H, y, s)

m, n = size(H)
Du = ComplexDF64.(similar(y, m))
DH_x = ComplexDF64.(similar(y, m, n))
evaluate_and_jacobian!(Du, DH_x, H, ComplexDF64.(y), s)
evaluate!(Du, H, ComplexDF64.(y), s)



prec = 53
ACBu = AcbMatrix(similar(y, m, 1); prec = prec)p
ACBH_x = AcbMatrix(similar(y, m, n); prec = prec)
ACBy = AcbVector(y; prec = prec)
evaluate_and_jacobian!(ACBu, ACBH_x, H, ACBy, s)
ACBu


Du
DH_x - H_x




@time evaluate!(u, H, y, 0.0)


track(Tracker(H), xs[4])



A0 = (randn(ComplexF64, 8, 8))
b0 = (randn(ComplexF64, 8))

A = A0 .^ 2
b = b0 .^ 2
A1 = ComplexDF64.(A0) .^ 2
b1 = ComplexDF64.(b0) .^ 2
A2 = Arblib.AcbMatrix(Arblib.AcbMatrix(A0) .^ 2)
b2 = Arblib.AcbMatrix(Arblib.AcbMatrix(b0) .^ 2)



c = A \ b
c1 = A1 \ b1

c2 = Arblib.AcbMatrix(8, 1)
Arblib.solve!(c2, A2, b2)


Arblib.Acb(x::ComplexDF64; prec = 256) = Acb(
    Arb(real(x).hi, prec = prec) + Arb(real(x).lo, prec = prec),
    Arb(imag(x).hi, prec = prec) + Arb(imag(x).lo, prec = prec),
)

c2 - c1



c - c0



function gen_taylor_code(M, T)
    taylor_field = Symbol(:taylor_, T)
    order_field = Symbol(:order_, M)
    quote
        I = F.$(taylor_field).$(order_field)
        if isnothing(I)
            I′ = interpreter(TruncatedTaylorSeries{$(M + 1),$T}, F.$(Symbol(:eval_, T)))
            F.$(taylor_field).$(order_field) = I′
            execute_taylor!(
                u,
                Order,
                I′,
                x,
                p;
                assign_highest_order_only = assign_highest_order_only,
            )
        else
            execute_taylor!(
                u,
                Order,
                I,
                x,
                p;
                assign_highest_order_only = assign_highest_order_only,
            )
        end
        u
    end
end

gen_taylor_code()