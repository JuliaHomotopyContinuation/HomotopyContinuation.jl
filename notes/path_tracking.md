# Path tracking algorithm

## Theory

A naive approach to path tracking requires in each step that for the at most $N$-th Newton
update $\Delta x$ it holds $\Vert \Delta x \Vert / \Vert x \Vert < \varepsilon$ for a given (fixed) tolerance $\varepsilon$.
The motivation to limit the number of Newton iterations to $N$ is to avoid _path jumping_.

However, there is one crucial flaw with this approach. If the region of convergence is smaller than $\varepsilon$ at some point along a path, then the path tracking algorithm will get most likely get stuck.

From the Newton–Kantorovich theorem follows that there exists a Lipschitz constant $\omega$ such that for the Newton iterates
$\Vert \Delta x ^{(1)}  \Vert$, $\Vert \Delta x ^{(2)}  \Vert$, ... holds

$$
\Vert \Delta x^{(k+1)} \Vert < \omega \Vert \Delta x ^{(k)}  \Vert^2
$$

(under some mild assumptions).

It roughly holds that $1 / \omega$ is the radius of the region of convergence of Newton's method.
From this follows that $\varepsilon$ needs to be smaller than $1/ \omega$ at each step.

In floating point arithemtic we have the additional difficulty that there exists a maximal achievable accuracy $\mu$ depending on the evaluation error of our system, the conditioning of the Jacobian and accuracy of the machine arithmetic.

We can decrase $\mu$ by computing with mixed precision arithmetic. Let's call this limit accuracy $\mu^*$. In all but the most degenerate cases this will be $\approx \textnormal{eps}$.

At each step we need to have $1 / \omega > \varepsilon > \mu > \mu^*$.
Note that $\mu$ and $\omega$ change at every step.

If we require

$$
\Vert  \Delta x^{(N-1)}  \Vert  \leq
\exp\left(\frac12 \left[\log(1/\omega_{\textnormal{prev}}) + \log(\mu^*)\right]\right)  =
\sqrt{\frac{\mu^*}{\omega_{\textnormal{prev}}}}
=: \tau(\omega_{\textnormal{prev}}, \mu^*),
$$

then follows

$$
\Vert  \Delta x^{(N)}  \Vert  \leq \omega \cdot \omega_{\textnormal {prev}} \cdot \mu^*  \approx \mu^* .
$$

Thus $x^{(N)}$ is almost as accurate as possible.

We need to achieve $\Vert  \Delta x^{(N-1)}  \Vert  \leq  \tau(\omega_{\textnormal{prev}}, \mu^*)$.
Considering

$$
\Vert  \Delta x^{(N-1)}  \Vert \le \omega^{2^{(N-1)}-1} \Vert \Delta x^{(0)} \Vert^{2^{(N-1)}} \le \tau(\omega_{\textnormal{prev}}, \mu^*)
$$

we get the sufficient condition

$$
\Vert \Delta x^{(0)} \Vert \le
\left( \omega^{-2^{(N-1)}+1} \tau(\omega_{\textnormal{prev}}, \mu^*) \right)^{1/2^{(N-1)}} =: \varphi(\omega_{\textnormal{prev}}, \mu^*, N).
$$

For our step size control this implies that the optimal step size $\delta^*$ is

$$
\delta^* := \left[ \;   \,\varphi(\omega_{\textnormal{prev}}, \mu^*, N) / \eta_p  \; \right]^{1/p}
$$

where $\eta_p$ is the prediction error of the predictor with local order $p$.

### Tracking with extended precision

We track with extended precision if the required $\Vert \Delta x^{(N-1)} \Vert$ is smaller than $\beta \mu$ where $\beta$ is a constant buffer.

## Computational testing

The true testing of any path tracking algorithm is in extreme cases. We test the following scenarios

1. Passing near a pole
2. Tracking a diverging path
3. Tracking towards a singular isolated point
4. Tracking towards a higher dimensional component
5. Tracking towards a non-reduced higher dimensional component

### Passing near a pole

```julia
F = let
    @var x a y b x_hat a_hat y_hat b_hat
    gamma = [Variable(Symbol(:gamma, i)) for i = 1:8]
    delta = [Variable(Symbol(:delta, i)) for i = 1:8]
    gamma_hat = [Variable(Symbol(:gamma_hat, i)) for i = 1:8]
    delta_hat = [Variable(Symbol(:delta_hat, i)) for i = 1:8]
    #system of polynomials
    D1 = [
        (a_hat * x - delta_hat[i] * x) * gamma[i] +
        (a * x_hat - delta[i] * x_hat) * gamma_hat[i] +
        (a_hat - x_hat) * delta[i] +
        (a - x) * delta_hat[i] - delta[i] * delta_hat[i] for i = 1:8
    ]
    D2 = [
        (b_hat * y - delta_hat[i] * y) * gamma[i] +
        (b * y_hat - delta[i] * y_hat) * gamma_hat[i] +
        (b_hat - y_hat) * delta[i] +
        (b - y) * delta_hat[i] - delta[i] * delta_hat[i] for i = 1:8
    ]
    D3 = [gamma[i] * gamma_hat[i] + gamma[i] + gamma_hat[i] for i = 1:8]
    System(
        [D1; D2; D3],
        [x; a; y; b; x_hat; a_hat; y_hat; b_hat; gamma; gamma_hat],
        [delta; delta_hat],
    )
end;

s = Complex{Float64}[
    -2.6370887+0.35329306im, -3.7223306+2.5267087im, -0.70272787-0.74162847im, -0.52112188-1.5249im,
    -0.29444549+0.2698171im, 0.23877446+1.1593521im, -0.2226372+0.24523739im, -0.23666974+0.87507397im,
    2.5277031-6.3353057im, -0.88010494+1.3133313im, -0.73606124+0.15544107im, 0.80744564-3.0613903im,
    -0.82122921+0.80545901im, 0.68534751-2.5000172im, 0.65872999-2.3719828im, -0.87767826-0.35435132im,
    -0.93290889+0.12048708im, -0.93106365-0.75512927im, 1.8130785-1.6567022im, -0.85699423+0.24221833im,
    -0.73738109-1.1832401im, -0.81460307+0.2750148im, -0.80200623+0.28313097im, -0.12955283+2.5215805im,
];

p = Complex{Float64}[
    7.297148065259444-3.7969299471707454im, 0.06125202773860465+0.9251492821389162im, 1.4586512320807385+0.641527879044937im, 3.3162292880980537-2.979942348702623im,
    -0.6613217458488018+0.8321586519581641im, 2.6212360706406703-2.6766427578876666im, 2.4522388848243555-2.5938408784246976im, -1.7265529790073868-0.29337325907998707im,
    0.24011427062306998+1.375412535460542im, -0.19336068638852033+0.5669754468190339im, 0.147943349839693-0.2784406507316879im, -0.01641513182688631+1.5963735195205522im,
    -0.33880136827427765+0.2772103200265199im, -0.19788506286822416+1.6632361622942196im, -0.26097545384785903+1.6764055030939693im, 0.44012684131549085+1.0212717350758722im,
];
q = Complex{Float64}[
    -0.35263257737537096+1.2598825875492277im, -0.5353479512983561-0.5963278468670689im, 1.4096839876427854-1.2227632666172732im, -0.2550909822371234-0.641270846350799im,
    0.24180730168363096-0.46611650665985527im, 0.29953332004980354-0.8469219548644632im, 0.35386303195558433-0.3683866291866925im, 1.2021910424897007+1.1148794002014533im,
    -0.35263257737537096-1.2598825875492277im, -0.5353479512983561+0.5963278468670689im, 1.4096839876427854+1.2227632666172732im, -0.2550909822371234+0.641270846350799im,
    0.24180730168363096+0.46611650665985527im, 0.29953332004980354+0.8469219548644632im, 0.35386303195558433+0.3683866291866925im, 1.2021910424897007-1.1148794002014533im,
];

tracker = Tracker(ParameterHomotopy(F, p, q; compile = false));
path_info(tracker, s, 1, 0)
```

Reduced the amount of iterative refinement since not much additional value was provided

### Tracking a singular isolated point

```julia
@var x z
y = 1
# This has two roots of multiplicity 6 at the hyperplane z=0
# each root has winding number 3
F = [
    0.75 * x^4 + 1.5 * x^2 * y^2 - 2.5 * x^2 * z^2 + 0.75 * y^4 - 2.5 * y^2 * z^2 +
    0.75 * z^4
    10 * x^2 * z + 10 * y^2 * z - 6 * z^3
]


Random.seed!(0x123412)
ET, start = total_degree(System(F));
S = collect(start)

tracker = ET.tracker

path_info(tracker, S[1], 1, 1e-8)
```

Decreased the safety factor for the trust region from 0.8 to 0.5

With 0.8

```
┌───┬──────────┬──────────┬───────┬─────────┬─────────┬─────────┬─────────┬─────────┬───────┬─────────┬───────┐
│   │        s │       Δs │     ω │   |Δx₀| │      h₀ │     acc │       μ │       τ │  Δx_t │  Δpred  │   |x| │
├───┼──────────┼──────────┼───────┼─────────┼─────────┼─────────┼─────────┼─────────┼───────┼─────────┼───────┤
│ ✓ │        1 │     -0.2 │     1 │ 0.00044 │ 0.00044 │ 1.2e-16 │ 2.2e-16 │    0.79 │   NaN │ 0.00046 │     1 │
│ ✓ │    0.802 │    -0.29 │     2 │   0.012 │   0.024 │ 8.1e-17 │ 2.2e-16 │    0.37 │   NaN │   0.012 │  0.91 │
│ ✗ │    0.509 │    -0.22 │   2.7 │   0.026 │    0.07 │ 3.7e-15 │ 3.7e-15 │    0.28 │   NaN │    0.33 │  0.91 │
│ ✓ │    0.509 │    -0.11 │   2.7 │  0.0012 │  0.0032 │ 3.7e-15 │ 3.7e-15 │    0.28 │   NaN │  0.0013 │   0.8 │
│ ✓ │      0.4 │    -0.13 │   2.7 │  0.0032 │  0.0087 │ 6.5e-17 │ 2.2e-16 │    0.38 │   NaN │  0.0037 │  0.61 │
│ ✓ │    0.271 │    -0.13 │   2.7 │   0.011 │   0.031 │   1e-16 │ 2.2e-16 │     0.2 │   NaN │   0.013 │  0.57 │
│ ✓ │    0.145 │     -0.1 │   1.7 │  0.0039 │  0.0067 │ 2.5e-16 │ 2.5e-16 │    0.13 │   NaN │  0.0054 │  0.85 │
│ ✗ │   0.0417 │   -0.042 │   4.5 │    0.12 │    0.54 │ 2.5e-16 │ 2.5e-16 │   0.053 │   NaN │     0.3 │  0.85 │
│ ✓ │   0.0417 │   -0.021 │   4.5 │ 8.1e-05 │ 0.00037 │ 2.5e-16 │ 2.5e-16 │   0.053 │   NaN │ 0.00011 │  0.91 │
│ ✗ │   0.0209 │   -0.021 │    26 │    0.11 │     2.9 │ 2.2e-15 │ 2.2e-15 │   0.028 │   NaN │    0.35 │  0.91 │
│ ✓ │   0.0209 │    -0.01 │    26 │ 0.00014 │  0.0036 │ 2.2e-15 │ 2.2e-15 │   0.028 │   NaN │ 0.00017 │  0.95 │
│ ✗ │   0.0104 │    -0.01 │   1.9 │    0.12 │    0.22 │ 7.2e-15 │ 7.2e-15 │   0.015 │   NaN │    0.39 │  0.95 │
│ ✓ │   0.0104 │  -0.0052 │   1.9 │ 0.00019 │ 0.00035 │ 7.2e-15 │ 7.2e-15 │   0.015 │   NaN │ 0.00023 │  0.97 │
│ ✗ │  0.00521 │  -0.0052 │   1.8 │    0.12 │    0.21 │ 6.2e-14 │ 6.2e-14 │  0.0075 │   NaN │    0.42 │  0.97 │
│ ✓ │  0.00521 │  -0.0026 │   1.8 │ 0.00021 │ 0.00038 │ 6.2e-14 │ 6.2e-14 │  0.0075 │   NaN │ 0.00025 │  0.98 │
│ ✗ │  0.00261 │  -0.0026 │     2 │    0.11 │    0.22 │ 1.7e-13 │ 1.7e-13 │  0.0038 │   NaN │    0.45 │  0.98 │
│ ✓ │  0.00261 │  -0.0013 │     2 │ 0.00022 │ 0.00044 │ 1.7e-13 │ 1.7e-13 │  0.0038 │   NaN │ 0.00027 │  0.99 │
│ ✗ │   0.0013 │  -0.0013 │   2.1 │    0.11 │    0.23 │   2e-13 │   2e-13 │  0.0019 │   NaN │    0.46 │  0.99 │
│ ✓ │   0.0013 │ -0.00065 │   2.1 │ 0.00023 │ 0.00048 │   2e-13 │   2e-13 │  0.0019 │   NaN │ 0.00027 │  0.99 │
│ ✗ │ 0.000652 │ -0.00065 │   2.3 │     0.1 │    0.23 │ 5.6e-13 │ 5.6e-13 │ 0.00096 │   NaN │    0.47 │  0.99 │
│ ✓ │ 0.000652 │ -0.00033 │   2.3 │ 0.00023 │ 0.00052 │ 5.6e-13 │ 5.6e-13 │ 0.00096 │   NaN │ 0.00027 │     1 │
│ ✗ │ 0.000326 │ -0.00033 │   2.4 │   0.097 │    0.24 │   6e-13 │   6e-13 │ 0.00048 │   NaN │    0.47 │     1 │
│ ✓ │ 0.000326 │ -0.00016 │   2.4 │ 0.00022 │ 0.00055 │   6e-13 │   6e-13 │ 0.00048 │   NaN │ 0.00027 │     1 │
│ ✗ │ 0.000163 │ -0.00016 │   2.6 │   0.094 │    0.24 │ 3.2e-12 │ 3.2e-12 │ 0.00024 │   NaN │    0.47 │     1 │
│ ✓ │ 0.000163 │ -8.1e-05 │   2.6 │ 0.00022 │ 0.00057 │ 3.2e-12 │ 3.2e-12 │ 0.00024 │   NaN │ 0.00027 │     1 │
│ ✗ │ 8.15e-05 │ -8.1e-05 │   2.7 │   0.091 │    0.24 │   1e-11 │   1e-11 │ 0.00012 │   NaN │    0.47 │     1 │
│ ✓ │ 8.15e-05 │ -4.1e-05 │   2.7 │ 0.00022 │ 0.00059 │   1e-11 │   1e-11 │ 0.00012 │   NaN │ 0.00027 │     1 │
│ ✗ │ 4.07e-05 │ -4.1e-05 │   2.8 │   0.089 │    0.24 │ 3.6e-11 │ 3.6e-11 │ 6.1e-05 │   NaN │    0.47 │     1 │
│ ✓ │ 4.07e-05 │   -2e-05 │   2.8 │ 0.00022 │  0.0006 │ 3.6e-11 │ 3.6e-11 │ 6.1e-05 │   NaN │ 0.00027 │     1 │
│ ✗ │ 2.04e-05 │   -2e-05 │   2.8 │   0.087 │    0.25 │ 4.5e-11 │ 4.5e-11 │   3e-05 │   NaN │    0.47 │     1 │
│ ✓ │ 2.04e-05 │   -1e-05 │   2.8 │ 0.00022 │ 0.00061 │ 4.5e-11 │ 4.5e-11 │   3e-05 │   NaN │ 0.00027 │     1 │
│ ✗ │ 1.02e-05 │   -1e-05 │   2.9 │   0.085 │    0.25 │ 1.3e-10 │ 1.3e-10 │ 1.5e-05 │   NaN │    0.47 │     1 │
│ ✓ │ 1.02e-05 │ -5.1e-06 │   2.9 │ 0.00021 │ 0.00062 │ 1.3e-10 │ 1.3e-10 │ 1.5e-05 │   NaN │ 0.00027 │     1 │
│ ✗ │  5.1e-06 │ -5.1e-06 │   2.9 │   0.083 │    0.24 │ 4.4e-10 │ 4.4e-10 │ 7.6e-06 │   NaN │    0.47 │     1 │
│ ✓ │  5.1e-06 │ -2.5e-06 │   2.9 │ 0.00021 │ 0.00062 │ 4.4e-10 │ 4.4e-10 │ 7.6e-06 │   NaN │ 0.00026 │     1 │
│ ✗ │ 2.56e-06 │ -2.5e-06 │     3 │    0.08 │    0.24 │ 2.7e-11 │ 2.7e-11 │ 3.8e-06 │   NaN │    0.47 │     1 │
│ ✓ │ 2.56e-06 │ -1.3e-06 │     3 │ 0.00021 │ 0.00062 │ 2.7e-11 │ 2.7e-11 │ 3.8e-06 │   NaN │ 0.00026 │     1 │
│ ✗ │ 1.28e-06 │ -1.3e-06 │     3 │   0.076 │    0.23 │ 2.1e-09 │ 2.2e-16 │ 1.9e-06 │   NaN │    0.46 │     1 │
│ ✓ │ 1.28e-06 │ -6.4e-07 │     3 │  0.0002 │  0.0006 │ 2.1e-09 │ 2.2e-16 │ 1.9e-06 │   NaN │ 0.00025 │     1 │
│ ✗ │ 6.47e-07 │ -6.4e-07 │     3 │   0.069 │    0.21 │ 3.3e-16 │ 3.3e-16 │ 9.7e-07 │   NaN │    0.45 │     1 │
│ ✓ │ 6.47e-07 │ -3.2e-07 │     3 │ 0.00019 │ 0.00058 │ 3.3e-16 │ 3.3e-16 │ 9.7e-07 │   NaN │ 0.00024 │     1 │
│ ✗ │ 3.28e-07 │ -3.2e-07 │     3 │   0.057 │    0.17 │ 2.3e-16 │ 2.3e-16 │ 4.9e-07 │   NaN │    0.44 │     1 │
│ ✓ │ 3.28e-07 │ -1.6e-07 │     3 │ 0.00017 │ 0.00053 │ 2.3e-16 │ 2.3e-16 │ 4.9e-07 │   NaN │ 0.00022 │     1 │
│ ✗ │ 1.69e-07 │ -1.6e-07 │     3 │   0.039 │    0.12 │ 2.1e-16 │ 2.2e-16 │ 2.5e-07 │   NaN │    0.41 │     1 │
│ ✓ │ 1.69e-07 │   -8e-08 │     3 │ 0.00014 │ 0.00043 │ 2.1e-16 │ 2.2e-16 │ 2.5e-07 │   NaN │ 0.00018 │     1 │
│ ✗ │ 8.96e-08 │   -8e-08 │     3 │    0.02 │   0.061 │ 1.2e-16 │ 2.2e-16 │ 1.3e-07 │   NaN │    0.37 │     1 │
│ ✓ │ 8.96e-08 │   -4e-08 │     3 │  0.0001 │  0.0003 │ 1.2e-16 │ 2.2e-16 │ 1.3e-07 │   NaN │ 0.00012 │     1 │
│ ✓ │ 4.98e-08 │   -4e-08 │   2.9 │  0.0073 │   0.021 │ 4.4e-17 │ 2.2e-16 │ 7.5e-08 │   NaN │    0.01 │     1 │
└───┴──────────┴──────────┴───────┴─────────┴─────────┴─────────┴─────────┴─────────┴───────┴─────────┴───────┘
```

This ends in the pattern fail -> correct (reducing the previous step size) -> fail (using step size control) -> correct

With 0.5

```
┌───┬──────────┬──────────┬───────┬─────────┬─────────┬─────────┬─────────┬─────────┬───────┬─────────┬───────┐
│   │        s │       Δs │     ω │   |Δx₀| │      h₀ │     acc │       μ │       τ │  Δx_t │  Δpred  │   |x| │
├───┼──────────┼──────────┼───────┼─────────┼─────────┼─────────┼─────────┼─────────┼───────┼─────────┼───────┤
│ ✓ │        1 │     -0.2 │     1 │ 0.00044 │ 0.00044 │       0 │ 2.2e-16 │    0.79 │   NaN │ 0.00046 │     1 │
│ ✓ │    0.802 │    -0.18 │     2 │  0.0011 │  0.0021 │ 5.1e-17 │ 2.2e-16 │    0.37 │   NaN │  0.0012 │  0.98 │
│ ✓ │    0.619 │    -0.14 │   2.7 │  0.0011 │  0.0029 │ 5.4e-17 │ 2.2e-16 │    0.29 │   NaN │  0.0013 │  0.88 │
│ ✓ │    0.476 │    -0.14 │   3.3 │  0.0067 │   0.022 │ 4.5e-17 │ 2.2e-16 │    0.29 │   NaN │  0.0077 │  0.71 │
│ ✓ │    0.332 │    -0.12 │   2.8 │  0.0039 │   0.011 │ 7.6e-17 │ 2.2e-16 │    0.33 │   NaN │  0.0047 │  0.54 │
│ ✓ │    0.211 │   -0.065 │   2.4 │ 0.00046 │  0.0011 │ 4.4e-17 │ 2.2e-16 │    0.13 │   NaN │ 0.00051 │  0.57 │
│ ✓ │    0.146 │   -0.066 │   1.5 │ 0.00032 │ 0.00049 │ 1.3e-16 │ 2.2e-16 │    0.13 │   NaN │ 0.00029 │  0.73 │
│ ✓ │   0.0805 │    -0.05 │   5.3 │  0.0005 │  0.0027 │   2e-16 │ 2.2e-16 │   0.099 │   NaN │ 0.00062 │  0.88 │
│ ✓ │   0.0308 │    -0.02 │    11 │ 0.00087 │  0.0092 │ 2.6e-15 │ 2.6e-15 │    0.04 │   NaN │  0.0011 │  0.95 │
│ ✓ │   0.0107 │  -0.0075 │   4.5 │  0.0022 │  0.0098 │ 3.4e-15 │ 3.4e-15 │   0.015 │   NaN │  0.0028 │  0.98 │
│ ✓ │  0.00319 │  -0.0023 │   2.3 │   0.003 │  0.0069 │ 3.3e-14 │ 3.3e-14 │  0.0046 │   NaN │  0.0039 │  0.99 │
│ ✓ │ 0.000882 │ -0.00065 │     3 │  0.0032 │  0.0095 │ 2.9e-15 │ 2.9e-15 │  0.0013 │   NaN │  0.0044 │     1 │
│ ✓ │ 0.000234 │ -0.00013 │   3.6 │ 0.00045 │  0.0016 │ 2.1e-12 │ 2.1e-12 │ 0.00035 │   NaN │ 0.00061 │     1 │
└───┴──────────┴──────────┴───────┴─────────┴─────────┴─────────┴─────────┴─────────┴───────┴─────────┴───────┘
```

### Tracking a divering path

```julia
Random.seed!(0x123412)
@var x y
f = [2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3, 2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5]
ET, start = total_degree(System(f))
S = collect(start)
tracker = ET.tracker
r = path_info(tracker, S[3], 1, 1e-4)
```

This get's penalized by a too tight trust region. Here it would be better to have the safety factor 0.8 again (here it is just too many too small steps)
