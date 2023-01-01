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
\exp\left(\frac12 \left[\log(1/\omega_{\textnormal{prev}}) + \log(\mu^*)\right]\right) =: \tau(\omega_{\textnormal{prev}}, \mu^*),
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
\delta^* := \left[ \; \eta_p  \,\varphi(\omega_{\textnormal{prev}}, \mu^*, N)  \; \right]^{1/p}
$$

where $\eta_p$ is the prediction error of the predictor with local order $p$.

## Computational testing

The true testing of any path tracking algorithm is in extreme cases. We test the following scenarios

1. Passing near a singular point
2. Tracking a diverging path
3. Tracking towards a singular isolated point
4. Tracking towards a higher dimensional component
5. Tracking towards a non-reduced higher dimensional component
