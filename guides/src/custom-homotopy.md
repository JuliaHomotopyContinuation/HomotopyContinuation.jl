<h3 class="section-head" id="h-introduction"><a href="#h-introduction">Introduction</a></h3>

To track solutions from a start system $G$ to the target system $F$ we use by default the straight-line homotopy
$$H(x,t) := (1-t)F+tG\;.$$
But this is in general *not* the best choice since you usually leave the solution space of your problem.
Therefore we support the ability to define arbitrary homotopies where you have the full power of Julia available.

In the following we will illustrate how to setup a custom homotopy on the following example.
For polynomial systems $F$ and $G$ we want to define the homotopy
$$H(x,t) = (1 - t) F( U(t) x ) +  tG( U(t) x )$$
where $U(t)$ is a random path in the space of [unitary matrices](https://en.wikipedia.org/wiki/Unitary_matrix) with $U(0) = U(1) = I$ and $I$ is the identity matrix.
Such a random path can be constructed by
$$U(t) =
U
\begin{bmatrix}\cos(2\pi t) & -\sin(2\pi t) & 0 &\cdots & 0 \\\\ 
\sin(2\pi t) & \cos(2\pi t) & 0 &\cdots & 0 \\\\ 0 & 0 & 1 &\cdots & 0\\\\0 & 0 & 0 &\ddots & 0\\\\0 & 0 & 0 &\cdots & 1
\end{bmatrix} U^T.
$$
with a random unitary matrix $U$.

<h3 class="section-head" id="h-math"><a href="#h-math">Figuring out the math</a></h3>

To define a homotopy we have to know how to compute for all $x \in \mathbb{C}^n$, $t \in \mathbb{C}$
$$H(x,t), \quad \frac{\partial H}{\partial x}(x,t) \quad \text{ and } \quad \frac{\partial H}{\partial t}(x,t)\;.$$
We denote the partial derivative of $H$ w.r.t. $x$ as the *Jacobian* of $H$.
For simplification (in the math as well as in the implementation) we introduce the helper homotopy
$$\tilde{H}(y, t) := (1 - t) F( y ) +  tG(y)\;.$$
Note $H(x,t) = \tilde{H}(U(t)x, t)$. Using the chain rule we get for the partial derivatives
$$\frac{\partial H}{\partial x}(x,t) = \frac{\partial \tilde{H}}{\partial y}(U(t)x,t) U(t)$$
and
$$\frac{\partial H}{\partial t}(x,t) = \frac{\partial \tilde{H}}{\partial y}(U(t)x,t) U'(t) x + \frac{\partial \tilde{H}}{\partial t}(U(t)x,t) $$
where
$$U'(t)= U
\begin{bmatrix}-2\pi\sin(2\pi t) & -2\pi\cos(2\pi t) & 0 &\cdots & 0 \\\\ 
2\pi\cos(2\pi t) & -2\pi\sin(2\pi t) & 0 &\cdots & 0 \\\\ 0 & 0 & 0 &\cdots & 0\\\\0 & 0 & 0 &\ddots & 0\\\\0 & 0 & 0 &\cdots & 0
\end{bmatrix} U^T.$$

<h3 class="section-head" id="h-data-structure"><a href="#h-data-structure">Constructing the homotopy data structures </a></h3>

A custom homotopy has to satisfy a [certain interface](https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/latest/homotopies.html#Interface-for-custom-homotopies-1). We start with the data structure for the homotopy.
A homotopy is represented by a [`struct`](https://docs.julialang.org/en/stable/manual/types/#Composite-Types-1)
which is a subtype of [`Homotopies.AbstractHomotopy`](https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/latest/homotopies.html#HomotopyContinuation.HomotopiesBase.AbstractHomotopy).
Since $\tilde{H}$ is the standard straight-line homotopy we can reuse this implementation to safe us some work since
homotopies compose easily.
```julia
using HomotopyContinuation

struct RandomUnitaryPath{Start,Target} <: Homotopies.AbstractHomotopy
    straightline::StraightLineHomotopy{Start, Target}
    U::Matrix{Complex128}
end
```

To get good performance it is important to be careful about memory allocations. It is much much better
to initialize a chunk of memory *once* and to reuse this memory. To support this optimization we have the concept
of a *cache*. This is a `struct` with supertype [`Homotopies.AbstractHomotopyCache`](https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/latest/homotopies.html#HomotopyContinuation.HomotopiesBase.AbstractHomotopyCache) where we allocate all memory necessary to evaluate and differentiate our homotopy.
This is an optimization and not necessary to have at the beginning, but for the best it is necessary to implement it.
To illustrate how to do this, we will implement here a cache. Don't look with too much detail on the exact type definition for now, we just allocate a bunch of stuff which will make much more sense later.
As a constructor for the cache we have to define the [`Homotopies.cache`](https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/latest/homotopies.html#HomotopyContinuation.HomotopiesBase.cache) method.
```julia
struct RandomUnitaryPathCache{C, T1, T2} <: Homotopies.AbstractHomotopyCache
    straightline::C
    U_t::Matrix{Complex128}
    y::Vector{T1}
    # More temporary storage necessary to avoid allocations
    jac::Matrix{T2} # holds a jacobian
    dt::Vector{T2} # holds a derivative w.r.t. t
    U::Matrix{Complex128} # holds something like U
end

# A cache is always constructed by this method.
function Homotopies.cache(H::RandomUnitaryPath, x, t)
    U_t = copy(H.U)
    y = U_t * x
    straightline = Homotopies.cache(H.straightline, y, t)

    jac = Homotopies.jacobian(H.straightline, y, t, straightline)
    dt = jac * y
    U = copy(U_t)
    RandomUnitaryPathCache(straightline, U_t, y, jac, dt, U)
end
```

<h3 class="section-head" id="h-implementation"><a href="#h-implementation">Implementing the evaluation and derivatives</a></h3>

We start with implementing subroutines to evaluate and differentiate $U(t)$ as well as computing $U(t)x$ and $U'(t)x$.
We use the `U_t` and `y` fields in `cache` to store the values $U(t)$ and $U(t)x$ (resp. $U'(t)$ and $U'(t)x$)
```julia
# U(t)x
function Ut_mul_x!(cache, U, x, t)
    # We start with U * (the 2x2 sin-cos block + I)
    cache.U .= U
    s, c = sin(2π*t), cos(2π*t)
    for i=1:size(U, 1)
        cache.U[i, 1] = U[i,2] * s + U[i,1] * c
        cache.U[i, 2] = U[i,2] * c - U[i,1] * s
    end
    # U(t) = cache.U * U'
    # y = cache.y = U(t) * x
    A_mul_B!(cache.y, A_mul_Bc!(cache.U_t, cache.U, U), x)
end

# U'(t)x
function U_dot_t_mul_x!(cache, U, x, t)
    # We start with U * (the derivative of the 2x2 sin-cos block + 0)
    cache.U .= zero(eltype(U))
    s, c = 2π*sin(2π*t), 2π*cos(2π*t)
    for i=1:size(U, 1)
        cache.U[i, 1] =  U[i,2] * c - U[i,1] * s
        cache.U[i, 2] = -U[i,2] * s - U[i,1] * c
    end

    # U'(t) = cache.U * U'
    # y' = cache.y = U'(t) * x
    A_mul_B!(cache.y, A_mul_Bc!(cache.U_t, cache.U, U), x)
end

```

Now we are ready to implement $H(x,t)$, its Jacobian and the derivative w.r.t. $t$.
```julia
function Homotopies.evaluate!(out, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    Homotopies.evaluate!(out, H.straightline, y, t, cache.straightline)
end

function Homotopies.jacobian!(out, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    Homotopies.jacobian!(cache.jac, H.straightline, y, t, cache.straightline)
    A_mul_B!(out, cache.jac, cache.U_t) # out = J_H(y, t) * U(t)
end

function Homotopies.dt!(out, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    # chain rule
    Homotopies.jacobian_and_dt!(cache.jac, out, H.straightline, y, t, cache.straightline)
    y_dot = U_dot_t_mul_x!(cache, H.U, x, t) # y_dot = U'(t)x
    A_mul_B!(cache.dt, cache.jac, y_dot) # dt = J_H(y, t) * y_dot
    out .+= cache.dt
end
```

We also support to compute `evaluate!` and `jacobian!`  simultaneously as well as to compute
`jacobian!` and `dt!` simultaneously. This can be very beneficial for the performance, so let's implement this here
since this mostly involve copy-paste.
```julia
function Homotopies.evaluate_and_jacobian!(val, jac, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    Homotopies.evaluate_and_jacobian!(val, cache.jac, H.straightline, y, t, cache.straightline)
    A_mul_B!(jac, cache.jac, cache.U_t)
end

function Homotopies.jacobian_and_dt!(jac, dt, H::RandomUnitaryPath, x, t, cache)
    y = Ut_mul_x!(cache, H.U, x, t)
    Homotopies.jacobian_and_dt!(cache.jac, dt, H.straightline, y, t, cache.straightline)
    A_mul_B!(jac, cache.jac, cache.U_t) # jac = J_H(y, t) * U(t)
    y_dot = U_dot_t_mul_x!(cache, H.U, x, t) # y_dot = U'(t)x
    A_mul_B!(cache.dt, cache.jac, y_dot) # dt = J_H(y, t) * y_dot
    dt .+= cache.dt
end
```

<h3 class="section-head" id="h-testing"><a href="#h-testing">Testing our implementation</a></h3>

Implementing these methods without ever testing them is ... not a good idea. Also just throwing the homotopy into `solve`
can result in confusing error messages. Therefore we provide a set of tests against which we can test our implementation.
Although they do not check that we implemented the math correctly (or that our math is correct in the first place), they will
catch any runtime errors.
```julia-repl
# Let us construct some test systems
julia> @polyvar x y z;

julia> F = SPSystem([x^2*y-3x*z, z^2*x+3y^2]);
julia> G = SPSystem([z*x^2-3x*y^2, z^3*x-2x*y*z^2]);

# Here we can test that our implementation does not produce an error
julia> InterfaceTest.homotopy(RandomUnitaryPath(G, F))
Test Passed
```

<h3 class="section-head" id="h-using"><a href="#h-using">Using our new homotopy</a></h3>
We are now ready to use our new homotopy!
```julia-repl
julia> solve([x^2 - y, y^3*x-x], homotopy=RandomUnitaryPath)
-----------------------------------------------
Paths tracked: 8
# non-singular finite solutions:  7
# singular finite solutions:  0
# solutions at infinity:  1
# failed paths:  0
Random seed used: 487296
-----------------------------------------------
```
Alternatively we could also construct the homotopy directly and give it to `solve` together with start solutions.
Note that in this case we have to ensure that our homotopy is already homogenous.

The source code of this example is available [here](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl/blob/master/examples/custom-homotopy.jl).
