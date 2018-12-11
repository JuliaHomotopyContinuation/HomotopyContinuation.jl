# <h3 class="section-head" id="overdetermined_systems"><a href="#overdetermined_systems">Overdetermined systems</a></h3>
#
# A system of polynomial equations $f=(f_1(x_1,\ldots, x_m),\ldots,  f_n(x_1,\ldots,x_m))$ is called *overdetermined*, if it has more equations than variables; i.e., when $n>m$. HomotopyContinuation.jl features [Newtons method for overdetermined systems](https://www.ams.org/journals/mcom/2000-69-231/S0025-5718-99-01115-1/S0025-5718-99-01115-1.pdf) for tracking solutions. We show in an example how it can be used.
#
# <h3 class="section-head" id="Rational_normal_curve"><a href="#Rational_normal_curve">Example: the rational normal curve</a></h3>
#
# The [rational normal curve](https://en.wikipedia.org/wiki/Rational_normal_curve) is a 1-dimensional algebraic variety within the 3-dimensional complex space:
# ```math
# C = \\{(x,y,z) \in \mathbb{C}^3 \mid xz-y^2 = 0,\, y-z^2=0, \, x-yz = 0\\}.
# ```
# Since $C$ is cut out by 3 equations one might expect it to be 0-dimensional, but it is not. In fact, the rational normal curve is a first example of a [non-complete intersection](https://en.wikipedia.org/wiki/Complete_intersection).
#
# <h3 class="section-head" id="tracking_overdetermined"><a href="#tracking_overdetermined">Intersecting the rational normal curve with a linear space</a></h3>
# Let $L$ be a random linear space of dimension 2. Then, almost surely, the intersection $C\cap L$ is finite. The goal of this example is to compute one of the points in $C\cap L$. First, let us generate the equations of $C$ and $L$.

using HomotopyContinuation, LinearAlgebra
@polyvar x y z
f = [x*z-y^2, y-z^2, x-y*z]
ℓ = randn(ComplexF64, 4) ⋅ [x, y, z, 1];

# One is tempted to compute $C\cap L$ by executing `solve([f;ℓ])` and simply compute *all* solutions. But this command will return an error message.
# ```julia
# solve([f;ℓ])
# ```
#
# ```
# ERROR: The input system is overdetermined. Therefore it is necessary to provide an explicit startsystem.
# See
#     https://www.JuliaHomotopyContinuation.org/guides/latest/overdetermined_tracking/
# for details.
# ```
# The problem here is that for overdetermined system we can't construct generic starting systems.
# The reason is simple: a generic overdetermined system has no solutions at all. Nevertheless, if we have a starting system and some of its solutions at hand, we can track them towards $[f, \ell]$.
#
# For instance, the point $p=(1, 1, 1)$ lies both on $C$ and $\\{(x,y,z)\in \mathbb{C}^3 \mid x-y+z -1= 0\\}$. We can track this solution towards $[f, \ell]$ by executing the following:

ℓ₁ = [1, -1, 1, -1] ⋅ [x, y, z, 1];
p = [1, 1, 1, 1];
#md solve([f; ℓ₁], [f; ℓ], [p]) #hide
solve([f; ℓ₁], [f; ℓ], [p])

# The syntax `[p]` is because solve takes as input an *array* of solutions.
