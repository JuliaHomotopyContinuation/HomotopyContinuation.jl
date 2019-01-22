# <h3 class="section-head" id="monodromy"><a href="#monodromy">The monodromy method</a></h3>
#
# An alternative to using the [solve function](/guides/examples.md) is solving a polynomial system ``F=(f_1,\ldots,f_n)`` by monodromy. This approach is more efficient, but requires the user to provide at least one solution of ``F=0``. Here is the basic idea:
#
# Suppose ``x`` is a solution ``F(x)=0`` and that ``F=F_{u_0}`` is a point in a family of polynomial systems ``\mathcal{F}=\\{F_u : u\in \mathbb{C}^k\\}`` which is defined with ``k\geq 1`` parameters. The monodromy method consists in moving around ``u`` in a loop starting and ending at ``u_0`` while tracking ``x`` along that loop. After one iteration usually one has found a new solution ``y\in \mathbb{C}^n``. This process then is repeated until some stopping criterion is fulfilled.
#
# The general syntax for this is:
# ```julia
# monodromy_solve(F_u, [x], u₀, parameters = u)
# ```
#
# <h3 class="section-head" id="monodromyexample"><a href="#monodromyexample">Example: computing the closest point on a variety</a></h3>
#
# Consider the problem of computing the point on the variety
# ```math
# V = \\{x=(x_1,x_2)^T\in \mathbb{R}^2 : f(x) = 0\\}, \text{ for } f(x) = x_1^2 + x_2^2 - (x_1^2 + x_2^2 + x_2)^2,
# ```
# which minimizes the distance to the point ``u₀ = (-3,-2)``. The situation looks like this:
#
# <p style="text-align:center;"><img src="/images/cardioid0.png" maxwidth="400px"/></p>
#
# The minimizer ``x^\star`` is a solution to the system
# ```math
# F_u = \begin{bmatrix}\det(\begin{bmatrix}x-u & \nabla_x(f)\end{bmatrix})\\\ f(x)\end{bmatrix} =0,
# ```
#
# where ``\nabla_x(f)`` is the gradient of ``f`` at ``x``. Let's set up this system in Julia.
#

using DynamicPolynomials, LinearAlgebra

u₀ = [-2; -1]

@polyvar x[1:2]
f = x[1]^2 + x[2]^2 - (x[1]^2 + x[2]^2 + x[1])^2
∇ = differentiate(f, x)

F = [det([x-u₀ ∇]); f]
#md nothing #hide

# Now, ``x^\star`` is a zero of `F`, which has totaldegree equal to 12. However, the actual number of solutions is only 3 as was shown in [this article](https://arxiv.org/pdf/1309.0049.pdf). For avoiding computing all 12 paths, we use monodromy instead of totaldegree traicking.
#
#
# An initial solution to `F` is ``x_0=(0,1)^T``. Let's use this initial solution for monodromy:

using HomotopyContinuation

x₀ = [0; 1]
@polyvar u[1:2]
F_u = [det([x-u ∇]); f]

monodromy_solve(F_u, [x₀], u₀, parameters = u, target_solutions_count=4)


# We get the three solutions. The following picture shows them:
#
# <p style="text-align:center;"><img src="/images/cardioid.png" maxwidth="400px"/></p>
#
# The minimizer is ``x^\star \approx (-1.68, -0.86)^T``.
