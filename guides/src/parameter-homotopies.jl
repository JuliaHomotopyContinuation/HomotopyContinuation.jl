# <h3 class="section-head" id="parameter_homotopies"><a href="#parameter_homotopies">Parameter Homotopies</a></h3>
#
# Consider the situation in which one has to solve a specific instance of a *parametrized family of polynomial systems*
# ```math
# P = \\{f(x_1,\ldots,x_n,a) = (f_1(x_1,\ldots,x_n,a), \ldots, f_n(x_1\ldots,x_n,a)) \mid a \in \mathbb{R}^m\\}.
# ```
# Often, there is a number ``N``, such that a generic member ``f\in P`` has exactly ``N`` solutions ``x\in\mathbb{R}^n`` with ``f(x)=0``. This ``N`` might be very considerably smaller than the number of solutions of an arbitrary polynomial system not in ``P``. To not destroy the solution structure it is desirable to not leave ``P`` during the homotopy.
#
# The basic `solve` of HomotopyContinuation.jl constructs a straight-line homotopy between the start system ``g`` and the target system ``f``; i.e. ``H(x,t)  = tg + (1-t)f``. When ``P`` is not convex, ``H(x,t)`` might leave the family ``P``. For this reason, we implemented *parametrized homotopies* into HomotopyContinuation.jl. The next example explains its usage.
#
# <h3 class="section-head" id="ellipses"><a href="#ellipses">Example: When are two ellipses tangent?</a></h3>
# The following example is inspired by topological data analysis: suppose that you have a point sample from a manifold ``M\subset \mathbb{R}^n``. An approach to estimate topological features of ``M`` from the sample is by [persistent homology](https://en.wikipedia.org/wiki/Persistent_homology). The idea is as follows. Around each point one puts a ball of radius ``r``. Then one computes the [Čech complex](https://en.wikipedia.org/wiki/Čech_complex) of the union of those balls. [It was argued](https://arxiv.org/pdf/1802.09436.pdf) that it could be beneficial to replace balls by *ellipses*. The obstacle in this approach is to compute when two growing ellipses first meet. This problem can be solved by using homotopy continuation.
#
# In dimension 2 the computational problem is as follows. Let the two ellipses be centered at ``p_1,p_2``, respectively, and be given by two symmetric matrices ``Q_1, Q_2``:
# ```math
# E_i( r ) = \\{x\in \mathbb{R}^2 \mid (x-p_i)^T Q_i^TQ_i(x-p_i) = r^2\\},\; i=1,2.
# ```
# We wish to find the smallest radius ``r`` for which ``E_1( r )\cap E_2( r )`` is not empty. Let ``r^\star`` be the solution for this optimization problem. For a generic choice of ``Q_1`` and ``Q_2`` we have that ``\vert E_1(r^\star)\cap E_2(r^\star) \vert =1`` and ``E_1(r^\star)``, ``E_2(r^\star)`` are tangent. In Julia we translate this into a polynomial system:
#

using HomotopyContinuation, LinearAlgebra
## generate the variables
@polyvar Q₁[1:2, 1:2] Q₂[1:2, 1:2] p₁[1:2] p₂[1:2]
@polyvar x[1:2] r
z₁ = x - p₁
z₂ = x - p₂
## initialize the equations for E₁ and E₂
f₁ = (Q₁ * z₁) ⋅ (Q₁ * z₁) - r^2
f₂ = (Q₂ * z₂) ⋅ (Q₂ * z₂) - r^2
## initialize the equation for E₁ and E₂ being tangent
@polyvar λ
g = (Q₁' * Q₁) * z₁ - λ .* (Q₂' * Q₂) * z₂
## gather everything in one system
F = [f₁; f₂; g];



# An initial solution is given by two circles, each of radius 1,  centered at ``(1,0)`` and ``(-1,0)``, respectively.

map(F) do f
    f(vec(Q₁) => [1,0,0,1], vec(Q₂) => [1,0,0,1], x => [0,0], p₁ => [1,0], p₂ => [-1,0], λ => -1, r => 1)
end


# Let us track this solution to the system given by ``p_1 = [7,5], p_2 = [1,2], Q_1 = \begin{pmatrix} 1 & 2 \\\ 2 & 4 \end{pmatrix}, Q_2 = \begin{pmatrix} 0 & 3 \\\ 3 & 1 \end{pmatrix}``.
#
# That is, the *parameters* are ``p_1, p_2, Q_1, Q_2`` and the *variables* are ``x,r,λ``.
# Now we track the starting solution towards the target system

params = [vec(Q₁); vec(Q₂); p₁; p₂]
startparams = [1, 0, 0, 1, 1, 0, 0, 1, 1, 0, -1, 0]
targetparams = [vec([1 2; 2 5]); vec([0 3; 3 1]); [7, 5]; [1, 2]]
#md solve(F, [[0, 0, 1, -1]], parameters=params, startparameters=startparams, targetparameters=targetparams) #hide
S = solve(F, [[0, 0, 1, -1]], parameters=params, startparameters=startparams, targetparameters=targetparams)


# The computation reveals that ``r^\star \approx 10.89``. We can plot the two ellipses:
# ```julia
# r = solution(S[1])[3]
# r = real(r)
# E₁ = [r .* (inv([1 2; 2 5]) * [cos(2π*t); sin(2π*t)]) + [7; 5] for t in 0:0.01:1]
# E₂ = [r .* (inv([0 3; 3 1]) * [cos(2π*t); sin(2π*t)]) + [1; 2] for t in 0:0.01:1]
#
# # convert E₁, E₂ into matrices
# E₁, E₂ = hcat(E₁...), hcat(E₂...)
#
# # Plot. The Plots package must be installed for this
# using Plots
# plot(E₁[1,:], E₁[2,:], label="Ellipse 1")
# plot!(E₂[1,:], E₂[2,:], label="Ellipse 2")
# ```
# This gives the following picture.
#
# ![img](/images/ellipse.png)
