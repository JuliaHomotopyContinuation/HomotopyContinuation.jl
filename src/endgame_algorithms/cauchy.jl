export CauchyEndgame

"""
    CauchyEndgame(;kwargs...)

The main idea of the Cauchy Endgame is to use [Cauchy's integral formula](https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula)
to predict the solution of the path ``x(t)``, i.e. ``x(0)``.
At each iteration we are at some point ``(x, t)``. We then track the polygon defined
by ``te^{i2πk/n}`` until we end again at ``x``. Here ``n`` is the number of samples we take
per loop.

The following options are available:
* `samples_per_loop=8`: The number of samples we take at one loop.
* `loopclosed_tolerance=1e-5`: The tolerance when a loop is considered closed.
* `L=0.75` and ``K=0.5``: These are paramters for heuristics. For more details see "A Parallel Endgame " by Bates, Hauenstein and Sommese [^1],
    page 8 and 9.

[^1]: Bates, Daniel J., Jonathan D. Hauenstein, and Andrew J. Sommese. "A Parallel Endgame." Contemp. Math 556 (2011): 25-35.
"""
struct CauchyEndgame <: AbstractEndgameAlgorithm
    samples_per_loop::Int
    loopclosed_tolerance::Float64
    # parameter for the first heuristic
    L::Float64
    # parameter for the second heuristic
    #β::Float64
    K::Float64
end

function CauchyEndgame(;
    samples_per_loop=8,
    loopclosed_tolerance=1e-5,
    L=0.75,
    #β=1e-8,
    K=0.5)
    CauchyEndgame(
        samples_per_loop,
        loopclosed_tolerance,
        L,
        #β,
        K)
end
