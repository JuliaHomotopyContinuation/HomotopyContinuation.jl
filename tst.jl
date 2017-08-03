using Revise
using HomotopyContinuation

x = MPoly.generator(Complex128, :x)
f = (x - 2.0) * (x - (2.5+ 4.0im))
g = (x - 4.3im) * (x + (2.1 - 4im))

F = MPoly.system([f])
G = MPoly.system([g])

H = StraightLineHomotopy(G, F)

res = solve(H, [[-2.1 + 4.0im], [4.3im]], PredictorCorrector.Spherical())