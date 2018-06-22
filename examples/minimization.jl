using HomotopyContinuation

# We want so solve the problem
#
#   minimize J = 3x^3y+y^2z^2-2x*y-4z^3 s.t. x^2+y^2+z^2= 1

@polyvar x y z
J = 3x^3*y+y^2*z^2-2x*y-x*4z^3
g = x^2+y^2+z^2-1
# Introduce auxillary variable for Lagrangian
@polyvar λ
# define Lagrangian
L = J - λ * g
# compute the gradient
∇L = map(var -> differentiate(L, var), [x, y, z, λ])
# Now we solve the polynomial system
result = solve(∇L)
# We are only interested in the real solutions
reals = realsolutions(result)
# Now we simply evaluate the objective J and find the minimum
minval, minindex = findmin(map(s -> J(s[1:3]), reals))
argmin = reals[minindex][1:3]
# And we see that the minimum is attained at
#     [ 0.496876, 0.0946083, 0.862649]
# with value -1.3284
