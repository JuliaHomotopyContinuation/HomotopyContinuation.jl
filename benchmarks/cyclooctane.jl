using HomotopyContinuation, LinearAlgebra, DynamicPolynomials

# See https://www.juliahomotopycontinuation.org/examples/cyclooctane/
c² = 2
@polyvar z[1:3, 1:6]
z_vec = vec(z)[1:17] # the 17 variables in a vector
Z = [zeros(3) z[:,1:5] [z[1,6]; z[2,6]; 0] [√c²; 0; 0]] # the eight points in a matrix

# define the functions for cyclooctane:
F1 = [(Z[:, i] - Z[:, i+1]) ⋅ (Z[:, i] - Z[:, i+1]) - c² for i in 1:7]
F2 = [(Z[:, i] - Z[:, i+2]) ⋅ (Z[:, i] - Z[:, i+2]) - 8c²/3 for i in 1:6]
F3 = (Z[:, 7] - Z[:, 1]) ⋅ (Z[:, 7] - Z[:, 1]) - 8c²/3
F4 = (Z[:, 8] - Z[:, 2]) ⋅ (Z[:, 8] - Z[:, 2]) - 8c²/3
f = [F1; F2; F3; F4]

n = 2 # dimension of the cyclooctane variety
N = 17 # ambient dimension
@polyvar Aᵥ[1:n, 1:N] bᵥ[1:n] # variables for the linear equations
p = [vec(Aᵥ); bᵥ] # parameters
F = [f; Aᵥ * z_vec - bᵥ] # the polynomial system we have to solve

# now we solve one particular instance for A,b complex. we use this as start system
A₀ = randn(ComplexF64, n, N)
b₀ = randn(ComplexF64, n)
p₀ = [vec(A₀); b₀]

F₀ = [subs(Fᵢ, p => [vec(A₀); b₀]) for Fᵢ in F]
complex_result = solve(F₀)
println("Degree of the variety: ", nfinite(complex_result)) # should be 1408
