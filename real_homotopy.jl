#  Needs Julia 1.1, HomotopyContinuation 1.0 and MixedSubdivisions 0.3
#  You can install this with
#  using Pkg; Pkg.add("HomotopyContinuation"); Pkg.add("MixedSubdivisions")


# using MixedSubdivisions, HomotopyContinuation, LinearAlgebra
# const HC = HomotopyContinuation
# using Test

# Hacky example to simulate real homotopy
# Using the example from p.382 from Bernds
# Viro Patchworking for Complete Intersections paper

# define support
f_support = [
    0 1 2 3 0 1 2 0 1 0
    3 2 1 0 2 1 0 1 0 0
]
f_lifting = [0, 1, 5, 12, 1, 4, 9, 5, 9, 12]
f_signs = [1, -1, -1, 1, -1, 1, -1, -1, -1, 1]

g_support = [
    0 1 2 0 1 0
    2 1 0 1 0 0
]
g_lifting = [8, 6, 6, 3, 2, 0]
g_signs = [1, -1, 1, -1, -1, 1]

@polyvar x y

# Choose t small enough such that the system has 6 real roots
t = 0.45

coeffs_f = f_signs .* t.^f_lifting
coeffs_g = g_signs .* t.^g_lifting

f = [x^f_support[1,i]*y^f_support[2,i] for i in 1:size(f_support,2)] ⋅ coeffs_f
g = [x^g_support[1,i]*y^g_support[2,i] for i in 1:size(g_support,2)] ⋅ coeffs_g

supp = [f_support,g_support]
lift = [f_lifting,g_lifting]

cells = mixed_cells(supp, lift)

# Check that we will degenerate to binmial cells
@test all(c -> is_fully_mixed_cell(c, supp, lift), cells) == true

# Reuse polyhedra homotopy implementation to simulate the
# real polyehdral homotopy. Arithmetic is over the complex numbers
# but we stay on the real line, so this is the same as doing things over the reals.
supports = [Int32.(f_support), Int32.(g_support)]
complex_coeffs = [complex.(coeffs_f), complex.(coeffs_g)]
start_solutions = PolyhedralStartSolutionsIteratorIterator(supports, complex_coeffs)
start_solutions.mixed_cells = cells
start_solutions.lifting = [f_lifting,g_lifting]

# setup toric homotopy
toric_homotopy = ToricHomotopy(supports, complex_coeffs)
HC.update_lifting!(toric_homotopy, start_solutions.lifting)

# setup path tracker
tracker = coretracker(toric_homotopy, rand(ComplexF64, 2); predictor=Pade21())

S = Vector{Float64}[]
# For each of the cell track from toric infinity
for (cell, X) in start_solutions
    update_cell!(toric_homotopy, cell)
    # (assume each cell has volume 1)
    x∞ = @view X[:,1]
    # We degenerate using exponential coordinates
    # e^s
    # The value e^(-50) is sufficiently small for our purposes here
    res = track(tracker, x∞, -50, 0.0)
    # check that tracking succeeded and store result
    if res.returncode == CoreTrackerStatus.success
        push!(S, real.(res.x))
    end
end

println("Real solutions: ")
for s in S
    println(s)
end
