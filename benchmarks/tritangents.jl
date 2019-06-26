using PolynomialTestSystems

# see https://www.juliahomotopycontinuation.org/examples/tritangents/
f = equations(tritangents())
res = solve(f)
println("# nonsingular: ", nnonsingular(res)) # should be 720

res = solve(f; start_system=:polyhedral)
println("# nonsingular: ", nnonsingular(res)) # should be 720
