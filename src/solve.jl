export solve

"""
    solve(H::AbstractHomotopy, startvalues_s, [algorithm]; kwargs...)

Solve the homotopy `H` via homotopy continuation with the given `startvalues_s` and the given
`algorithm`.
"""
function solve end

solve(H::AbstractHomotopy, startvalue_s; kwargs...) = solve(H, startvalue_s, SphericalPredictorCorrector(); kwargs...)
