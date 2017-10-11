export solve

import MultivariatePolynomials
const MP = MultivariatePolynomials

"""
    solve(H::AbstractHomotopy, startvalues_s, [algorithm]; kwargs...)

Solve the homotopy `H` via homotopy continuation with the given `startvalues_s` and the given
`algorithm`.
"""
function solve end

solve(H::AbstractHomotopy, startvalue_s; kwargs...) = solve(H, startvalue_s, SphericalPredictorCorrector(); kwargs...)


"""

    solve(f::MP.AbstractPolynomial{T}; kwargs...)

Solves the polynomial system f.
As optional arguments you can declare the homotopytype (default = GeodesicOnTheSphere) and algorithm (default = SphericalPredictorCorrector()).
"""
function solve(f::MP.AbstractPolynomial{T}; homotopytype=GeodesicOnTheSphere, algorithm=SphericalPredictorCorrector()) where {T<:Number}
      H,s = totaldegree(homotopytype,[f])
      solve(H, s, algorithm)
end

function solve(f::Vector{<:MP.AbstractPolynomial{T}}; homotopytype=GeodesicOnTheSphere, algorithm=SphericalPredictorCorrector()) where {T<:Number}
     H,s = totaldegree(homotopytype,f)
     solve(H, s, algorithm)
end
