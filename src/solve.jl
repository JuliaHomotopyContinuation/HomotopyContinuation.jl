export solve

"""
    solve(H::AbstractHomotopy, startvalues_s, [algorithm]; kwargs...)

Solve the homotopy `H` via homotopy continuation with the given `startvalues_s` and the given
`algorithm`.


    solve(f::Vector{<:MP.AbstractPolynomial{T}}; homotopytype=GeodesicOnTheSphere, algorithm=SphericalPredictorCorrector(), kwargs...)
    solve(f::MP.AbstractPolynomial{T}; homotopytype=GeodesicOnTheSphere, algorithm=SphericalPredictorCorrector(), kwargs...)

Solves the polynomial system `f` via homotopy continuation. This uses a totaldegree homotopy
of the given `homotopytype` and uses the given `algorithm` for pathtracking.
"""
function solve end

function solve(H::AbstractHomotopy, startvalue_s; kwargs...)
    solve(H, startvalue_s, SphericalPredictorCorrector(); kwargs...)
end

function solve(
    f::MP.AbstractPolynomial{T};
    homotopytype=GeodesicOnTheSphere,
    algorithm=SphericalPredictorCorrector(),
    kwargs...) where {T<:Number}
      H, s = totaldegree(homotopytype, [f])
      solve(H, s, algorithm; kwargs...)
end
function solve(f::Vector{<:MP.AbstractPolynomial{T}};
    homotopytype=GeodesicOnTheSphere,
    algorithm=SphericalPredictorCorrector(),
    kwargs...) where {T<:Number}
     H, s = totaldegree(homotopytype, f)
     solve(H, s, algorithm; kwargs...)
end
