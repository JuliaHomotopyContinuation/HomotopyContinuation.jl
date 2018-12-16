export AbstractProblem, Projective, homotopy, homogenization, embed, homvars

abstract type AbstractProblem end

"""
    Projective(H::AbstractHomotopy, homogenization::AbstractHomogenization, seed::Int)

Construct a `ProjectiveProblem`. The homotopy `H` needs to be homogenous.
"""
struct Projective{H<:AbstractHomotopy, N} <: AbstractProblem
    homotopy::H
    vargroups::VariableGroups{N}
    seed::Int
end
function Projective(G::AbstractSystem, F::AbstractSystem,
        vargroups::VariableGroups, seed::Int; homotopy=DEFAULT_HOMOTOPY)
    Projective(homotopy(G, F), vargroups, seed)
end

Base.broadcastable(P::AbstractProblem) = Ref(P)

"""
    homotopy(prob::Projective)

Get the homotopy stored in the problem `prob`.
"""
homotopy(prob::Projective) = prob.homotopy

"""
    homvars(prob::Projective)

Get the homogenization variables of the problem. Returns `nothing` if there are no.
"""
function homvars(prob::Projective{H, N}) where {H,N}
    if prob.vargroups.dedicated_homvars
        map(last, prob.vargroups.groups)
    else
        nothing
    end
end

"""
    embed(prob::Projective, x)

Embed the solution `x` into projective space if necessary.
"""
function embed(prob::Projective{<:AbstractHomotopy, N}, x) where {N}
    dims = projective_dims(prob.vargroups)
    if sum(dims) == length(x)
        ProjectiveVectors.embed(x, dims)
    else
        PVector(x, dims)
    end
end
embed(prob::Projective{<:AbstractHomotopy, N}, x::PVector) where {N} = x
