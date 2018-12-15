export AbstractProblem, Projective, homotopy, homogenization, embed

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
    homogenization(prob::Projective)

Get the [`HomogenizationStrategy`](@ref) stored in the problem `prob`.
"""
homogenization(prob::Projective) = prob.homogenization

"""
    embed(prob::Projective, x)

Embed the solution `x` into projective space if necessary.
"""
function embed(prob::Projective{<:AbstractHomotopy, N}, x) where {N}
    # M, N = size(homotopy(prob))
    dims = projective_dims(prob.vargroups)
    if sum(dims) == length(x)
        ProjectiveVectors.embed(x, dims)
    else
        PVector(x, dims)
    end
end
# function embed(prob::Projective{<:AbstractHomotopy, Homogenization}, x)
#     M, N = size(homotopy(prob))
#     if length(x) == N
#         return PVector(x, homogenization(prob).homvaridx)
#     elseif length(x) == N - 1
#         return ProjectiveVectors.embed(x, homogenization(prob).homvaridx)
#     end
#     throw(error("The length of the initial solution is $(length(x)) but expected length $N or $(N-1)."))
# end
# embed(prob::Projective{<:AbstractHomotopy, Homogenization}, x::PVector) = x
# embed(prob::Projective{<:AbstractHomotopy, NullHomogenization}, x::PVector) = x
