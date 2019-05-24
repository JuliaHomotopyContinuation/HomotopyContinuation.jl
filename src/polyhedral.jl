import MixedSubdivisions: MixedCell, MixedCellIterator, RegenerationTraverser
import ElasticArrays: ElasticArray

"""
    PolyhedralStartSolutionsIterator(support, lifting, start_coefficients)

Creates an iterator over all mixed cells and their solutions induced by the given `lifting` of the
`support` with the given `start_coefficients`.
The iterator returns a tuple `(cell, X)` where `cell` is a `MixedCell` and
`X` is a `n × D` matrix containing the in each column a solution of the binomial
system corresponding to the mixed cell.

    PolyhedralStartSolutionsIterator(support)

Construct a polyhedral start solutions iterator by picking a random lifting
and random start coefficients.
"""
mutable struct PolyhedralStartSolutionsIterator
    support::Vector{Matrix{Int32}}
    lifting::Vector{Vector{Int32}}
    start_coefficients::Vector{Vector{ComplexF64}}
    mixed_cell_iter::MixedCellIterator{Int32, Int64, RegenerationTraverser{Int32, Int64}}
    X::ElasticArray{ComplexF64, 2, 1} # solutions
    # Cache
    D::Matrix{Int} # binomial system lhs
    b::Vector{ComplexF64} # binomial system rhs
    H::Matrix{Int} # hermite normal form
    U::Matrix{Int} # trafo matrix for hnf
    γ::Vector{Float64} # angle of the elements in b
    μ::Vector{Float64} # γ^U
    S::ElasticArray{Int, 2, 1}
    αs::Vector{Float64} # partial angles in the triangular solve
    DT::Matrix{Float64} # transposed of D to solve coords
end

function PolyhedralStartSolutionsIterator(f::MPPolys)
    PolyhedralStartSolutionsIterator(MixedSubdivisions.support(f))
end
function PolyhedralStartSolutionsIterator(support::Vector{Matrix{Int32}})
    n = size(first(support), 1)
    # Our random lifting strategy is to pick random integers
    # in the range of 0:nterms where nterms is the total number of terms in the support
    lifting_range = Int32(0):Int32(2_000) # ??????????????? #Int32(10 * sum(A -> size(A, 2), support))
    lifting = map(A -> rand(lifting_range, size(A,2)), support)
    # As coefficients do we sample random gaussian numbers with norm one.
    # start_coefficients = map(A -> map(_ -> cis(2π * rand()), 1:size(A,2)), support)
    start_coefficients = map(A -> randn(ComplexF64, size(A,2)), support)
    PolyhedralStartSolutionsIterator(support, lifting, start_coefficients)
end
function PolyhedralStartSolutionsIterator(
            support::Vector{Matrix{Int32}},
            lifting::Vector{Vector{Int32}},
            start_coefficients::Vector{Vector{ComplexF64}})

    n = size(first(support), 1)
    mixed_cell_iter = MixedCellIterator(support, lifting)

    D = zeros(Int, n, n)
    b = zeros(ComplexF64, n)
    H = zeros(Int, n, n)
    U = zeros(Int, n, n)
    γ = zeros(Float64, n)
    μ = zeros(Float64, n)
    X = ElasticArray{ComplexF64}(undef, n, 0)
    S = ElasticArray{Int}(undef, n, 0)
    αs = zeros(Float64, n)
    DT = zeros(Float64, n, n)

    PolyhedralStartSolutionsIterator(support, lifting, start_coefficients,
        mixed_cell_iter, X, D, b, H, U, γ, μ, S, αs, DT)
end

Base.IteratorSize(::Type{<:PolyhedralStartSolutionsIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:PolyhedralStartSolutionsIterator}) = Base.HasEltype()
Base.eltype(iter::PolyhedralStartSolutionsIterator) = Tuple{MixedCell, typeof(iter.X)}

function Base.iterate(iter::PolyhedralStartSolutionsIterator)
    cell_cellstate = iterate(iter.mixed_cell_iter)
    cell_cellstate !== nothing || return nothing

    cell, cellstate = cell_cellstate
    cell_solutions!(iter, cell)

    return (cell, iter.X), cellstate
end

function Base.iterate(iter::PolyhedralStartSolutionsIterator, cellstate)
    cell_cellstate = iterate(iter.mixed_cell_iter, cellstate)
    cell_cellstate !== nothing || return nothing

    cell, cellstate = cell_cellstate
    cell_solutions!(iter, cell)

    return (cell, iter.X), cellstate
end

function cell_solutions!(iter::PolyhedralStartSolutionsIterator, cell::MixedCell)
    @unpack support, start_coefficients, D, b, H, U, γ, μ, X, S, αs, DT = iter
    n = length(iter.support)

    # setup system and compute hermite normal form
    construct_binomial_system!(D, b, support, start_coefficients, cell.indices)
    hnf!(H, U, D)

    # resize memory
    d̂ = cell.volume
    resize!(S, n, d̂)
    resize!(X, n, d̂)

    # setup a table of all combinations of the types of roots of unity we need
    d, e = d̂, 1
    @inbounds for I in 1:n
        dᵢ = H[I,I]
        d = d ÷ dᵢ
        k = 1
        for _ in 1:e, j in 0:(dᵢ-1), _ in 1:d
            S[I, k] = j
            k += 1
        end
        e *= dᵢ
    end

    # We split the computation of the solution in 2 stages
    # 1) Compute solutions for the angular part -> d̂ many solutions
    # 2) Solve for the absolute value

    # 1)
    γ .= angle.(b) ./ (2π)
    apply_trafo_angles!(μ, γ, U)
    @inbounds for i in 1:d̂, j in n:-1:1
        α = μ[j] + S[j, i]
        for k in n:-1:(j+1)
            α -= αs[k] * H[k, j]
        end
        α /= H[j,j]
        X[j, i] = cis(2π * α)
        αs[j] = α
    end
    # 2)
    μ .= log.(abs.(b))
    LinearAlgebra.transpose!(DT, D)
    LinearAlgebra.ldiv!(LinearAlgebra.lu!(DT), μ)
    μ .= exp.(μ)
    for j in 1:d̂, i in 1:n
        X[i, j] *= μ[i]
    end

    nothing
end


function apply_trafo_angles!(μ, γ, U)
    n = length(γ)
    μ .= zero(eltype(μ))
    @inbounds for j in 1:n, i in 1:n
        μ[j] += U[i,j] * γ[i]
    end
    μ
end

function construct_binomial_system!(D, b, support, coeffs, indices)
    n = length(indices)
    @inbounds for i in 1:n
        aᵢ, bᵢ = indices[i]
        for j in 1:n
            D[j, i] = support[i][j, aᵢ] - support[i][j, bᵢ]
        end
        b[i] = -coeffs[i][bᵢ] / coeffs[i][aᵢ]
    end
    nothing

end

# This is just a prototype
export polyhedral_solve

function polyhedral_solve(f)
    support, coeffs = HC.support_coefficients(f)

    cell_iter = HC.PolyhedralStartSolutionsIterator(support)

    polyhedral_homotopy = HC.PolyhedralHomotopy(cell_iter.support, cell_iter.lifting, cell_iter.start_coefficients)
    generic_homotopy = HC.CoefficientHomotopy(support, cell_iter.start_coefficients, coeffs)

    polyhedral_tracker = let
        x = randn(ComplexF64, size(support[1], 1))
        first(coretracker_startsolutions(polyhedral_homotopy, [x], affine_tracking=true, predictor=Pade21()))
    end

    generic_target_tracker = let
        x = randn(ComplexF64, size(support[1], 1))
        pathtracker(generic_homotopy, x, affine_tracking=true)
    end

    path_results = PathResult{Vector{ComplexF64}}[]
    for (cell, X) in cell_iter
        HC.update_cell!(polyhedral_homotopy, cell)
        for k in 1:size(X,2)
            x∞ = @view X[:,k]
            retcode = track!(polyhedral_tracker, x∞, -20, 0.0)
            if retcode == CoreTrackerStatus.success
                result = track(generic_target_tracker, HC.currx(polyhedral_tracker))
                push!(path_results, result)
            end
        end
    end
    path_results
end
