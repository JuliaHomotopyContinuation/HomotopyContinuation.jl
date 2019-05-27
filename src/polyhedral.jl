export PolyhedralTracker, update_cell!
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
    lifting_range = Int32(-2^15):Int32(2^15)
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


struct PolyhedralTracker{CT<:CoreTracker, PT<:PathTracker, H<:ToricHomotopy}
    toric_homotopy::H
    toric_tracker::CT
    generic_tracker::PT
    # state
    s₀::Base.RefValue{Float64}
end

function PolyhedralTracker(support, coeffs, cell_iter::PolyhedralStartSolutionsIterator; kwargs...)
    toric_homotopy = ToricHomotopy(cell_iter.support, cell_iter.lifting, cell_iter.start_coefficients)
    generic_homotopy = CoefficientHomotopy(support, cell_iter.start_coefficients, coeffs)

    x = randn(ComplexF64, size(support[1], 1))
    toric_tracker =
        coretracker(toric_homotopy, [x], affine_tracking=true, predictor=Pade21())
    generic_tracker =
        pathtracker(generic_homotopy, x, affine_tracking=true, predictor=Heun(), kwargs...)

    PolyhedralTracker(toric_homotopy, toric_tracker, generic_tracker, Ref(NaN))
end

seed(PT::PolyhedralTracker) = PT.generic_tracker.problem.seed

function track!(PT::PolyhedralTracker, x∞)
    retcode = track!(PT.toric_tracker, x∞, PT.s₀[], 0.0)
    if retcode == CoreTrackerStatus.terminated_invalid_startvalue
        PT.s₀[] *= 2
        retcode = track!(PT.toric_tracker, x∞, PT.s₀[], 0.0)
    end
    if retcode != CoreTrackerStatus.success

        return PathTrackerStatus.status(retcode)
    end
    track!(PT.generic_tracker, currx(PT.toric_tracker))
end

@inline function PathResult(PT::PolyhedralTracker, args...; kwargs...)
    PathResult(PT.generic_tracker, args...; kwargs...)
end

function track(tracker::PolyhedralTracker, x₁; path_number::Int=1, details::Symbol=:default, kwargs...)
    track!(tracker, x₁; kwargs...)
    PathResult(tracker, x₁, path_number; details=details)
end

"""
    update_cell!(PT::PolyhedralTracker, cell::MixedCell)

Update the polyhedral tracker `PT` to track paths coming from the mixed cell `cell`.
"""
function update_cell!(PT::PolyhedralTracker, cell::MixedCell)
    min_weight = update_cell!(PT.toric_homotopy, cell)
    # We construct s₀ such that
    #   min_{1 ≤ i ≤ n} min_{aᵢ ∈ Aᵢ} exp(s₀ * w(aᵢ)+aᵢ⋅γ-βᵢ) = 10^-8
    # From this follows:
    PT.s₀[] = -8.0*log(10) / min_weight
    PT
end

function track_paths(PT::PolyhedralTracker, start_solutions::PolyhedralStartSolutionsIterator;
                threading=false, show_progress=true,
                path_result_details::Symbol=:default, save_all_paths=false)
    results = Vector{result_type(PT.generic_tracker)}()

    if show_progress
        progress = ProgressMeter.ProgressUnknown("Tracking paths... ")
    else
        progress = nothing
    end

    stats = SolveStats()
    ntracked = 1
    try
        nthreads = Threads.nthreads()
        if threading && nthreads > 1
            error("TODO")
            # batch_size = 32 * nthreads
            # batches = BatchIterator(start_solutions, batch_size)
            # batch_tracker = BatchTracker(tracker, batches, path_result_details, save_all_paths)
            # k = 0
            # for batch in batches
            #     prepare_batch!(batch_tracker, batch)
            #     ccall(:jl_threading_run, Ref{Cvoid}, (Any,), batch_tracker)
            #
            #     for R in batch_tracker.results
            #         if R !== nothing
            #             push!(results, R)
            #             issuccess(R) && update!(stats, R)
            #         end
            #     end
            #     k += length(batch)
            #     if batch_tracker.interrupted
            #         return results
            #     end
            #
            #     update_progress!(progress, k, stats)
            # end
        else
            for (cell, X) in start_solutions
                update_cell!(PT, cell)
                for i in 1:size(X, 2)
                    x∞ = @view X[:,i]

                    return_code = track!(PT, x∞)
                    if save_all_paths || return_code == PathTrackerStatus.success
                        R = PathResult(PT, x∞, ntracked; details=path_result_details)
                        push!(results, R)

                        if return_code == PathTrackerStatus.success
                            update!(stats, R)
                        end
                    end
                    ntracked += 1
                end
                if ntracked % 32 == 0
                    if progress !== nothing
                        update_progress!(progress, ntracked, stats)
                    end
                end
            end
            if progress !== nothing
                update_progress!(progress, ntracked - 1, stats; finished=true)
            end
        end
    catch e
        if isa(e, InterruptException)
            return results
        else
            rethrow(e)
        end
    end
    results, ntracked - 1
end

function tracker_startsolutions(prob::PolyhedralProblem, startsolutions::PolyhedralStartSolutionsIterator; kwargs...)
    x = randn(ComplexF64, size(prob.toric_homotopy)[2])
    toric_tracker =
        coretracker(prob.toric_homotopy, [x], affine_tracking=true, predictor=Pade21())
    generic_tracker =
        pathtracker(prob.generic_homotopy, x, affine_tracking=isa(prob.tracking_type, AffineTracking), kwargs...)
    tracker = PolyhedralTracker(prob.toric_homotopy, toric_tracker, generic_tracker, Ref(NaN))
    (tracker=tracker, startsolutions=startsolutions)
end
