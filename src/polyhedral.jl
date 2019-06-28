export PolyhedralStartSolutionsIterator, update_cell!, mixed_volume

import MixedSubdivisions: MixedCell, MixedCellIterator, RegenerationTraverser, mixed_volume
import ElasticArrays: ElasticArray


#####################
## START SOLUTIONS ##
#####################

"""
    PolyhedralStartSolutionsIterator(f)

Creates an iterator over all mixed cells and their solutions corresponding to the polynomial
system `f`.
The iterator returns a tuple `(cell, X)` where `cell` is a `MixedCell` and
`X` is a `n × D` matrix containing the in each column a solution of the binomial
system corresponding to the mixed cell.

    PolyhedralStartSolutionsIterator(support)

Construct a polyhedral start solutions iterator by picking a random lifting
and random start coefficients.
"""
mutable struct PolyhedralStartSolutionsIterator{CT<:CoreTracker}
    support::Vector{Matrix{Int32}}
    start_coefficients::Vector{Vector{ComplexF64}}
    lifting::Union{Nothing,Vector{Vector{Int32}}}
    mixed_cells::Union{Nothing, Vector{MixedCell}} # all mixed_cells
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
    # binomial homotopy fallback solution
    γU::Vector{Int}
    x::Vector{ComplexF64}
    αs_rational::Vector{Rational{Int}}
    binomial_tracker::CT
end

function PolyhedralStartSolutionsIterator(f::MPPolys)
    PolyhedralStartSolutionsIterator(MixedSubdivisions.support(f))
end
function PolyhedralStartSolutionsIterator(support::Vector{<:Matrix{<:Integer}})
    # Intuitively we want to sample the coefficients from the normal distribution
    # because we like this one very much.
    # But this yields to a bad surpise for the start solutions in the case that
    # the ratios of the coefficients in the binomial system is very different in between
    # the equations
    start_coefficients = map(A -> rand_approx_unit(size(A,2)), support)
    PolyhedralStartSolutionsIterator(convert(Vector{Matrix{Int32}}, support), start_coefficients)
end

"""
    rand_approx_unit(n::Integer)::Vector{ComplexF64}

This samples uniformly from the rectangle ``[0.9,1.1] × [0,2π]`` and transforms the sampled
values with the map ``(r, φ) ↦ r * e^{i φ}``.
"""
rand_approx_unit(n) = map(_ -> (0.9 + 0.2 * rand()) * cis(2π * rand()), 1:n)

function PolyhedralStartSolutionsIterator(
            support::Vector{Matrix{Int32}},
            start_coefficients::Vector{Vector{ComplexF64}})

    n = size(first(support), 1)
    D = zeros(Int, n, n)
    b = ones(ComplexF64, n)
    H = zeros(Int, n, n)
    U = zeros(Int, n, n)
    γ = ones(Float64, n)
    μ = zeros(Float64, n)
    X = ElasticArray{ComplexF64}(undef, n, 0)
    S = ElasticArray{Int}(undef, n, 0)
    αs = zeros(Float64, n)
    DT = zeros(Float64, n, n)

    γU = zeros(Int, n)
    x = randn(ComplexF64, n)
    αs_rational = Vector{Rational{Int64}}(undef, n)
    binomial_tracker = coretracker(BinomialHomotopy(D, b, γ), x;
                        affine_tracking=true,
                        initial_step_size=0.2,
                        predictor=Pade21())

    PolyhedralStartSolutionsIterator(support, start_coefficients, nothing, nothing,
            X, D, b, H, U, γ, μ, S, αs, DT, γU, x, αs_rational, binomial_tracker)
end

Base.IteratorSize(::Type{<:PolyhedralStartSolutionsIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:PolyhedralStartSolutionsIterator}) = Base.HasEltype()
Base.eltype(iter::PolyhedralStartSolutionsIterator) = Tuple{MixedCell, typeof(iter.X)}
function Base.length(iter::PolyhedralStartSolutionsIterator)
    compute_mixed_cells!(iter)
    sum(MixedSubdivisions.volume, iter.mixed_cells)
end

function compute_mixed_cells!(iter::PolyhedralStartSolutionsIterator)
    if isnothing(iter.mixed_cells) || isnothing(iter.lifting)
        res = MixedSubdivisions.fine_mixed_cells(iter.support)
        if isnothing(res)
            throw(OverflowError("Cannot compute a start system due to an overflow in the mixed subdivision algorithm."))
        end
        mixed_cells, lifting = res
        iter.mixed_cells = mixed_cells
        iter.lifting = lifting
    end
    iter
end
function Base.iterate(iter::PolyhedralStartSolutionsIterator)
    compute_mixed_cells!(iter)

    isempty(iter.mixed_cells) && return nothing

    cell = iter.mixed_cells[1]
    cell_solutions!(iter, cell)

    return (cell, iter.X), 1
end

function Base.iterate(iter::PolyhedralStartSolutionsIterator, cell_number::Int)
    cell_number ≥ length(iter.mixed_cells) && return nothing

    cell = iter.mixed_cells[cell_number + 1]
    cell_solutions!(iter, cell)

    return (cell, iter.X), cell_number + 1
end

function cell_solutions!(iter::PolyhedralStartSolutionsIterator, cell::MixedCell)
    @unpack support, start_coefficients, D, b, H, U = iter
    n = length(iter.support)

    # setup system and compute hermite normal form
    construct_binomial_system!(D, b, support, start_coefficients, cell.indices)

    # This beauty tries first to compute the hnf using Int64,
    # then falls back to Int128 and finally gives up and just uses
    # BigInt
    try
        hnf!(H, U, D)
        cell_solutions!(iter, cell, H, U)
    catch e
        isa(e, OverflowError) || rethrow(e)
        try
            H128, U128 = hnf(D, Int128)
            cell_solutions!(iter, cell, H128, U128)
        catch e128
            isa(e128, OverflowError) || rethrow(e128)

            H_big, U_big = hnf(D, BigInt)
            cell_solutions!(iter, cell, H_big, U_big)
        end
    end

    nothing
end

function cell_solutions!(iter::PolyhedralStartSolutionsIterator, cell::MixedCell, H, U)
    @unpack D, b, γ, μ, X, S, DT = iter
    n = size(H, 1)
    # resize memory
    d̂ = cell.volume
    resize!(S, n, d̂)
    resize!(X, n, d̂)

    # setup a table of all combinations of the types of roots of unity we need
    fill_unit_roots_combinations!(S, H, d̂)

    cell_solutions_direct!(iter, cell, H, U)
    if !verify_solutions(iter)
        cell_solutions_homotopy!(iter, cell, H, U)
    end
end

function verify_solutions(iter)
    @unpack D, b, X = iter
    n = size(D,1)
    for k in 1:size(X,2)
        for i in 1:n
            yᵢ = one(X[1,k])
            for j in 1:n
                if D[j, i] != 0
                    yᵢ *= X[j,k]^D[j, i]
                end
            end
            abs(yᵢ - b[i]) < 1e-10 || return false
        end
    end
    true
end

function cell_solutions_direct!(iter::PolyhedralStartSolutionsIterator, cell::MixedCell, H, U)
    @unpack D, b, μ, X, DT = iter
    n = size(D,1)
    # We split the computation of the solution in 2 stages
    # 1) Compute solutions for the angular part -> d̂ many solutions
    compute_angular_part!(iter, cell, H, U)
    # 2) Solve for the absolute value
    μ .= log.(abs.(b))
    LinearAlgebra.transpose!(DT, D)
    LinearAlgebra.ldiv!(LinearAlgebra.lu!(DT), μ)
    μ .= exp.(μ)
    for j in 1:cell.volume, i in 1:n
        X[i, j] *= μ[i]
    end

    nothing
end

function fill_unit_roots_combinations!(S, H, d̂)
    n = size(H, 1)
    d, e = d̂, 1
    @inbounds for I in 1:n
        dᵢ = convert(Int64, H[I,I])
        d = d ÷ dᵢ
        k = 1
        for _ in 1:e, j in 0:(dᵢ-1), _ in 1:d
            S[I, k] = j
            k += 1
        end
        e *= dᵢ
    end
end

function compute_angular_part!(iter::PolyhedralStartSolutionsIterator, cell::MixedCell, H::Matrix{I}, U::Matrix{I}) where {I<:Union{Int64,Int128}}
    @unpack γ, b, μ, X, S, αs = iter
    n = length(iter.support)
    γ .= angle.(b) ./ (2π)
    apply_trafo_angles!(μ, γ, U)
    @inbounds for i in 1:cell.volume, j in n:-1:1
        α = μ[j] + S[j, i]
        for k in n:-1:(j+1)
            α -= αs[k] * H[k, j]
        end
        α /= H[j,j]
        X[j, i] = cis(2π * α)
        αs[j] = α
    end
end


function compute_angular_part!(iter::PolyhedralStartSolutionsIterator, cell::MixedCell, H::Matrix{BigInt}, U::Matrix{BigInt})
    @unpack b, γ, X, S = iter
    n = length(iter.support)
    μ = Vector{BigFloat}(undef, n)
    αs = Vector{BigFloat}(undef, n)
    γ .= angle.(b) ./ (2π)
    apply_trafo_angles!(μ, γ, U)
    @inbounds for i in 1:cell.volume, j in n:-1:1
        # This is always small in our applications
        α = μ[j] + S[j, i]
        for k in n:-1:(j+1)
            α -= αs[k] * H[k, j]
        end
        α /= H[j,j]
        X[j, i] = cis(2π * convert(Float64, mod(α, 1.0)))
        αs[j] = α
    end
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


function cell_solutions_homotopy!(iter::PolyhedralStartSolutionsIterator, cell::MixedCell, H, U)
    @unpack D, b, γ, μ, X, γU, x, S, αs_rational, binomial_tracker = iter

    γ .= sign.(real.(b))
    apply_trafo!(γU, γ, U)

    n = size(D,1)
    for i in 1:cell.volume
        αs_denom = 1
        for j in n:-1:1
            α = S[j, i]
            denom = lcm(αs_denom, convert(Int, H[j,j]))
            if γU[j] == -1
                α += 1//2
                denom *= 2
            end
            for k in n:-1:(j+1)
                α -= αs_rational[k] * convert(Int, H[k, j] % denom)
            end
            α //= H[j,j]
            αs_denom = lcm(αs_denom, denominator(α))
            αs_rational[j] = (numerator(α) % denominator(α)) // denominator(α)
            x[j] = cis(2π * αs_rational[j])
        end

        ret = track!(binomial_tracker, x, 1.0, 0.0)
        x̄ = current_x(binomial_tracker)
        for j in 1:n
            X[j,i] = x̄[j]
        end
    end

    nothing
end

function apply_trafo!(γU, γ, U)
    for i in 1:length(γ)
        γUᵢ = 1
        for j in 1:size(U,2)
            if γ[j] == -1 && isodd(U[j, i])
                γUᵢ = -γUᵢ
            end
        end
        γU[i] = γUᵢ
    end
    γU
end


#######################
# POLYHEDRAL TRACKER ##
#######################

"""
    PolyhedralTracker{CT<:CoreTracker, PT<:PathTracker, H<:ToricHomotopy}

The polyhedral tracker combines the tracking from toric infinity toward the target system
by a two-stage approach.
"""
struct PolyhedralTracker{CT<:CoreTracker, PT<:PathTracker, H<:ToricHomotopy}
    toric_homotopy::H
    toric_tracker::CT
    generic_tracker::PT
    # state
    s₀::Base.RefValue{Float64}
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
    track!(PT.generic_tracker, current_x(PT.toric_tracker))
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

    compute_mixed_cells!(start_solutions)
    update_lifting!(PT.toric_homotopy, start_solutions.lifting)

    n = length(start_solutions)
    if show_progress
        progress = ProgressMeter.Progress(n; dt=0.1, desc="Tracking $n paths... ",
                                    delay=0.3, clear_output_ijulia=true)
    else
        progress = nothing
    end

    stats = SolveStats()
    ntracked = 1
    try
        nthreads = Threads.nthreads()
        if threading && nthreads > 1
            @warn "Multithreading is currently not supported with polyhedral homotopy, the computation runs only with a single thread."
        end
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
                if ntracked % 32 == 0
                    if progress !== nothing
                        update_progress!(progress, ntracked, stats)
                    end
                end
            end
        end
        if progress !== nothing
            update_progress!(progress, ntracked - 1, stats; finished=true)
        end
    catch e
        if isa(e, InterruptException)
            return results, ntracked
        else
            rethrow(e)
        end
    end
    results, ntracked - 1
end

function tracker_startsolutions(prob::PolyhedralProblem, startsolutions::PolyhedralStartSolutionsIterator; kwargs...)
    x = randn(ComplexF64, size(prob.toric_homotopy)[2])
    toric_tracker =
        coretracker(prob.toric_homotopy, [x]; seed=prob.seed, affine_tracking=true, predictor=Pade21())
    generic_prob = target_problem(prob)
    generic_tracker, _ =
        tracker_startsolutions(generic_prob, startsolutions; kwargs...)
    tracker = PolyhedralTracker(prob.toric_homotopy, toric_tracker, generic_tracker, Ref(NaN))
    (tracker=tracker, startsolutions=startsolutions)
end

function target_problem(P::PolyhedralProblem)
	Problem(P.tracking_type, P.generic_homotopy, P.vargroups,
				P.seed, P.startsolutions_need_reordering, P.regauging_factors)
end

function start_solution_sample(iter::PolyhedralStartSolutionsIterator)
    randn(ComplexF64, size(iter.X, 1))
end
