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
The iterator returns a tuple `(cell, x)` where `cell` is a `MixedCell` and
`x` is a vector matrix containing  a solution of the binomial system corresponding to
the mixed cell.

    PolyhedralStartSolutionsIterator(support)

Construct a polyhedral start solutions iterator by picking a random lifting
and random start coefficients.
"""
mutable struct PolyhedralStartSolutionsIterator{CT<:CoreTracker}
    support::Vector{Matrix{Int32}}
    start_coefficients::Vector{Vector{ComplexF64}}
    lifting::Union{Nothing,Vector{Vector{Int32}}}
    mixed_cells::Union{Nothing,Vector{MixedCell}} # all mixed_cells
    X::ElasticArray{ComplexF64,2,1} # solutions
    # Cache
    D::Matrix{Int} # binomial system lhs
    b::Vector{ComplexF64} # binomial system rhs
    H::Matrix{Int} # hermite normal form
    U::Matrix{Int} # trafo matrix for hnf
    γ::Vector{Float64} # angle of the elements in b
    μ::Vector{Float64} # γ^U
    S::ElasticArray{Int,2,1}
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
    start_coefficients = map(A -> rand_approx_unit(size(A, 2)), support)
    PolyhedralStartSolutionsIterator(
        convert(Vector{Matrix{Int32}}, support),
        start_coefficients,
    )
end

"""
    rand_approx_unit(n::Integer)::Vector{ComplexF64}

This samples uniformly from the rectangle ``[0.9,1.1] × [0,2π]`` and transforms the sampled
values with the map ``(r, φ) ↦ r * e^{i φ}``.
"""
rand_approx_unit(n) = map(_ -> (0.9 + 0.2 * rand()) * cis(2π * rand()), 1:n)

function PolyhedralStartSolutionsIterator(
    support::Vector{Matrix{Int32}},
    start_coefficients::Vector{Vector{ComplexF64}},
)
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
    binomial_tracker = coretracker(BinomialHomotopy(D, b, γ), x; affine_tracking = true)

    PolyhedralStartSolutionsIterator(
        support,
        start_coefficients,
        nothing,
        nothing,
        X,
        D,
        b,
        H,
        U,
        γ,
        μ,
        S,
        αs,
        DT,
        γU,
        x,
        αs_rational,
        binomial_tracker,
    )
end

Base.IteratorSize(::Type{<:PolyhedralStartSolutionsIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{<:PolyhedralStartSolutionsIterator}) = Base.HasEltype()
Base.eltype(iter::PolyhedralStartSolutionsIterator) = Tuple{MixedCell,Vector{ComplexF64}}
function Base.length(iter::PolyhedralStartSolutionsIterator)
    compute_mixed_cells!(iter)
    sum(MixedSubdivisions.volume, iter.mixed_cells)
end

function compute_mixed_cells!(iter::PolyhedralStartSolutionsIterator)
    if isnothing(iter.mixed_cells) || isnothing(iter.lifting)
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

    el = (cell, iter.X[:, 1])
    next_state = size(iter.X, 2) == 1 ? (2, 1) : (1, 2)
    el, next_state
end

function Base.iterate(iter::PolyhedralStartSolutionsIterator, (i, j)::Tuple{Int,Int})
    i > length(iter.mixed_cells) && return nothing

    cell = iter.mixed_cells[i]
    j == 1 && cell_solutions!(iter, cell)
    el = (cell, iter.X[:, j])
    next_state = size(iter.X, 2) == j ? (i + 1, 1) : (i, j + 1)
    return el, next_state
end

start_solution_sample(iter::PolyhedralStartSolutionsIterator) =
    randn(ComplexF64, size(iter.X, 1))

###########################
## SOLVE BINOMIAL SYSTEM ##
###########################
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
        isa(e, OverflowError) || rethrow(e)
        try
            H128, U128 = hnf(D, Int128)
            cell_solutions!(iter, cell, H128, U128)
        catch e128
            isa(e128, OverflowError) || rethrow(e128)

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
    n = size(D, 1)
    for k = 1:size(X, 2), i = 1:n
        yᵢ = one(X[1, k])
        for j = 1:n
            if D[j, i] != 0
                yᵢ *= X[j, k]^D[j, i]
            end
        end
        abs(yᵢ - b[i]) < 1e-10 || return false
    end
    true
end

function cell_solutions_direct!(
    iter::PolyhedralStartSolutionsIterator,
    cell::MixedCell,
    H,
    U,
)
    @unpack D, b, μ, X, DT = iter
    n = size(D, 1)
    # We split the computation of the solution in 2 stages
    # 1) Compute solutions for the angular part -> d̂ many solutions
    compute_angular_part!(iter, cell, H, U)
    # 2) Solve for the absolute value
    μ .= log.(abs.(b))
    LinearAlgebra.transpose!(DT, D)
    LinearAlgebra.ldiv!(LinearAlgebra.lu!(DT), μ)
    μ .= exp.(μ)
    for j = 1:cell.volume, i = 1:n
        X[i, j] *= μ[i]
    end

    nothing
end

function fill_unit_roots_combinations!(S, H, d̂)
    n = size(H, 1)
    d, e = d̂, 1
    @inbounds for I = 1:n
        dᵢ = convert(Int64, H[I, I])
        d = d ÷ dᵢ
        k = 1
        for _ = 1:e, j = 0:(dᵢ-1), _ = 1:d
            S[I, k] = j
            k += 1
        end
        e *= dᵢ
    end
end

function compute_angular_part!(
    iter::PolyhedralStartSolutionsIterator,
    cell::MixedCell,
    H::Matrix{I},
    U::Matrix{I},
) where {I<:Union{Int64,Int128}}
    @unpack γ, b, μ, X, S, αs = iter
    n = length(iter.support)
    γ .= angle.(b) ./ (2π)
    apply_trafo_angles!(μ, γ, U)
    @inbounds for i = 1:cell.volume, j = n:-1:1
        α = μ[j] + S[j, i]
        for k = n:-1:(j + 1)
            α -= αs[k] * H[k, j]
        end
        α /= H[j, j]
        X[j, i] = cis(2π * α)
        αs[j] = α
    end
end


function compute_angular_part!(
    iter::PolyhedralStartSolutionsIterator,
    cell::MixedCell,
    H::Matrix{BigInt},
    U::Matrix{BigInt},
)
    @unpack b, γ, X, S = iter
    n = length(iter.support)
    μ = Vector{BigFloat}(undef, n)
    αs = Vector{BigFloat}(undef, n)
    γ .= angle.(b) ./ (2π)
    apply_trafo_angles!(μ, γ, U)
    @inbounds for i = 1:cell.volume, j = n:-1:1
        # This is always small in our applications
        α = μ[j] + S[j, i]
        for k = n:-1:(j + 1)
            α -= αs[k] * H[k, j]
        end
        α /= H[j, j]
        X[j, i] = cis(2π * convert(Float64, mod(α, 1.0)))
        αs[j] = α
    end
end

function apply_trafo_angles!(μ, γ, U)
    n = length(γ)
    μ .= zero(eltype(μ))
    @inbounds for j = 1:n, i = 1:n
        μ[j] += U[i, j] * γ[i]
    end
    μ
end

function construct_binomial_system!(D, b, support, coeffs, indices)
    n = length(indices)
    @inbounds for i = 1:n
        aᵢ, bᵢ = indices[i]
        for j = 1:n
            D[j, i] = support[i][j, aᵢ] - support[i][j, bᵢ]
        end
        b[i] = -coeffs[i][bᵢ] / coeffs[i][aᵢ]
    end
    nothing
end


function cell_solutions_homotopy!(
    iter::PolyhedralStartSolutionsIterator,
    cell::MixedCell,
    H,
    U,
)
    @unpack D, b, γ, μ, X, γU, x, S, αs_rational, binomial_tracker = iter

    γ .= sign.(real.(b))
    apply_trafo!(γU, γ, U)

    n = size(D, 1)
    for i = 1:cell.volume
        αs_denom = 1
        for j = n:-1:1
            α = S[j, i]
            denom = lcm(αs_denom, convert(Int, H[j, j]))
            if γU[j] == -1
                α += 1 // 2
                denom *= 2
            end
            for k = n:-1:(j + 1)
                α -= αs_rational[k] * convert(Int, H[k, j] % denom)
            end
            α //= H[j, j]
            αs_denom = lcm(αs_denom, denominator(α))
            αs_rational[j] = (numerator(α) % denominator(α)) // denominator(α)
            x[j] = cis(2π * αs_rational[j])
        end

        ret = track!(binomial_tracker, x, 1.0, 0.0)
        x̄ = current_x(binomial_tracker)
        for j = 1:n
            X[j, i] = x̄[j]
        end
    end

    nothing
end

function apply_trafo!(γU, γ, U)
    for i = 1:length(γ)
        γUᵢ = 1
        for j = 1:size(U, 2)
            if γ[j] == -1 && isodd(U[j, i])
                γUᵢ = -γUᵢ
            end
        end
        γU[i] = γUᵢ
    end
    γU
end


########################
## POLYHEDRAL TRACKER ##
########################

"""
    PolyhedralTracker{CT<:CoreTracker, PT<:PathTracker, H<:ToricHomotopy}

The polyhedral tracker combines the tracking from toric infinity toward the target system
by a two-stage approach.
"""
struct PolyhedralTracker{
    CT<:CoreTracker,
    PT<:PathTracker,
    H<:ToricHomotopy,
} <: AbstractPathTracker
    toric_homotopy::H
    toric_tracker::CT
    generic_tracker::PT
    # state
    s₀::Base.RefValue{Float64}
    curr_mixed_cell::Base.RefValue{Union{Nothing,MixedCell}}
end

function construct_tracker(
    prob::PolyhedralProblem,
    startsolutions::PolyhedralStartSolutionsIterator;
    kwargs...,
)
    x = randn(ComplexF64, size(prob.toric_homotopy)[2])
    toric_tracker = coretracker(
        prob.toric_homotopy,
        [x];
        seed = prob.seed,
        affine_tracking = true,
        logarithmic_time_scale = true,
        from_infinity = true,
        track_cond = false,
    )
    generic_prob = Problem(
        prob.tracking_type,
        prob.generic_homotopy,
        prob.vargroups,
        prob.seed,
        prob.startsolutions_need_reordering,
        prob.regauging_factors,
    )
    generic_tracker = construct_tracker(generic_prob, startsolutions; kwargs...)

    PolyhedralTracker(
        prob.toric_homotopy,
        toric_tracker,
        generic_tracker,
        Ref(NaN),
        Base.RefValue{Union{Nothing,MixedCell}}(nothing),
    )
end

seed(PT::PolyhedralTracker) = seed(PT.generic_tracker)
function PathResult(
    PT::PolyhedralTracker,
    (cell, x)::Tuple{MixedCell,Vector{ComplexF64}},
    path_number = nothing;
    kwargs...,
)
    PathResult(PT.generic_tracker, nothing, path_number; kwargs...)
end


result_type(PT::PolyhedralTracker) = result_type(PT.generic_tracker)
solution(tracker::PolyhedralTracker) = solution(tracker.generic_tracker)
LA.norm(tracker::PolyhedralTracker) = LA.norm(tracker.generic_tracker)

function prepare!(PT::PolyhedralTracker, S::PolyhedralStartSolutionsIterator)
    if S.lifting === nothing
        compute_mixed_cells!(S)
    end
    update_lifting!(PT.toric_homotopy, S.lifting)
    PT
end
function track!(
    PT::PolyhedralTracker,
    (cell, x∞)::Tuple{MixedCell,Vector{ComplexF64}};
    accuracy::Union{Nothing,Float64} = nothing,
    max_corrector_iters::Union{Nothing,Int} = nothing,
)
    update_cell!(PT, cell)
    retcode = track!(PT.toric_tracker, x∞, PT.s₀[], 0.0)
    if is_invalid_startvalue(retcode)
        track!(PT.toric_tracker, x∞, 2 * PT.s₀[], 0.0)
    end
    track!(
        PT.generic_tracker,
        current_x(PT.toric_tracker);
        accuracy = accuracy,
        max_corrector_iters = max_corrector_iters,
    )
end

function track(
    tracker::PolyhedralTracker,
    cell_x::Tuple{MixedCell,Vector{ComplexF64}};
    details::Symbol = :default,
    path_number::Union{Int,Nothing} = nothing,
    accuracy::Union{Nothing,Float64} = nothing,
    max_corrector_iters::Union{Nothing,Int} = nothing,
)
    track!(tracker, cell_x; accuracy = accuracy, max_corrector_iters = max_corrector_iters)
    PathResult(tracker, cell_x, path_number; details = details)
end


"""
    update_cell!(PT::PolyhedralTracker, cell::MixedCell)

Update the polyhedral tracker `PT` to track paths coming from the mixed cell `cell`.
"""
function update_cell!(PT::PolyhedralTracker, cell::MixedCell)
    curr_cell = PT.curr_mixed_cell[]
    if curr_cell !== nothing
        curr_cell === cell && return PT
    end

    PT.curr_mixed_cell[] = cell
    min_weight = update_cell!(PT.toric_homotopy, cell)
    # We construct s₀ such that
    #   min_{1 ≤ i ≤ n} min_{aᵢ ∈ Aᵢ} exp(s₀ * w(aᵢ)+aᵢ⋅γ-βᵢ) = 10^-8
    # From this follows:
    PT.s₀[] = -8.0 * log(10) / min_weight
    PT
end
