# Start Solutions
import MixedSubdivisions: MixedCell

struct PolyhedralStartSolutions
    support::Vector{Matrix{Int32}}
    start_coefficients::Vector{Vector{ComplexF64}}
    lifting::Vector{Vector{Int32}}
    mixed_cells::Vector{MixedCell}
    BSS::BinomialSystemSolver
end

function PolyhedralStartSolutions(
    support::AbstractVector{<:AbstractMatrix{<:Integer}},
    coeffs::AbstractVector{<:AbstractVector{<:Number}},
)
    BSS = BinomialSystemSolver(length(support))
    PolyhedralStartSolutions(
        convert(Vector{Matrix{Int32}}, support),
        coeffs,
        map(c -> zeros(Int32, length(c)), coeffs),
        MixedCell[],
        BSS,
    )
end

Base.IteratorSize(::Type{<:PolyhedralStartSolutions}) = Base.HasLength()
Base.IteratorEltype(::Type{<:PolyhedralStartSolutions}) = Base.HasEltype()
Base.eltype(iter::PolyhedralStartSolutions) = Tuple{MixedCell,Vector{ComplexF64}}
function Base.length(iter::PolyhedralStartSolutions)
    compute_mixed_cells!(iter)
    sum(MixedSubdivisions.volume, iter.mixed_cells)
end

function compute_mixed_cells!(iter::PolyhedralStartSolutions)
    if isempty(iter.mixed_cells)
        res = MixedSubdivisions.fine_mixed_cells(iter.support)
        if isnothing(res) || isempty(res[1])
            throw(OverflowError("Cannot compute a start system."))
        end
        mixed_cells, lifting = res
        empty!(iter.mixed_cells)
        append!(iter.mixed_cells, mixed_cells)

        for (i, w) in enumerate(lifting)
            empty!(iter.lifting[i])
            append!(iter.lifting[i], w)
        end
    end
    iter
end

function Base.iterate(iter::PolyhedralStartSolutions)
    compute_mixed_cells!(iter)

    isempty(iter.mixed_cells) && return nothing

    cell = iter.mixed_cells[1]
    solve(iter.BSS, iter.support, iter.start_coefficients, cell)

    el = (cell, iter.BSS.X[:, 1])
    next_state = size(iter.BSS.X, 2) == 1 ? (2, 1) : (1, 2)
    el, next_state
end

function Base.iterate(iter::PolyhedralStartSolutions, (i, j)::Tuple{Int,Int})
    i > length(iter.mixed_cells) && return nothing

    cell = iter.mixed_cells[i]
    j == 1 && solve(iter.BSS, iter.support, iter.start_coefficients, cell)
    el = (cell, iter.BSS.X[:, j])
    next_state = size(iter.BSS.X, 2) == j ? (i + 1, 1) : (i, j + 1)
    return el, next_state
end

# Tracker
function polyehdral_system(support)
    n = length(support)
    m = sum(A -> size(A, 2), support)
    @var x[1:n] c[1:m]
    k = 1
    System(map(support) do A
        fi = Expression(0)
        for a in eachcol(A)
            ModelKit.add!(fi, fi, c[k] * prod(x .^ a))
            k += 1
        end
        fi
    end, x, c)
end

"""
    PolyhedralTracker
The polyhedral tracker combines the tracking from toric infinity toward the target system
by a two-stage approach.
"""
struct PolyhedralTracker{H1<:ToricHomotopy,H2,N,M<:AbstractMatrix{ComplexF64}}
    toric_tracker::Tracker{H1,N,Vector{ComplexF64},Vector{ComplexDF64},M}
    generic_tracker::PathTracker{H2,N,Vector{ComplexF64},Vector{ComplexDF64},M}
    support::Vector{Matrix{Int32}}
    lifting::Vector{Vector{Int32}}
end

Base.broadcastable(T::PolyhedralTracker) = Ref(T)

function polyhedral(f::ModelKit.System)
    supp_coeffs = exponents_coefficients.(f, Ref(variables(f)))
    support = first.(supp_coeffs)
    target_coeffs = map(ci -> float.(ModelKit.to_number.(ci)), last.(supp_coeffs))
    start_coeffs =
        map(c -> exp.(randn.(ComplexF64) .* 0.1 .+ log.(complex.(c))), target_coeffs)
    F = ModelKit.compile(polyehdral_system(support))

    H₁ = ToricHomotopy(F, start_coeffs)
    toric_tracker = Tracker(H₁)

    H₂ = begin
        p = reduce(append!, start_coeffs; init = ComplexF64[])
        q = reduce(append!, target_coeffs; init = ComplexF64[])
        ParameterHomotopy(ModelKitSystem(F), p, q)
    end
    generic_tracker = PathTracker(Tracker(H₂))

    S = PolyhedralStartSolutions(support, start_coeffs)
    tracker = PolyhedralTracker(toric_tracker, generic_tracker, S.support, S.lifting)

    tracker, S
end

function track(PT::PolyhedralTracker, (cell, x∞)::Tuple{MixedCell,Vector{ComplexF64}})
    H = PT.toric_tracker.homotopy
    # The tracker works in two stages
    # 1) Revert the totric degeneration by tracking 0 to 1
    #    Here we use the following strategy
    #    a) Reparameterize the path such that lowest power of t occuring is 1
    #    b)
    #      i) If maximal power is less than 10 track from 0 to 1
    #      ii) If maximal power is larger than 10 we can get some problems for values close to 1
    #          (since t^k is << 1 for t < 1 and k large)
    #          Therefore split in two stages:
    #           1) track from 0 to 0.9
    #           2) Second, reperamerize path s.t. max power of t is 10. This shifts the
    #              problem towards 0 again (where we have more accuracy available.)
    #              And then track to 1
    # 2) Perform coefficient homotopy with possible endgame
    min_weight, max_weight =
        update_weights!(H, PT.support, PT.lifting, cell, min_weight = 1.0)

    if max_weight < 10
        retcode = track!(
            PT.toric_tracker,
            x∞,
            0.0,
            1.0;
            ω = 20.0,
            μ = 1e-12,
            max_initial_step_size = 0.2,
        )
        @unpack μ, ω = PT.toric_tracker.state
    else
        retcode = track!(
            PT.toric_tracker,
            x∞,
            0.0,
            0.9;
            ω = 20.0,
            μ = 1e-12,
            max_initial_step_size = 0.2,
        )
        @unpack μ, ω = PT.toric_tracker.state
        # TODO: check retcode?
        min_weight, max_weight =
            update_weights!(H, PT.support, PT.lifting, cell, max_weight = 10.0)

        if is_success(retcode)
            t_restart = 0.9^(1 / min_weight)
            # set min_step_size to 0.0 to not accidentally loose a solution
            min_step_size = PT.toric_tracker.options.min_step_size
            PT.toric_tracker.options.min_step_size = 0.0

            retcode = track!(
                PT.toric_tracker,
                PT.toric_tracker.state.x,
                t_restart,
                1.0;
                ω = ω,
                μ = μ,
                τ = 0.1 * t_restart,
                keep_steps = true,
            )
            @unpack μ, ω = PT.toric_tracker.state
            PT.toric_tracker.options.min_step_size = min_step_size
        end
    end
    if !is_success(retcode)
        state = PT.toric_tracker.state
        return PathResult(
            return_code = PathTrackerCode.polyhedral_failed,
            solution = copy(state.x),
            t = real(state.t),
            accuracy = state.accuracy,
            winding_number = nothing,
            last_path_point = nothing,
            valuation = nothing,
            ω = state.ω,
            μ = state.μ,
            extended_precision = state.extended_prec,
            accepted_steps = state.accepted_steps,
            rejected_steps = state.rejected_steps,
            extended_precision_used = state.used_extended_prec,
        )
    end

    r = track(PT.generic_tracker, PT.toric_tracker.state.x; ω = ω, μ = μ)
    # report accurate steps
    r.accepted_steps += PT.toric_tracker.state.accepted_steps
    r.rejected_steps += PT.toric_tracker.state.rejected_steps

    r
end
