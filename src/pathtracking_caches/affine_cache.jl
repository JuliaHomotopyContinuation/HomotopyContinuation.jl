struct AffineCache{T} <: AbstractPathtrackerCache{T}
    A::Matrix{T} # this is a NxN matrix
    b::Vector{T} # this is a N-Vector
end

function alg_cache(alg::AffinePredictorCorrector, H::AbstractHomotopy, x::AbstractVector{T}) where T
    n = length(x)
    A = zeros(T, n, n)
    b = zeros(T, n)
    AffineCache(A, b)
end


function perform_step!(tracker, values::PathtrackerPrecisionValues{T}, cache::AffineCache{Complex{T}}) where T
    @unpack s, ds = tracker
    @unpack H, cfg, x, xnext = values
    @unpack A, b = cache

    m = size(A,2)

    # PREDICT
    # put jacobian in A
    jacobian!(A, H, x, s, cfg)
    # put Hdt in b
    dt!(b, H, x, s, cfg, true)

    # this computes A x = b and stores the result x in b
    LU = lufact!(A)
    # there is a bug in v0.6.0 see patches.jl
    my_A_ldiv_B!(LU, b)

    xnext .= x .- ds .* b

    # CORRECT
    @unpack abstol, corrector_maxiters = tracker.options
    tracker.step_sucessfull = correct!(xnext, s + ds, H, cfg, abstol, corrector_maxiters, cache)
    nothing
end

function correct!(xnext, s, H, cfg,
    abstol::Float64,
    maxiters::Int,
    cache::AffineCache{Complex{T}}) where T
    @unpack A, b = cache
    m = size(A,2)
    k = 0
    while true
        k += 1
        evaluate!(b, H, xnext, s, cfg)

        if norm(b, Inf) < abstol
            return true
        elseif k > maxiters
            return false
        end

        # put jacobian in A
        jacobian!(A, H, xnext, s, cfg, true)

        # this computes A x = b and stores the result x in b
        LU = lufact!(A)
        # there is a bug in v0.6.0 see patches.jl
        my_A_ldiv_B!(LU, b)
        xnext .= xnext .- b
    end
end
