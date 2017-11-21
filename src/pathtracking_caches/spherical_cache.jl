mutable struct SphericalCache{T, aType<:AbstractMatrix{T}, asub, bType<:AbstractVector{T}, bsub} <: AbstractPathtrackerCache{T}
    A::aType # this is a NxN matrix
    # we store the views since the allocate currently
    A_sub::asub # this is the Nx(N-1) submatrix where the jacobian is stored
    # we store the views since the allocate currently
    b::bType # this is a N-Vector
    b_sub::bsub # this is the (N-1) subvector to store H(x,t) and Hdt(x,t)
end

function alg_cache(alg::SphericalPredictorCorrector, H::AbstractHomotopy, x::AbstractVector{T}) where T
    n = length(x)
    A = zeros(x, T, n, n)
    A_sub = @view A[1:end-1,:]
    b = zeros(x, T, n)
    b_sub = @view b[1:end-1]
    SphericalCache{T, typeof(A), typeof(A_sub), typeof(b), typeof(b_sub)}(A, A_sub, b, b_sub)
end

function precondition!(tracker, values::PathtrackerPrecisionValues{T}, cache::SphericalCache{Complex{T}}) where T
    normalize!(values.x)
    nothing
end

function perform_step!(tracker, values::PathtrackerPrecisionValues{T}, cache::SphericalCache{Complex{T}}) where T
    @unpack s, ds = tracker
    @unpack H, cfg, x, xnext = values
    @unpack A, A_sub, b, b_sub = cache

    m = size(A,2)

    # PREDICT
    # put jacobian in A
    jacobian!(A_sub, H, x, s, cfg)
    for j=1:m
        A[end, j] = conj(x[j])
    end
    # put Hdt in b
    dt!(b_sub, H, x, s, cfg, true)
    b[end] = zero(T)

    # this computes A x = b and stores the result x in b
    LU = lufact!(A)
    # there is a bug in v0.6.0 see patches.jl
    my_A_ldiv_B!(LU, b)

    xnext .= x .- ds .* b
    normalize!(xnext)

    # CORRECT
    @unpack abstol, corrector_maxiters = tracker.options
    s += ds

    tracker.step_sucessfull = correct!(xnext, s, H, cfg, abstol, corrector_maxiters, cache)
    nothing
end

function correct!(xnext,
    s,
    H,
    cfg,
    abstol::Float64,
    maxiters::Int,
    cache::SphericalCache{Complex{T}}) where T
    @unpack A, A_sub, b, b_sub = cache

    m = size(A,2)
    k = 0
    while true
        k += 1
        evaluate!(b_sub, H, xnext, s, cfg)
        b[end] = zero(T)

        if norm(b, Inf) < abstol
            return true
        elseif k > maxiters
            return false
        end

        # put jacobian in A
        jacobian!(A_sub, H, xnext, s, cfg, true)
        for j=1:m
            A[end, j] = conj(xnext[j])
        end
        # this computes A x = b and stores the result x in b
        LU = lufact!(A)
        # there is a bug in v0.6.0 see patches.jl
        my_A_ldiv_B!(LU, b)
        xnext .= xnext .- b
        normalize!(xnext)
    end
end


function setup_workers(cache::SphericalCache)
    # We need this to get pmap working:
    # Otherwise the references A_sub and b_sub would point
    # to the initial array and we would get singular exceptions since
    # the top of A and b will not get filled
    #
    # We can get rid of all of that if @view doesn't allocate anymore
    cache.A_sub = @view cache.A[1:end-1,:]
    cache.b_sub = @view cache.b[1:end-1]
    nothing
end
