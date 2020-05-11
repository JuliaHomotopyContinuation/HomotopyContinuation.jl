#
# We follow the Numerical differentiation scheme derived in Mackens '89

mutable struct NumericalDifferentiation
    xh::Vector{ComplexF64}
    u₁::Vector{ComplexF64}
    u₂::Vector{ComplexF64}
    u₃::Vector{ComplexF64}
    logabs_norm::Vector{Float64}
end

function NumericalDifferentiation(m::Int, n::Int)
    xh = zeros(ComplexF64, n)
    NumericalDifferentiation(xh, (zeros(ComplexF64, m) for i = 1:3)..., zeros(4))
end

function g!(u, H, tx::TaylorVector{1}, t, h, xh)
    @inbounds for i = 1:length(tx)
        xh[i] = first(tx[i])
    end
    evaluate!(u, H, xh, t + h)
end
function g!(u, H, tx::TaylorVector{2}, t, h, xh)
    @inbounds for i = 1:length(tx)
        x, x¹ = tx[i]
        xh[i] = muladd(h, x¹, x)
    end
    evaluate!(u, H, xh, t + h)
end
function g!(u, H, tx::TaylorVector{3}, t, h, xh)
    @inbounds for i = 1:length(tx)
        x, x¹, x² = tx[i]
        xh[i] = muladd(h, muladd(h, x², x¹), x)
    end
    evaluate!(u, H, xh, t + h)
end

function g!(u, H, tx::TaylorVector{4}, t, h, xh)
    @inbounds for i = 1:length(tx)
        x, x¹, x², x³ = tx[i]
        xh[i] = muladd(h, muladd(h, muladd(h, x³, x²), x¹), x)
    end
    evaluate!(u, H, xh, t + h)
end

function finite_diff!(
    u,
    H::AbstractHomotopy,
    x::TaylorVector,
    t,
    ND::NumericalDifferentiation,
    h;
    order::Int,
)
    @unpack u₁, u₂, xh = ND

    g!(u₁, H, x, t, h, xh)
    g!(u₂, H, x, t, -h, xh)

    hk = h^order
    if iseven(order)
        u .= 0.5 .* (u₁ .+ u₂) ./ hk
    else
        u .= 0.5 .* (u₁ .- u₂) ./ hk
    end
end


function finite_diff_taylor!(
    u,
    ::Val{N},
    H::AbstractHomotopy,
    x::TaylorVector{N},
    t,
    ND::NumericalDifferentiation;
    incremental::Bool = false,
    dist_to_target::Float64,
) where {N}
    a = ND.logabs_norm
    tx = vectors(x)
    if incremental
        ND.logabs_norm[N] = log(LA.norm(tx[N], InfNorm()))
    else
        for i = 1:N
            ND.logabs_norm[i] = log(LA.norm(tx[i], InfNorm()))
        end
    end

    # We minimize the sum of the squared distance to the mean of the values
    # Use not machine precision to account for some error in the evaluation
    ε = 1.4210854715202004e-14

    # λ should be a scaling factor such that the derivatives of x(t/λ) are of of the same
    # order of magnitude. We can either compute this by balancing the derivatives or
    # re-using the trust region size of the previous step
    if N == 1
        logλ = 0.0
    elseif N == 2
        logλ = a[1] - a[2]
    elseif N == 3
        logλ = 0.5 * (a[1] - a[3])
    elseif N == 4
        logλ = 0.1 * (3 * a[1] + a[2] - a[3] - 3 * a[4])
    end
    λ = exp(logλ)
    # truncation err = (1/λ)^2*h^2
    # round-off err   = ε/h^N
    # -->
    # Compute optimal ĥ by
    #      round-off err = trunc err
    # <=>  ελ^-N/ĥ^3 = λ^(-N-2)*ĥ^2
    # <=>  (ε λ^2)^(1/(N+2)) = ĥ
    ĥ = (ε * λ^2)^(1 / (N + 2))
    # Check that truncation and round-off error are both acceptable
    trunc_err = ĥ^2 / λ^2
    rnd_err = ε / ĥ^N
    if trunc_err + rnd_err > 0.05^(4 - N) || !isfinite(ĥ)
        return false
    end

    if dist_to_target < ĥ
        finite_diff!(u, H, x, t, ND, im * ĥ; order = N)
    else
        finite_diff!(u, H, x, t, ND, ĥ; order = N)
    end
    return true
end

## Default handling ignores incremental
taylor!(u, v::Val, H::AbstractHomotopy, tx::TaylorVector, t, incremental::Bool) =
    taylor!(u, v, H, tx, t)

## Type dispatch on automatic differentiation or numerical differentiation
struct AD{N} end
function AD(N::Int)
    if !(0 ≤ N ≤ 4)
        throw(ArgumentError("`automatic_differentiation` has to be between 0 and 4."))
    end
    AD{N}()
end


@generated function taylor!(
    u,
    v::Val{M},
    H,
    tx,
    t,
    ::AD{N},
    ND::NumericalDifferentiation;
    dist_to_target::Float64,
    incremental::Bool = false,
) where {M,N}
    if M ≤ N
        quote
            taylor!(u, v, H, tx, t, incremental)
            true
        end
    else
        quote
            trust = finite_diff_taylor!(
                u,
                v,
                H,
                tx,
                t,
                ND;
                incremental = M > N + 1,
                dist_to_target = dist_to_target,
            )
            trust
        end
    end
end
