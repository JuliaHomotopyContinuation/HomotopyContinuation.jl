export StraightLineHomotopy

"""
    StraightLineHomotopy(G, F; gamma=exp(i * 2π*rand()))
Construct the homotopy ``H(x, t) = γtG(x) + (1-t)F(x)``.
"""
struct StraightLineHomotopy{S,T,P1,P2} <: AbstractHomotopy
    start::ModelKit.CompiledSystem{S}
    target::ModelKit.CompiledSystem{T}

    start_parameters::Vector{P1}
    target_parameters::Vector{P2}

    u::Vector{ComplexF64}
    ū::Vector{ComplexDF64}
    U::Matrix{ComplexF64}

    dv_start::NTuple{5,Vector{ComplexF64}}
    dv_target::NTuple{5,Vector{ComplexF64}}
end

function StraightLineHomotopy(
    start,
    target;
    start_parameters = ComplexF64[],
    target_parameters = ComplexF64[],
)
    size(start) == size(target) ||
    throw(ArgumentError("Start and target do not have the same size, got $(size(start)) and $(size(target))"))

    m, n = size(start)
    u = zeros(ComplexF64, m)
    ū = zeros(ComplexDF64, m)
    U = zeros(ComplexF64, m, n)

    dv_start = tuple((zeros(ComplexF64, m) for i = 0:4)...)
    dv_target = tuple((zeros(ComplexF64, m) for i = 0:4)...)

    StraightLineHomotopy(
        start,
        target,
        copy(start_parameters),
        copy(target_parameters),
        u,
        ū,
        U,
        dv_start,
        dv_target,
    )
end

Base.size(H::StraightLineHomotopy) = size(H.start)

function evaluate!(u, H::StraightLineHomotopy, x::AbstractVector{T}, t) where {T}
    ModelKit.evaluate!(u, H.start, x, H.start_parameters)

    if T isa ComplexDF64 || T isa DoubleF64
        ModelKit.evaluate!(H.ū, H.target, x, H.target_parameters)
        @inbounds u .= t .* u .+ (1.0 .- t) .* H.ū
    else
        ModelKit.evaluate!(H.u, H.target, x, H.target_parameters)
        @inbounds u .= t .* u .+ (1.0 .- t) .* H.u
    end
    u
end

function evaluate_and_jacobian!(
    u,
    U,
    H::StraightLineHomotopy,
    x::AbstractVector{T},
    t,
) where {T}
    ModelKit.evaluate_and_jacobian!(u, U, H.start, x, H.start_parameters)

    if T isa ComplexDF64 || T isa DoubleF64
        ModelKit.evaluate_and_jacobian!(H.ū, H.U, H.target, x, H.target_parameters)
        @inbounds u .= t .* u .+ (1.0 .- t) .* H.ū
    else
        ModelKit.evaluate_and_jacobian!(H.u, H.U, H.target, x, H.target_parameters)
        @inbounds u .= t .* u .+ (1.0 .- t) .* H.u
    end
    @inbounds U .= t .* U .+ (1.0 .- t) .* H.U
    nothing
end

function diff_t!(u, H::StraightLineHomotopy, x, t, ::Tuple{})
    ModelKit.evaluate!(u, H.start, x, H.start_parameters)
    ModelKit.evaluate!(H.u, H.target, x, H.target_parameters)
    @inbounds u .= u .- H.u
end
function diff_t!(u, H::StraightLineHomotopy, x, t, dx::Tuple)
    ModelKit.taylor!(H.dv_start, H.start, x, dx, H.start_parameters)
    ModelKit.taylor!(H.dv_target, H.target, x, dx, H.target_parameters)

    k = length(dx) + 1
    @inbounds v_start, v_start′ = H.dv_start[k], H.dv_start[k+1]
    @inbounds v_target, v_target′ = H.dv_target[k], H.dv_target[k+1]
    @inbounds for i = 1:size(H, 1)
        u[i] = t * v_start′[i] + v_start[i] + (1 - t) * v_target′[i] - v_target[i]
    end
    u
end
