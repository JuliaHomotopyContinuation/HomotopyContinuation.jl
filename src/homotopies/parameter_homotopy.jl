export ParameterHomotopy

"""
    ParameterHomotopy(F::Union{AbstractSystem,System}; start_parameters, target_parameters)
    ParameterHomotopy(F::Union{AbstractSystem,System}, start_parameters, target_parameters)

Construct the parameter homotopy ``H(x,t) = F(x; t p + (1 - t) q)`` where ``p`` is
`start_parameters` and ``q`` is `target_parameters`.
"""
struct ParameterHomotopy{T<:AbstractSystem} <: AbstractHomotopy
    F::T
    p::Vector{ComplexF64}
    q::Vector{ComplexF64}
    #cache
    t_cache::Base.RefValue{ComplexF64}
    taylor_pt::TaylorVector{2,ComplexF64}
end

function ParameterHomotopy(
    F;
    start_parameters::AbstractVector,
    target_parameters::AbstractVector,
)
    ParameterHomotopy(F, start_parameters, target_parameters)
end
function ParameterHomotopy(F::ModelKit.System, p::AbstractVector, q::AbstractVector)
    ParameterHomotopy(ModelKitSystem(F), p, q)
end
function ParameterHomotopy(F::AbstractSystem, p::AbstractVector, q::AbstractVector)
    @assert length(p) == length(q) == nparameters(F)

    p̂ = Vector{ComplexF64}(p)
    q̂ = Vector{ComplexF64}(q)
    taylor_pt = TaylorVector{2}(ComplexF64, length(q))

    ParameterHomotopy(F, p̂, q̂, Ref(complex(NaN)), taylor_pt)
end

Base.size(H::ParameterHomotopy) = size(H.F)

function start_parameters!(H::ParameterHomotopy, p)
    H.p .= p;
    # void cache
    H.t_cache[] = NaN
    H
end
function target_parameters!(H::ParameterHomotopy, q)
    H.q .= q
    H.t_cache[] = NaN
    H
end

function tp!(H::ParameterHomotopy, t::Union{ComplexF64,Float64})
    t == H.t_cache[] && return H.taylor_pt

    if imag(t) == 0
        let t = real(t)
            @inbounds for i = 1:length(H.taylor_pt)
                H.taylor_pt[i] = (t * H.p[i] + (1.0 - t) * H.q[i], H.p[i] - H.q[i])
            end
        end
    else
        @inbounds for i = 1:length(H.taylor_pt)
            H.taylor_pt[i] = (t * H.p[i] + (1.0 - t) * H.q[i], H.p[i] - H.q[i])
        end
    end
    H.t_cache[] = t

    H.taylor_pt
end

function evaluate!(u, H::ParameterHomotopy, x, t)
    evaluate!(u, H.F, x, first(vectors(tp!(H, t))))
end

function evaluate_and_jacobian!(u, U, H::ParameterHomotopy, x, t)
    evaluate_and_jacobian!(u, U, H.F, x, first(vectors(tp!(H, t))))
end

function taylor!(u, v::Val, H::ParameterHomotopy, tx, t)
    taylor!(u, v, H.F, tx, tp!(H, t))
    u
end
