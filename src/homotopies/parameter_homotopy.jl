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
    pt::Vector{ComplexF64}
    taylor_pt::TaylorVector{2,ComplexF64}
    t_cache_extended::Base.RefValue{ComplexDF64}
    pt_extended::Vector{ComplexDF64}
    taylor_pt_extended::TaylorVector{2,ComplexDF64}
end

function ParameterHomotopy(
    F;
    start_parameters::AbstractVector,
    target_parameters::AbstractVector,
)
    ParameterHomotopy(F, start_parameters, target_parameters)
end
function ParameterHomotopy(
    F::ModelKit.System,
    p::AbstractVector,
    q::AbstractVector;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
)
    ParameterHomotopy(fixed(F; compile = compile), p, q)
end
function ParameterHomotopy(F::AbstractSystem, p::AbstractVector, q::AbstractVector)
    @assert length(p) == length(q) == nparameters(F)

    p̂ = Vector{ComplexF64}(p)
    q̂ = Vector{ComplexF64}(q)
    taylor_pt = TaylorVector{2}(ComplexF64, length(q))
    pt = copy(p̂)

    ParameterHomotopy(
        F,
        p̂,
        q̂,
        Ref(complex(NaN)),
        pt,
        taylor_pt,
        Ref(ComplexDF64(NaN)),
        ComplexDF64.(p̂),
        TaylorVector{2}(ComplexDF64, length(q)),
    )
end

Base.size(H::ParameterHomotopy) = size(H.F)

function start_parameters!(H::ParameterHomotopy, p)
    H.p .= p
    # void cache
    H.t_cache[] = NaN
    H.t_cache_extended[] = ComplexDF64(NaN)
    H
end
function target_parameters!(H::ParameterHomotopy, q)
    H.q .= q
    H.t_cache[] = NaN
    H.t_cache_extended[] = ComplexDF64(NaN)
    H
end
function parameters!(H::ParameterHomotopy, p, q)
    H.p .= p
    H.q .= q
    H.t_cache[] = NaN
    H.t_cache_extended[] = ComplexDF64(NaN)
    H
end

function tp!(H::ParameterHomotopy, t::Union{ComplexF64,Float64})
    t == H.t_cache[] && return H.taylor_pt

    if imag(t) == 0
        let t = real(t)
            @inbounds for i = 1:length(H.taylor_pt)
                ptᵢ = t * H.p[i] + (1.0 - t) * H.q[i]
                H.pt[i] = ptᵢ
                H.taylor_pt[i] = (ptᵢ, H.p[i] - H.q[i])
            end
        end
    else
        @inbounds for i = 1:length(H.taylor_pt)
            ptᵢ = t * H.p[i] + (1.0 - t) * H.q[i]
            H.pt[i] = ptᵢ
            H.taylor_pt[i] = (ptᵢ, H.p[i] - H.q[i])
        end
    end
    H.t_cache[] = t

    H.taylor_pt
end

function ModelKit.evaluate!(u, H::ParameterHomotopy, x, t)
    T = promote_type(ModelKit._interpreter_precision(x), ModelKit._interpreter_precision(t))
    tp!(H, t, T)
    evaluate!(u, H.F, x, T === ComplexDF64 ? H.pt_extended : H.pt)
end

function ModelKit.evaluate_and_jacobian!(u, U, H::ParameterHomotopy, x, t)
    T = promote_type(ModelKit._interpreter_precision(x), ModelKit._interpreter_precision(t))
    tp!(H, t, T)
    evaluate_and_jacobian!(u, U, H.F, x, T === ComplexDF64 ? H.pt_extended : H.pt)
end

function ModelKit.taylor!(u, v::Val, H::ParameterHomotopy, tx, t)
    T = promote_type(
        ModelKit._interpreter_precision(tx),
        ModelKit._interpreter_precision(t),
    )
    taylor!(u, v, H.F, tx, tp!(H, t, T))
    u
end
function ModelKit.taylor!(u, v::Val, H::ParameterHomotopy, tx, t, incremental::Bool)
    T = promote_type(
        ModelKit._interpreter_precision(tx),
        ModelKit._interpreter_precision(t),
    )
    ModelKit.predictor_taylor!(u, v, H.F, tx, tp!(H, t, T), incremental)
end

tp!(H::ParameterHomotopy, t, ::Type{ComplexF64}) = tp!(H, ComplexF64(t))
function tp!(H::ParameterHomotopy, t, ::Type{ComplexDF64})
    τ = ComplexDF64(t)
    τ == H.t_cache_extended[] && return H.taylor_pt_extended
    from_target = abs2(ComplexF64(t)) <= abs2(ComplexF64(t) - 1)
    @inbounds for i in eachindex(H.p)
        p, q = ComplexDF64(H.p[i]), ComplexDF64(H.q[i])
        d = p - q
        value = from_target ? q + τ * d : p + (τ - 1) * d
        H.pt_extended[i] = value
        H.taylor_pt_extended[i] = (value, d)
    end
    H.t_cache_extended[] = τ
    H.taylor_pt_extended
end
