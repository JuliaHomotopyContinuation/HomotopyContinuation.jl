export ParameterHomotopy, nparameters, set_parameters!

"""
    ParameterHomotopy(F, parameters;
        variables=setdiff(MP.variables(F), parameters),
        p₁=randn(ComplexF64, length(parameters)),
        p₀=randn(ComplexF64, length(parameters)),
        γ₁=nothing, γ₀=nothing)

Construct the homotopy
```math
H(x, t) = F(x, (tγ₁p₁+(1-t)γ₀p₀) / (tγ₁+(1-t)γ₀)),
```
where `p₁` and `p₀` are a vector of parameter values for ``F`` and
`γ₁` and `γ₀` are complex numbers. If `γ₁` or `γ₀` is `nothing`, it is assumed
that `γ₁` and `γ₀` are ``1``.
The input ``parameters`` specifies the parameter variables of ``F``.
Neccessarily, `length(parameters) == length(p₁) == length(p₀)`.

Note that `p₁` and `p₀` are stored as a tuple `p` of `SVectors` and `γ₁` and `γ₀`
are stored as a tuple `γ` or as `γ=nothing`

    ParameterHomotopy(F, parameters;
        variables=setdiff(MP.variables(F), parameters),
        start_parameters=randn(ComplexF64, length(parameters)),
        target_parameters=randn(ComplexF64, length(parameters)),
        start_gamma=nothing, target_gamma=nothing)

This is a non-unicode variant where `γ₁=start_parameters`, `γ₀=target_parameters`,
`γ₁=start_gamma`, `γ₀=target_gamma`.
"""
mutable struct ParameterHomotopy{Sys<:AbstractSystem, T<:Number} <: AbstractHomotopy
    F::Sys
    p::NTuple{2, Vector{T}}
    γ::Union{Nothing, NTuple{2, ComplexF64}}
end

function ParameterHomotopy(F::AbstractSystem, p₁::AbstractVector, p₀::AbstractVector;
    start_gamma=nothing, γ₀ = start_gamma,
    target_gamma=nothing, γ₁ = target_gamma)

    length(p₁) == length(p₀) || error("Length of parameters provided doesn't match.")

    γ = (γ₁ === nothing || γ₀ === nothing) ? nothing : (γ₁, γ₀)
    p = Vector.(promote(p₁, p₀))

    ParameterHomotopy(F, p, γ)
end

function ParameterHomotopy(F::AbstractSystem;
    start_parameters=nothing, p₁ = start_parameters,
    target_parameters=nothing, p₀ = target_parameters,
    kwargs...)
    (p₁ !== nothing && p₀ !== nothing) || error("Parameters not provided for ParameterHomotopy")
    ParameterHomotopy(F, p₁, p₀; kwargs...)
end

function ParameterHomotopy(F::Vector{T},
    parameters::AbstractVector{V};
    variables=setdiff(MP.variables(F), parameters),
    start_parameters=randn(ComplexF64, length(parameters)),
    target_parameters=randn(ComplexF64, length(parameters)),
    kwargs...) where {T<:MP.AbstractPolynomialLike, V<:MP.AbstractVariable}
    G = SPSystem(F; variables=variables, parameters=parameters)
    ParameterHomotopy(G; start_parameters=start_parameters,
        target_parameters=target_parameters, kwargs...)
end

struct ParameterHomotopyCache{C<:AbstractSystemCache, T1<:Number, T2<:Number} <: AbstractHomotopyCache
    F_cache::C
    pt::Vector{T1}
    ∂p∂t::Vector{T1}
    J_p::Matrix{T2}
end

(H::ParameterHomotopy)(x, t, c=cache(H, x, t)) = evaluate(H, x, t, c)

function cache(H::ParameterHomotopy, x, t)
    if H.γ === nothing
        pt = Vector{typeof(t * H.p[1][1])}(undef, length(H.p[1]))
    else
        pt = Vector{typeof(H.γ[1] * t * H.p[1][1])}(undef, length(H.p[1]))
    end
    p!(pt, H, t)
    F_cache = cache(H.F, x, pt)
    ∂p∂t = Vector{eltype(pt)}(undef, length(pt))
    J_p = Matrix(differentiate_parameters(H.F, x, pt, F_cache))

    ParameterHomotopyCache(F_cache, pt, ∂p∂t, J_p)
end

Base.size(H::ParameterHomotopy) = size(H.F)

"""
    nparameters(H::ParameterHomotopy)

Returns the number of parameters of `H`.
"""
nparameters(H::ParameterHomotopy) = length(H.p[1])


"""
    set_parameters!(H::ParameterHomotopy, p::Tuple, γ)

Update the parameters `p` and `γ` of `H`.
"""
function set_parameters!(H::ParameterHomotopy, p::Tuple, γ=nothing)
    H.p[1] .= p[1]
    H.p[2] .= p[2]
    H.γ = γ
    H
end

"""
    set_start_parameters!(H::ParameterHomotopy, p)

Update the start parameters of `H`.
"""
function set_start_parameters!(H::ParameterHomotopy, p)
    H.p[1] .= p
    H
end

"""
    set_target_parameters!(H::ParameterHomotopy, p)

Update the target parameters of `H`.
"""
function set_target_parameters!(H::ParameterHomotopy, p)
    H.p[2] .= p
    H
end

"""
    set_parameters!(H::ParameterHomotopy, p₁, p₀, γ)

Update the parameters `p` and `γ` of `H`.
"""
function set_parameters!(H::ParameterHomotopy, p₁::AbstractVector, p₀::AbstractVector, γ=nothing)
    set_parameters!(H, (p₁, p₀), γ)
end


@inline function p!(pt, H::ParameterHomotopy, t)
    p₁, p₀ = H.p
    if H.γ === nothing
        for i in eachindex(pt)
            @inbounds pt[i] = t * p₁[i] + (1 - t) * p₀[i]
        end
    else
        # compute (tγ₁p₁+(1-t)γ₀p₀) / (tγ₁+(1-t)γ₀)
        γ₁, γ₀ = H.γ
        tγ₁, γ₀_₁₋t = t * γ₁, (1 - t) * γ₀
        γ = (tγ₁ + γ₀_₁₋t)
        a = (@fastmath tγ₁ / γ)
        b = (@fastmath γ₀_₁₋t / γ)
        for i in eachindex(pt)
            @inbounds pt[i] = a * p₁[i] + b * p₀[i]
        end
    end
    pt
end

@inline function ∂p∂t!(u, H::ParameterHomotopy, t, c::ParameterHomotopyCache)
    p₁, p₀ = H.p
    if H.γ === nothing
        for i in eachindex(p₁)
            u[i] = p₁[i] - p₀[i]
        end
    else
        γ₁, γ₀ = H.γ
        tγ₁, γ₀_₁₋t = t * γ₁, (1 - t) * γ₀
        γ = (tγ₁ + γ₀_₁₋t)
        λ = @fastmath γ₁ * γ₀ / (γ * γ)
        for i in eachindex(p₁)
            u[i] = λ * (p₁[i] - p₀[i])
        end
    end
end


function evaluate!(u, H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    evaluate!(u, H.F, x, p!(c.pt, H, t), c.F_cache)
end
function evaluate(H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    evaluate(H.F, x, p!(c.pt, H, t), c.F_cache)
end

function jacobian!(u, H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    jacobian!(u, H.F, x, p!(c.pt, H, t), c.F_cache)
end
function jacobian(H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    jacobian(H.F, x, p!(c.pt, H, t), c.F_cache)
end

function evaluate_and_jacobian!(u, U, H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    evaluate_and_jacobian!(u, U, H.F, x, p!(c.pt, H, t), c.F_cache)
end

function dt!(u, H::ParameterHomotopy, x, t, c::ParameterHomotopyCache)
    # apply chain rule to H(x, p(t))
    p!(c.pt, H, t)
    ∂p∂t!(c.∂p∂t, H, t, c)
    differentiate_parameters!(c.J_p, H.F, x, c.pt, c.F_cache)
    LinearAlgebra.mul!(u, c.J_p, c.∂p∂t)
end
