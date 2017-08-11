"""
    GammaTrickHomotopy(G, F[, γ])

This yields the homotopy (1-t)F + tγ(t)G, with t from 1 to 0.
The gamma trick is a special type of homotopy interpolation.
One chooses an ``γ`` uniformly from the (complex) unit circle.
This yields a more generic interpolation.

Possible to create it optionally with a seed.
"""
struct GammaTrickHomotopy{T<:Complex} <: AbstractHomotopy{T}
    start::PolySystem{T}
    target::PolySystem{T}
    γ::T

    function GammaTrickHomotopy{T}(start::PolySystem{T}, target::PolySystem{T}, γ::T) where {T<:Complex}
        N_start = nvars(start)
        N_target = nvars(target)
        n_start = nequations(start)
        n_target = nequations(target)

        if (N_start != N_target)
            return error("The number of unknowns of start and target system doesn't match." *
                "The start system of the homotopy has $N_start unknowns and the target system has $N_target unknowns.")
        end

        if (n_start != n_target)
            return error("The number of equations of start and target system doesn't match." *
                "The start system of the homotopy has $n_start equations and the target system has $n_target equations.")
        end

        new(start, target, γ)
    end
end

GammaTrickHomotopy(start::PolySystem{T}, target::PolySystem{T}, γ::T) where {T<:Complex} = GammaTrickHomotopy{T}(start, target, γ)
# GammaTrickHomotopy(start::Poly{T}, target::Poly{T}, γ::T) where {T<:Complex} = GammaTrickHomotopy([start], [target], γ)

GammaTrickHomotopy(start, target) = GammaTrickHomotopy(start, target, exp(im * (rand() * 2π - π)))
function GammaTrickHomotopy(start, target, seed::Int)
    srand(seed)
    θ = rand() * 2π - π
    GammaTrickHomotopy(start, target, exp(im * θ))
end

function evaluate(H::GammaTrickHomotopy{T}, x::Vector{T}, t::Float64) where {T<:Complex}
    (1 - t) * evaluate(H.target, x) + t * H.γ * evaluate(H.start, x)
end
function evaluate(H::GammaTrickHomotopy{T}, x::Vector{S}, t::Float64) where {T<:Complex, S<:Number}
    u = convert(Vector{T}, x)
    (1 - t) * evaluate(H.target, u) + t * H.γ * evaluate(H.start, u)
end

function jacobian(H::GammaTrickHomotopy{T}) where {T<:Number}
    J_s = jacobian(H.start)
    J_t = jacobian(H.target)

    return (x::Vector{T}, t::Float64) -> t * H.γ * J_s(x) + (1 - t) * J_t(x)
end
dt(H::GammaTrickHomotopy{T}) where {T<:Number} = (x::Vector{T}, ::Float64) -> evaluate(H.start, x) - evaluate(H.target, x)

function homogenize(H::GammaTrickHomotopy, var::MP.AbstractVariable)
    N = nvars(H)
    n = nequations(H)
    # currently square system
    if N == n
        return GammaTrickHomotopy(homogenize(H.start, var), homogenize(H.target, var), H.γ)
    elseif N == n + 1
        if is_homogenous(H.start) && is_homogenous(H.start)
           return H
        else
            return error("The Homotopy has $(N) unknowns and $(n) equations. Expected that the polynomials are already homogenous which is not the case.")
        end
    else
        return error("The Homotopy has $(N) unknowns and $(n) equations. Excepted $(n) equations and $(n) or $(n+1) unknowns.")
    end  
end
nvars(H::GammaTrickHomotopy) = nvars(H.start)
vars(H::GammaTrickHomotopy) = vars(H.start)
nequations(H::GammaTrickHomotopy) = nequations(H.start)
degrees(H::GammaTrickHomotopy) = max.(degrees(H.start), degrees(H.target))
startsystem(H::GammaTrickHomotopy) = H.start
targetsystem(H::GammaTrickHomotopy) = H.target

"""
    weyl_norm(H, t)

Computes the weyl_norm of the homotopy to the given time `t`.

## Explanation
For ``H = (1-t)F+tγG`` we have
```math
\begin{align*}
<H,H> &= <(1-t)F+tγG,(1-t)F+tγG> \\
      &= <(1-t)F,(1-t)F+tγG> + <tγG,(1-t)F+tγG> \\
      &= <(1-t)F,(1-t)F> + <(1-t)F,tγG> + <tγG,(1-t)F> + <tγG,tγG> \\
      &= <(1-t)F,(1-t)F> + 2*real(<(1-t)F,tγG>) + <tγG,tγG> \\
      &= (1-t)^2<F,F> + (1-t)t*2*real(γ)*real(<F,G>) + abs(γ)*t^2<G,G> \\
\end{align*}
```
"""
function weyl_norm(H::GammaTrickHomotopy{T}, t::Float64) where {T<:Complex}
    F = H.target
    G = H.start
    sqrt((1 - t)^2 * weyl_dot(F, F) + (1 - t) * t * 2 * real(H.γ) * real(weyl_dot(F, G)) + t^2 * abs2(H.γ) * weyl_dot(G, G))
end
