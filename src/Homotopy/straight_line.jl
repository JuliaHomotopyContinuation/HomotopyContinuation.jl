"""
```
StraightLineHomotopy(G, F)
```

The homotopy ``(1-t)F + tG``.
"""
struct StraightLineHomotopy{T<:Number} <: AbstractHomotopy{T}
    start::Vector{Poly{T}}
    target::Vector{Poly{T}}

    function StraightLineHomotopy{T}(start::Vector{Poly{T}}, target::Vector{Poly{T}}) where {T<:Number}
        N_start = nvars(start)
        N_target = nvars(target)
        n_start = nequations(start)
        n_target = nequations(target)

        if (N_start != N_target)
            return error( "The number of unknowns of start and target system doesn't match." *
                "The start system of the homotopy has $N_start unknowns and the target system has $N_target unknowns.")
        end

        if (n_start != n_target)
            return error( "The number of equations of start and target system doesn't match." *
                "The start system of the homotopy has $n_start equations and the target system has $n_target equations.")
        end

        new(start, target)
    end
end

function StraightLineHomotopy(start::Vector{Poly{T}},target::Vector{Poly{T}}) where {T<:Number}
    StraightLineHomotopy{T}(start,target)
end

function StraightLineHomotopy(start::Poly{T},target::Poly{T}) where {T<:Number}
    StraightLineHomotopy([start], [target])
end

function evaluate(H::StraightLineHomotopy{T}, x::Vector{T}, t::Float64) where {T<:Number}
    (1-t) * evaluate(H.target, x) + t * evaluate(H.start, x)
end

function jacobian(H::StraightLineHomotopy{T}) where {T<:Number}
    J_s = jacobian(H.start)
    J_t = jacobian(H.target)

    return (x::Vector{T}, t::Float64) -> t * J_s(x) + (1 - t) * J_t(x)
end
dt(H::StraightLineHomotopy{T}) where {T<:Number} = (x::Vector{T}, ::Float64) -> evaluate(H.start, x) - evaluate(H.target, x)

function homogenize(H::StraightLineHomotopy)
    N = nvars(H)
    n = nequations(H)
    # currently square system
    if N == n
        return StraightLineHomotopy(homogenize(H.start), homogenize(H.target))
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
nvars(H::StraightLineHomotopy) = nvars(H.start)
nequations(H::StraightLineHomotopy) = nequations(H.start)
degrees(H::StraightLineHomotopy) = max.(degrees(H.start), degrees(H.target))
startsystem(H::StraightLineHomotopy) = H.start
targetsystem(H::StraightLineHomotopy) = H.target

"""
    weyl_norm(H, t)

Computes the weyl_norm of the homotopy to the given time `t`.

## Explanation
For ``H = (1-t)F+tG`` we have
```math
\begin{align*}
<H,H> &= <(1-t)F+tG,(1-t)F+tG> \\
      &= <(1-t)F,(1-t)F+tG> + <tG,(1-t)F+tG> \\
      &= <(1-t)F,(1-t)F> + <(1-t)F,tG> + <tG,(1-t)F> + <tG,tG> \\
      &= <(1-t)F,(1-t)F> + 2*Real(<(1-t)F,tG>) + <tG,tG> \\
      &= (1-t)^2<F,F> + (1-t)t*2*Real(<F,G>) + t^2<G,G> \\
\end{align*}
```
"""
function weyl_norm(H::StraightLineHomotopy{T}, t::Float64) where {T<:Complex}
    F = H.target
    G = H.start

    sqrt((1-t)^2 * weyl_dot(F,F) + (1-t)*t*2*real(weyl_dot(F,G)) + t^2 * weyl_dot(G,G))
end
