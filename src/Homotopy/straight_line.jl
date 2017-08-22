"""
    StraightLineHomotopy(start::AbstractPolySystem, target::AbstractPolySystem)

Constructs the homotopy `t * start + (1-t) * target`.
"""
struct StraightLineHomotopy{T<:Number, Start<:AbstractPolySystem{T}, Target<:AbstractPolySystem{T}} <: AbstractHomotopy{T}
    start::Start
    target::Target

    function StraightLineHomotopy{T, Start, Target}(start::Start, target::Target) where {T<:Number, Start<:AbstractPolySystem{T}, Target<:AbstractPolySystem{T}}
        if (homogenized(start) != homogenized(target))
            return error("start and target have to be both either homogenized or not")
        end
        if (nvariables(start) != nvariables(target))
            return error("The number of unknowns of start and target system doesn't match." *
                "The start system of the homotopy has $(nvariables(start)) unknowns " *
                "and the target system has $(nvariables(target)) unknowns.")
        end

        if (length(target) != length(start))
            return error("The number of equations of start and target system doesn't match." *
                "The start system of the homotopy has $(length(start)) equations " *
                "and the target system has $(length(target)) equations.")
        end

        new(start, target)
    end
end

function StraightLineHomotopy(start::Start,target::Target) where {T<:Number, Start<:AbstractPolySystem{T}, Target<:AbstractPolySystem{T}}
    StraightLineHomotopy{T, Start, Target}(start,target)
end
# promote if they don't have the same coefficients
function StraightLineHomotopy(start::AbstractPolySystem,target::AbstractPolySystem)
    s, t = promote(start, target)
    StraightLineHomotopy(s, t)
end
StraightLineHomotopy(s, t) = StraightLineHomotopy(PolySystem(s), PolySystem(t))

function evaluate(H::StraightLineHomotopy, x, t)
    (1-t) .* evaluate(H.target, x) .+ t .* evaluate(H.start, x)
end
(H::StraightLineHomotopy)(x,t) = evaluate(H,x,t)

function substitute(H::StraightLineHomotopy, pair::Pair{Symbol,<:Number})
    StraightLineHomotopy(substitute(H.start, pair), substitute(H.target, pair))
end

function removepoly(H::StraightLineHomotopy, i::Int)
    StraightLineHomotopy(removepoly(H.start,i), removepoly(H.target,i))
end

function Base.show(io::IO, H::StraightLineHomotopy)
    print(io, typeof(H), ":\n")
    println(io, "* start:")
    println(io, H.start)
    println(io, "* target:")
    print(io, H.target)
end
"""
    differentiate(H::StraightLineHomotopy)

Differentiates `H` w.r.t. the place and returns an evaluation function
`(x,t) -> t * J_s(x) + (1 - t) * J_t(x)`
where `J_s` is the jacobian of the start system and
`J_t` is the jacobian of the target system.
"""
function differentiate(H::StraightLineHomotopy)
    J_s = differentiate(H.start)
    J_t = differentiate(H.target)

    (x, t) -> t .* J_s(x) .+ (1 - t) .* J_t(x)
end
"""
    dt(H::StraightLineHomotopy)

Differentiates `H` w.r.t. the the time and returns an evaluation function
`(x,t) -> T`
"""
dt(H::StraightLineHomotopy) = (x, ::Number) -> evaluate(H.start, x) .- evaluate(H.target, x)

"""
    homogenize(H::StraightLineHomotopy)

Make the start and target system of `H` homogenous.
"""
function homogenize(H::StraightLineHomotopy)
    StraightLineHomotopy(homogenize(H.start), homogenize(H.target))
end

"""
    homogenized(H::StraightLineHomotopy)

Checks whether `H` was homogenized.
"""
homogenized(H::StraightLineHomotopy) = homogenized(H.start)

ishomogenous(H::StraightLineHomotopy) = ishomogenous(H.start) && ishomogenous(H.target)

"""
    nvariables(H::StraightLineHomotopy)

The number variables of `H`.
"""
function nvariables(H::StraightLineHomotopy)
    # Constructor guarantees that start and target have the same number of variables
    nvariables(H.start)
end
function length(H::StraightLineHomotopy)
    # Constructor guarantees that start and target have the same length
    length(H.start)
end
"""
    degrees(H::StraightLineHomotopy)

The (total) degrees of the polynomials of the system `P`.
"""
degrees(H::StraightLineHomotopy) = max.(degrees(H.start), degrees(H.target))
"""
    startsystem(H::StraightLineHomotopy)

The startsystem of `H`.
"""
startsystem(H::StraightLineHomotopy) = H.start
"""
    targetsystem(H::StraightLineHomotopy)

The targetsystem of `H`.
"""
targetsystem(H::StraightLineHomotopy) = H.target

"""
    weylnorm(H, t)

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
function weylnorm(H::StraightLineHomotopy{T}, t::Float64) where {T<:Complex}
    F = H.target
    G = H.start

    sqrt((1-t)^2 * weyldot(F,F) + (1-t)*t*2*real(weyldot(F,G)) + t^2 * weyldot(G,G))
end
