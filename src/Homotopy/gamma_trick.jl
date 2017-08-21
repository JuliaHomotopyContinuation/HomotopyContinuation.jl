"""
    GammaTrickHomotopy(start, target [, γ])
    GammaTrickHomotopy(start, target [, seed])

Construct the homotopy (1-t)target + tγ(t)start, with t from 1 to 0.
The gamma trick is a special type of homotopy interpolation.
One chooses an ``γ`` uniformly from the (complex) unit circle.
This yields a more generic interpolation.

Possible to create it optionally with a seed.
"""
struct GammaTrickHomotopy{T<:Number, Start<:AbstractPolySystem{T}, Target<:AbstractPolySystem{T}} <: AbstractHomotopy{T}
    start::Start
    target::Target
    γ::Complex128

    function GammaTrickHomotopy{T, Start, Target}(start::AbstractPolySystem{T}, target::AbstractPolySystem{T}, γ::Complex128) where {T<:Number, Start<:AbstractPolySystem{T}, Target<:AbstractPolySystem{T}}
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

        new(start, target, γ)
    end
end

function GammaTrickHomotopy(start::Start,target::Target, γ::Complex128) where {T<:Number, Start<:AbstractPolySystem{T}, Target<:AbstractPolySystem{T}}
    GammaTrickHomotopy{T, Start, Target}(start,target,γ)
end
# promote if they don't have the same coefficients
function GammaTrickHomotopy(start::AbstractPolySystem,target::AbstractPolySystem, γ::Complex128)
    s, t = promote(start, target)
    GammaTrickHomotopy(s, t, γ)
end
GammaTrickHomotopy(s, t, γ::Complex128) = GammaTrickHomotopy(PolySystem(s), PolySystem(t), γ)

GammaTrickHomotopy(start, target) = GammaTrickHomotopy(start, target, exp(im * (rand() * 2π - π)))
function GammaTrickHomotopy(start, target, seed::Int)
    srand(seed)
    θ = rand() * 2π - π
    GammaTrickHomotopy(start, target, exp(im * θ))
end

function Base.show(io::IO, H::GammaTrickHomotopy)
    print(io, typeof(H), ":\n")
    print(io, "* γ: ", H.γ, "\n")
    println(io, "* start:")
    println(io, H.start)
    println(io, "* target:")
    print(io, H.target)
end


function evaluate(H::GammaTrickHomotopy, x, t)
    (1 - t) * evaluate(H.target, x) + t * H.γ * evaluate(H.start, x)
end
(H::GammaTrickHomotopy)(x,t) = evaluate(H,x,t)

function substitute(H::GammaTrickHomotopy, pair::Pair{Symbol,<:Number})
    GammaTrickHomotopy(substitute(H.start, pair), substitute(H.target, pair), H.γ)
end

"""
    differentiate(H::GammaTrickHomotopy)

Differentiates `H` w.r.t. the place and returns an evaluation function
`(x,t) -> t * H.γ * J_s(x) + (1 - t) * J_t(x)`
where `J_s` is the jacobian of the start system and
`J_t` is the jacobian of the target system.
"""
function differentiate(H::GammaTrickHomotopy)
    J_s = differentiate(H.start)
    J_t = differentiate(H.target)

    (x, t) -> t .* H.γ .* J_s(x) .+ (1 - t) .* J_t(x)
end

"""
    dt(H::GammaTrickHomotopy)

Differentiates `H` w.r.t. the the time and returns an evaluation function
`(x,t) -> T`
"""
dt(H::GammaTrickHomotopy) = (x, ::Number) -> H.γ * evaluate(H.start, x) .- evaluate(H.target, x)

"""
    homogenize(H::GammaTrickHomotopy)

Make the start and target system of `H` homogenous.
"""
function homogenize(H::GammaTrickHomotopy)
    GammaTrickHomotopy(homogenize(H.start), homogenize(H.target), H.γ)
end

"""
    homogenized(H::GammaTrickHomotopy)

Checks whether `H` was homogenized.
"""
homogenized(H::GammaTrickHomotopy) = homogenized(H.start)

ishomogenous(H::GammaTrickHomotopy) = ishomogenous(H.start) && ishomogenous(H.target)

"""
    nvariables(H::GammaTrickHomotopy)

The number variables of `H`.
"""
function nvariables(H::GammaTrickHomotopy)
    # Constructor guarantees that start and target have the same number of variables
    nvariables(H.start)
end
function length(H::GammaTrickHomotopy)
    # Constructor guarantees that start and target have the same length
    length(H.start)
end
"""
    degrees(H::GammaTrickHomotopy)

The (total) degrees of the polynomials of the system `P`.
"""
degrees(H::GammaTrickHomotopy) = max.(degrees(H.start), degrees(H.target))
"""
    startsystem(H::GammaTrickHomotopy)

The startsystem of `H`.
"""
startsystem(H::GammaTrickHomotopy) = H.start
"""
    targetsystem(H::GammaTrickHomotopy)

The targetsystem of `H`.
"""
targetsystem(H::GammaTrickHomotopy) = H.target


"""
    weylnorm(H, t)

Computes the weylnorm of the homotopy to the given time `t`.

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
function weylnorm(H::GammaTrickHomotopy{T}, t::Float64) where {T<:Complex}
    F = H.target
    G = H.start
    sqrt((1 - t)^2 * weyldot(F, F) + (1 - t) * t * 2 * real(H.γ) * real(weyldot(F, G)) + t^2 * abs2(H.γ) * weyldot(G, G))
end
