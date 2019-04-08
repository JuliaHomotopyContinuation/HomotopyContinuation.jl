export HomotopyNullCache,
    nvariables,
    cache,
    evaluate!,
    evaluate,
    jacobian!, jacobian,
    evaluate_and_jacobian!, evaluate_and_jacobian,
    dt!, dt,
    jacobian_and_dt!, jacobian_and_dt,
    basehomotopy

# Cache


"""
    HomotopyNullCache

The default `AbstractHomotopyCache` containing nothing.
"""
struct HomotopyNullCache <: AbstractHomotopyCache end
const HNC = HomotopyNullCache

"""
    cache(H::AbstractHomotopy, x, t)::AbstractHomotopyCache

Create a cache for the evaluation (incl. Jacobian) of `F` with elements of the type
of `x`. The default implementation returns [`HomotopyNullCache`](@ref).
"""
cache(H::AbstractHomotopy, x, t) = HomotopyNullCache()


# Homotopy API

"""
    evaluate!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` at `(x, t)` and store the result in `u`.
"""
evaluate!(u, H::AbstractHomotopy, args...) = error(MethodError(evaluate!, tuple(u, H, args...)))

"""
    evaluate(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` at `(x, t)`.
"""
evaluate(H::AbstractHomotopy, x, t) = evaluate(H, x, t, cache(H, x, t))
function evaluate(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)
    u = fill(zero(eltype(x)), size(H)[1])
    evaluate!(u, H, x, t, cache)
    u
end

"""
    dt!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` at `(x, t)` and store the result in `u`.
"""
function dt! end

"""
    dt(H::AbstractHomotopy, x::AbstractVector, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` at `(x, t)`.
"""
dt(H::AbstractHomotopy, x, t) = dt(H, x, t, cache(H, x, t))
function dt(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)
    u = fill(zero(eltype(x)), size(H)[1])
    dt!(u, H, x, t, cache)
    u
end

"""
    jacobian!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the Jacobian of the homotopy `H` at `(x, t)` and store the result in `u`.
"""
jacobian!(u, H::AbstractHomotopy, args...) = error(MethodError(jacobian!, tuple(u, H, args...)))

"""
    jacobian(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the Jacobian of the homotopy `H` at `(x, t)`.
"""
jacobian(H::AbstractHomotopy, x, t) = jacobian(H, x, t, cache(H, x, t))
function jacobian(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)
    U = fill(zero(eltype(x)), size(H))
    jacobian!(U, H, x, t, cache)
    U
end

# Optional
"""
    evaluate_and_jacobian!(u, U, F, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` and its Jacobian at `(x, t)` and store the results in `u` (evalution)
and `U` (Jacobian).
"""
function evaluate_and_jacobian!(u, U, H::AbstractHomotopy, x, t, c=cache(H, x, t))
    evaluate!(u, H, x, t, c)
    jacobian!(U, H, x, t, c)
    nothing
end

"""
    evaluate_and_jacobian(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` and its Jacobian at `(x, t)`.
"""
function evaluate_and_jacobian(H::AbstractHomotopy, x, t, c=cache(H, x, t))
    u = evaluate(H, x, t, c)
    U = jacobian(H, x, t, c)
    u, U
end

"""
    jacobian_and_dt!(U, u, H, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` and its derivative w.r.t. `t` at `(x, t)` and store the results in `u` (evalution)
and `v` (∂t).
"""
function jacobian_and_dt!(U, u, H::AbstractHomotopy, x, t, c=cache(H, x, t))
    jacobian!(U, H, x, t, c)
    dt!(u, H, x, t, c)
    nothing
end

"""
    jacobian_and_dt(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)

Evaluate the homotopy `H` and its derivative w.r.t. `t` at `(x, t)`.
"""
function jacobian_and_dt(H::AbstractHomotopy, x, t, c=cache(H, x, t))
    U = jacobian(H, x, t, c)
    u = dt(H, x, t, c)
    U, u
end

"""
    nvariables(H::AbstractHomotopy)

Returns the number of variables of the homotopy `H`.
"""
nvariables(H::AbstractHomotopy) = last(size(H))

Base.size(F::AbstractHomotopy, i::Integer) = size(F)[i]
Base.length(F::AbstractHomotopy) = size(F, 1)


"""
    precondition!(H::AbstractHomotopy, x, t, cache)

Prepare a homotopy for things like pathtracking starting at `x` and `t`.
This can modify `x` as well as `H` and anything in `cache`.
By default this is a no-op. If `H` wraps another homotopy this should call
`precondition!` on this as well.
"""
precondition!(H::AbstractHomotopy, x, t, cache) = nothing
"""
    update!(H::AbstractHomotopy, x, t, cache)

Update a homotopy for new values of `x` and `x`, i.e., update an affine patch.
This can modify `x` as well as `H` and anything in `cache`.
By default this is a no-op. If `H` wraps another homotopy this should call
`update!` on this as well.
"""
update!(H::AbstractHomotopy, x, t, cache) = nothing

"""
    basehomotopy(H::AbstractHomotopy)

Returns the 'proper' homotopy describing the problem. Any wrapper homotopy
recursively calls `wrappedhomotopy` with the wrapped homotopy as argument.
"""
basehomotopy(H::AbstractHomotopy) = H

"""
    Base.size(H::AbstractHomotopy)

Returns a tuple `(m, n)` indicating that `H` is a homotopy of `m` polynomials `m` in `n` variables.
"""
Base.size(::AbstractHomotopy) = error("Obligatory to define `Base.size($H)`")


include("homotopies/homotopy_witch_cache.jl")
include("homotopies/straight_line.jl")
include("homotopies/parameter_homotopy.jl")
include("homotopies/fixed_point.jl")
include("homotopies/patched_homotopy.jl")

function homotopy_interface_test(H::AbstractHomotopy, x=rand(Complex{Float64}, size(H, 2)) )
    m, n = size(H)
    t = rand()
    u = zeros(Complex{Float64}, m)
    U = zeros(Complex{Float64}, m, n)

    homotopy_cache = cache(H, x, t)

    evaluate!(u, H, x, t, homotopy_cache)
    @test evaluate(H, x, t, homotopy_cache) ≈ u atol=1e-14
    @test evaluate(H, x, t) ≈ u atol=1e-14

    dt!(u, H, x, t, homotopy_cache)
    @test dt(H, x, t, homotopy_cache) ≈ u atol=1e-14
    @test dt(H, x, t) ≈ u atol=1e-14

    jacobian!(U, H, x, t, homotopy_cache)
    @test jacobian(H, x, t, homotopy_cache) ≈ U atol=1e-14
    @test jacobian(H, x, t) ≈ U atol=1e-14

    evaluate_and_jacobian!(u, U, H, x, t, homotopy_cache)
    (v, V) = evaluate_and_jacobian(H, x, t, homotopy_cache)
    @test v ≈ u
    @test V ≈ U
    @test evaluate(H, x, t) ≈ u atol=1e-14
    @test jacobian(H, x, t) ≈ U atol=1e-14

    jacobian_and_dt!(U, u, H, x, t, homotopy_cache)
    (V, v) = jacobian_and_dt(H, x, t, homotopy_cache)
    @test V ≈ U atol=1e-14
    @test v ≈ u atol=1e-14
    @test jacobian(H, x, t) ≈ U atol=1e-14
    @test dt(H, x, t) ≈ u atol=1e-14
end
