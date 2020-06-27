export WitnessSet,
    witness_set, affine_subspace, linear_subspace, system, degree, dim, codim, trace_test

"""
    WitnessSet(F, L, S)

Store solutions `S` of the polynomial system `F(x) = L(x) = 0` into a witness set.
"""
struct WitnessSet{
    Sub<:AbstractSubspace,
    S<:AbstractSystem,
    R<:Union{Vector{ComplexF64},PathResult},
}
    F::S
    L::Sub
    # only non-singular results or solutions
    R::Vector{R}
end

"""
    system(W::WitnessSet)

Get the system stored in `W`.
"""
system(W::WitnessSet) = W.F

"""
    linear_subspace(W::WitnessSet)

Get the linear subspace stored in `W`.
"""
linear_subspace(W::WitnessSet) = W.L

"""
    affine_subspace(W::WitnessSet)

Get the affine linear subspace stored in `W`.
"""
affine_subspace(W::WitnessSet{<:LinearSubspace}) = W.L

"""
    solutions(W::WitnessSet)

Get the solutions stored in `W`.
"""
solutions(W::WitnessSet{A,B,PathResult}) where {A,B} = solutions(W.R)
solutions(W::WitnessSet{A,B,Vector{ComplexF64}}) where {A,B} = W.R

"""
    results(W::WitnessSet)

Get the results stored in `W`.
"""
results(W::WitnessSet{<:Any,<:Any,PathResult}) = W.R

"""
    degree(W::WitnessSet)

Returns the degree of the witness set `W`. This equals the number of solutions stored.
"""
ModelKit.degree(W::WitnessSet) = length(W.R)

"""
    dim(W::WitnessSet)

The dimension of the algebraic set encoded by the witness set.
"""
dim(W::WitnessSet) = codim(W.L)

"""
    codim(W::WitnessSet)

The dimension of the algebraic set encoded by the witness set.
"""
codim(W::WitnessSet) = dim(W.L)

function Base.show(io::IO, W::WitnessSet)
    print(io, "Witness set for dimension $(dim(W)) of degree $(degree(W))")
end

### Construct witness sets
"""
    witness_set(F; codim = nvariables(F) - length(F), dim = nothing, options...)

Compute a [`WitnessSet`](@ref) for `F` in the given dimension (resp. codimension)
by sampling a random (affine) linear subspace.
After constructing the system this calls [`solve`](@ref) with the provided `options`.

    witness_set(F, L; options...)

Compute [`WitnessSet`](@ref) for `F` and the (affine) linear subspace `L`.

    witness_set(W::WitnessSet, L; options...)

Compute a new [`WitnessSet`](@ref) with the (affine) linear subspace `L` by moving
the linear subspace stored in `W` to `L`.

### Example
```julia-repl
julia> @var x y;
julia> F = System([x^2 + y^2 - 5], [x, y])
System of length 1
 2 variables: x, y

 -5 + x^2 + y^2

julia> W = witness_set(F)
Witness set for dimension 1 of degree 2
```
"""
witness_set(F::System, args...; compile = COMPILE_DEFAULT[], kwargs...) =
    witness_set(fixed(F; compile = compile), args...; kwargs...)
function witness_set(F::AbstractSystem; dim = nothing, codim = nothing, options...)
    if isnothing(dim) && isnothing(codim)
        codim = size(F, 2) - size(F, 1)
    end
    f = System(F)
    if !is_homogeneous(f)
        L = rand_affine_subspace(size(F, 2); dim = dim, codim = codim)
    else
        error("TODO")
    end
    witness_set(F, L; options...)
end

function witness_set(F::AbstractSystem, L::LinearSubspace; options...)
    F_L = slice(F, L)
    res = solve(F_L; options...)
    R = results(res; only_nonsingular = true)
    WitnessSet(F, L, R)
end

### Move witness sets around
function witness_set(W::WitnessSet, L::LinearSubspace; options...)
    H = LinearSubspaceHomotopy(W.F, W.L, L)
    res = solve(H, W.R; options...)
    WitnessSet(W.F, L, results(res; only_nonsingular = true))
end

"""
    trace_test(W::WitnessSet; options...)

Performs a trace test to verify whether the given witness set `W` is complete.
Returns the trace of the witness set which should be theoretically be 0 if `W` is complete.
Due to floating point arithmetic this is not the case, thus is has to be manually checked
that the trace is sufficiently small.
Returns `nothing` if the trace test failed due to path tracking failures.
The `options` are the same as for calls to [`witness_set`](@ref).

```julia-repl
julia> @var x y;
julia> F = System([x^2 + y^2 - 5], [x, y])
System of length 1
 2 variables: x, y

 -5 + x^2 + y^2

julia> W = witness_set(F)
Witness set for dimension 1 of degree 2

julia> trace = trace_test(W)
9.981960497718987e-16
```
"""
function trace_test(W₀::WitnessSet; options...)
    L₀ = affine_subspace(W₀)

    v = randn(ComplexF64, codim(L₀))
    L₁ = translate(L₀, v)
    L₋₁ = translate(L₀, -v)

    W₁ = witness_set(W₀, L₁; options...)
    degree(W₁) == degree(W₀) || return nothing

    W₋₁ = witness_set(W₀, L₋₁; options...)
    degree(W₋₁) == degree(W₀) || return nothing

    s₀ = sum(solutions(W₀))
    s₁ = sum(solutions(W₁))
    s₋₁ = sum(solutions(W₋₁))
    Δs₁ = (s₁ - s₀)
    Δs₋₁ = (s₀ - s₋₁)
    trace = InfNorm()(Δs₁, Δs₋₁) / max(norm(Δs₁), norm(Δs₋₁))

    trace
end
