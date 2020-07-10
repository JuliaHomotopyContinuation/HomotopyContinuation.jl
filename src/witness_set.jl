export WitnessSet, witness_set, linear_subspace, system, dim, codim, trace_test

"""
    WitnessSet(F, L, S)

Store solutions `S` of the polynomial system `F(x) = L(x) = 0` into a witness set.
"""
struct WitnessSet{
    S<:AbstractSystem,
    Sub<:AbstractSubspace,
    R<:Union{Vector{ComplexF64},PathResult},
}
    F::S
    L::Sub
    # only non-singular results or solutions
    R::Vector{R}
    projective::Bool
end

function WitnessSet(
    F::AbstractSystem,
    L::LinearSubspace,
    R;
    projective::Bool = is_linear(A) && is_homogeneous(System(F)),
)
    WitnessSet(F, L, R, projective)
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
function witness_set(
    F::AbstractSystem;
    target_parameters = nothing,
    dim = nothing,
    codim = nothing,
    options...,
)
    f = System(F)
    projective = is_homogeneous(f)
    if isnothing(dim) && isnothing(codim)
        dim = corank(F) - projective
    elseif !isnothing(codim) && projective
        codim += 1
    end
    L = rand_subspace(size(F, 2); dim = codim, codim = dim, affine = !projective)
    witness_set(
        F,
        L;
        target_parameters = target_parameters,
        projective = projective,
        options...,
    )
end

function witness_set(F::AbstractSystem, L::LinearSubspace; projective = nothing, options...)
    res = solve(F; target_subspace = L, options...)
    if isnothing(projective)
        WitnessSet(F, L, results(res; only_nonsingular = true))
    else
        WitnessSet(F, L, results(res; only_nonsingular = true); projective = projective)
    end
end

corank(F::AbstractSystem; kwargs...) = nvariables(F) - LA.rank(F)
function LinearAlgebra.rank(F::AbstractSystem; target_parameters = nothing)
    m, n = size(F)
    u = zeros(ComplexF64, m)
    U = zeros(ComplexF64, m, n)
    x = randn(ComplexF64, n)
    evaluate_and_jacobian!(u, U, F, x, target_parameters)
    LA.rank(U)
end


### Move witness sets around
function witness_set(W::WitnessSet, L::LinearSubspace; options...)
    res = solve(W.F, W.R; start_subspace = W.L, target_subspace = L, options...)
    WitnessSet(W.F, L, results(res; only_nonsingular = true), is_linear(L) && W.projective)
end

"""
    trace_test(W::WitnessSet; options...)

Performs a trace test [^LRS18] to verify whether the given witness set `W` is complete.
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

[^LRS18]: Leykin, Anton, Jose Israel Rodriguez, and Frank Sottile. "Trace test." Arnold Mathematical Journal 4.1 (2018): 113-125.
APA

"""
function trace_test(W₀::WitnessSet; options...)
    L₀ = linear_subspace(W₀)
    F = system(W₀)

    # if we are in the projective setting, we need to make sure that
    # all solutions are on the same affine chart
    # Therefore make the affine chart now
    if W₀.projective
        F = on_affine_chart(F)
        s₀ = sum(s -> set_solution!(s, F, s), solutions(W₀))
    else
        s₀ = sum(solutions(W₀))
    end

    v = randn(ComplexF64, codim(L₀))
    L₁ = translate(L₀, v)
    L₋₁ = translate(L₀, -v)

    R₁ = solve(F, solutions(W₀); start_subspace = L₀, target_subspace = L₁, options...)
    nsolutions(R₁) == degree(W₀) || return nothing

    R₋₁ = solve(F, solutions(W₀); start_subspace = L₀, target_subspace = L₋₁, options...)
    nsolutions(R₋₁) == degree(W₀) || return nothing

    s₁ = sum(solutions(R₁))
    s₋₁ = sum(solutions(R₋₁))
    Δs₁ = (s₁ - s₀)
    Δs₋₁ = (s₀ - s₋₁)
    trace = InfNorm()(Δs₁, Δs₋₁) / max(InfNorm()(Δs₁), InfNorm()(Δs₋₁))

    trace
end
