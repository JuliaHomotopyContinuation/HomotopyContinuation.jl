export RandomizedSystem, square_up

"""
    RandomizedSystem(F::Union{System,AbstractSystem}, k::Integer;
            identity_block = true,
            compile = $(COMPILE_DEFAULT[]))

Given a ``n × N`` system `F`` this constructs the system
``\\mathfrak{R}(F; k)(x) = A⋅F(x)`` where ``A`` is a ``k × N`` matrix.
If `identity_block = true` then the first `k` columns of ``A`` form an identity matrix.
The other entries are random complex numbers. See Chapter 13.5 in [^SW05] for more details.

    RandomizedSystem(F::Union{System,AbstractSystem}, A::Matrix{ComplexF64};
        identity_block = true,
        compile = $(COMPILE_DEFAULT[]))

Explicitly provide the used randomization matrix. If `identity = true` then `A` has to be `k × (n - k)` matrix.
Otherwise `A` has to be a `k × n` matrix..

[^SW05]: Sommese, A. J., & Wampler, C. W. (2005). The Numerical Solution of Systems of Polynomials Arising in Engineering and Science. World Scientific.
"""
struct RandomizedSystem{S<:AbstractSystem} <: AbstractSystem
    system::CompositionSystem{S,LinearSystem}
end

function RandomizedSystem(
    F::Union{AbstractSystem,System},
    k::Integer;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
)
    n, = size(F)
    # Create a random matrix unit upper triangular
    A = [LA.UnitUpperTriangular(rand_unitary_matrix(k)) randn(ComplexF64, k, n - k)]
    RandomizedSystem(F, A; compile = compile)
end

RandomizedSystem(
    F::System,
    A::Matrix{ComplexF64};
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    kwargs...,
) = RandomizedSystem(fixed(F; compile = compile), A; kwargs...)

function RandomizedSystem(
    F::AbstractSystem,
    A::Matrix{ComplexF64};
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
)
    L = LinearSystem(A, zeros(ComplexF64, size(A, 1)); variables = variables(F))
    system = CompositionSystem(L, F)

    RandomizedSystem(system)
end

function Base.show(io::IO, mime::MIME"text/plain", F::RandomizedSystem)
    show(io, mime, F.system)
end

Base.size(F::RandomizedSystem) = size(F.system)
ModelKit.variables(F::RandomizedSystem) = variables(F.system)
ModelKit.parameters(F::RandomizedSystem) = parameters(F.system)
ModelKit.variable_groups(F::RandomizedSystem) = variable_groups(F.system)

"""
    square_up(F::Union{System, AbstractSystem}; identity_block = true, compile = compile = $(COMPILE_DEFAULT[]))

Creates the [`RandomizedSystem`](@ref) ``\\mathfrak{R}(F(x); N)`` where ``N`` is the number
of variables of `F`.
"""
square_up(F::Union{AbstractSystem,System}; kwargs...) =
    RandomizedSystem(F, last(size(F)); kwargs...)


(F::RandomizedSystem)(x, p = nothing) = F.system(x, p)
ModelKit.evaluate!(u, F::RandomizedSystem, x::AbstractVector, p = nothing) =
    evaluate!(u, F.system, x, p)
ModelKit.evaluate_and_jacobian!(u, U, F::RandomizedSystem, x, p = nothing) =
    evaluate_and_jacobian!(u, U, F.system, x, p)
ModelKit.taylor!(u, v::Val, F::RandomizedSystem, tx, p = nothing) =
    taylor!(u, v, F.system, tx, p)
