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
    system::S
    A::Matrix{ComplexF64}
    identity_block::Bool
    u::Vector{ComplexF64}
    ū::Vector{ComplexDF64}
    U::Matrix{ComplexF64}
    taylor_u::TaylorVector{4,ComplexF64}
end

function RandomizedSystem(
    F::Union{AbstractSystem,System},
    k::Integer;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
    identity_block = true,
)
    n, N = size(F)
    if identity_block
        A = randn(ComplexF64, k, n - k)
    else
        # Create a random matrix unit upper triangular
        A = [LA.UnitUpperTriangular(rand_unitary_matrix(k)) randn(ComplexF64, k, n - k)]
    end
    RandomizedSystem(F, A; identity_block = identity_block, compile = compile)
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
    identity_block = true,
)
    n, N = size(F)
    is_compatible = identity_block ? size(A, 1) + size(A, 2) == n : size(A, 2) == n
    if !is_compatible
        throw(ArgumentError("Randomization matrix has incompatible size."))
    end

    u = zeros(ComplexF64, n)
    ū = zeros(ComplexDF64, n)
    U = zeros(ComplexF64, n, N)
    taylor_u = TaylorVector{4}(ComplexF64, n)
    RandomizedSystem(F, A, identity_block, u, ū, U, taylor_u)
end

function Base.show(io::IO, mime::MIME"text/plain", F::RandomizedSystem)
    println(io, typeof(F), ":")
    println(io, "A:")
    show(io, mime, F.A)
    println(io, "\n\nF:")
    show(io, mime, F.system)
end

Base.size(F::RandomizedSystem) = (size(F.A, 1), last(size(F.system)))
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


@inline function randomize!(u, A, v::AbstractVector; identity_block::Bool)
    n, m = size(A, 1), length(v)

    if identity_block
        for i = 1:n
            u[i] = v[i]
        end
        for j = 1:(m-n)
            vnj = v[n+j]
            for i = 1:n
                u[i] += A[i, j] * vnj
            end
        end
    else
        for i = 1:n
            u[i] = 0
        end
        for j = 1:m, i = 1:n
            u[i] += A[i, j] * v[j]
        end
    end

    u
end

@inline function randomize!(
    u::Union{AbstractMatrix,TaylorVector},
    A,
    v::Union{AbstractMatrix,TaylorVector},
    ncols::Int = size(A, 1);
    identity_block::Bool,
)
    n, m = size(A, 1), size(v, 1)

    if identity_block
        for j = 1:ncols, i = 1:n
            u[i, j] = v[i, j]
        end
        for j = 1:ncols, k = 1:(m-n), i = 1:n
            u[i, j] += A[i, k] * v[n+k, j]
        end
    else
        for j = 1:ncols, i = 1:n
            u[i, j] = 0
        end
        for j = 1:ncols, k = 1:m, i = 1:n
            u[i, j] += A[i, k] * v[k, j]
        end
    end

    u
end

function (F::RandomizedSystem)(x, p = nothing)
    if F.identity_block
        [LA.I F.A] * F.system(x, p)
    else
        F.A * F.system(x, p)
    end
end

function ModelKit.evaluate!(u, F::RandomizedSystem, x::Vector{ComplexDF64}, p = nothing)
    evaluate!(F.ū, F.system, x, p)
    randomize!(u, F.A, F.ū; identity_block = F.identity_block)
    u
end
function ModelKit.evaluate!(u, F::RandomizedSystem, x, p = nothing)
    evaluate!(F.u, F.system, x, p)
    randomize!(u, F.A, F.u; identity_block = F.identity_block)
    u
end

function ModelKit.evaluate_and_jacobian!(u, U, F::RandomizedSystem, x, p = nothing)
    evaluate_and_jacobian!(F.u, F.U, F.system, x, p)
    randomize!(u, F.A, F.u; identity_block = F.identity_block)
    randomize!(U, F.A, F.U; identity_block = F.identity_block)
end

function ModelKit.taylor!(
    u::AbstractVector,
    v::Val{N},
    F::RandomizedSystem,
    tx,
    p = nothing,
) where {N}
    taylor!(F.u, v, F.system, tx, p)
    randomize!(u, F.A, F.u; identity_block = F.identity_block)
end
function ModelKit.taylor!(
    u::TaylorVector,
    v::Val{N},
    F::RandomizedSystem,
    tx,
    p = nothing,
) where {N}
    taylor!(F.taylor_u, v, F.system, tx, p)
    randomize!(u, F.A, F.taylor_u, N + 1; identity_block = F.identity_block)
end
