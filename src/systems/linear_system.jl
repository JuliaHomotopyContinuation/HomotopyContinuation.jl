export LinearSystem

"""
The system `Ax-b`.
"""
struct LinearSystem <: ModelKit.AbstractSystem
    A::Matrix{ComplexF64}
    b::Vector{ComplexF64}
    ū::Vector{ComplexDF64}
    variables::Vector{Variable}
end

function LinearSystem(
    A::AbstractMatrix,
    b::AbstractVector;
    variables = first(@var x[1:size(A, 2)]),
)
    @assert size(A, 1) == length(b)
    LinearSystem(
        Matrix{ComplexF64}(A),
        Vector{ComplexF64}(b),
        Vector{ComplexDF64}(undef, size(A, 1)),
        variables,
    )
end
function LinearSystem(
    L::LinearSubspace,
    ::Coordinates{:Intrinsic};
    variables = first(@var x[1:codim(L)]),
)
    LinearSystem(intrinsic(L).A, -intrinsic(L).b₀; variables = variables)
end
function LinearSystem(
    L::LinearSubspace,
    ::Coordinates{:Extrinsic} = Extrinsic;
    variables = first(@var x[1:dim(L)]),
)
    LinearSystem(extrinsic(L).A, extrinsic(L).b; variables = variables)
end


(F::LinearSystem)(x, p = nothing) = F.A * x - F.b

ModelKit.variables(F::LinearSystem) = F.variables
ModelKit.parameters(::LinearSystem) = Variable[]

Base.size(F::LinearSystem) = size(F.A)

function ModelKit.evaluate!(u, F::LinearSystem, x::AbstractVector, p = nothing)
    LA.mul!(u, F.A, x)
    u .-= F.b
    u
end

function ModelKit.evaluate!(u, F::LinearSystem, x::Vector{ComplexDF64}, p = nothing)
    LA.mul!(F.ū, F.A, x)
    u .= F.ū .- F.b

    u
end

function ModelKit.evaluate_and_jacobian!(u, U, F::LinearSystem, x, p = nothing)
    evaluate!(u, F, x)
    U .= F.A
    nothing
end

function ModelKit.taylor!(
    u::AbstractVector,
    v::Val{K},
    F::LinearSystem,
    tx,
    p = nothing,
) where {K}
    m, n = size(F.A)
    u .= 0
    for j = 1:n, i = 1:m
        u[i] += F.A[i, j] * tx[j, K+1]
    end
    u
end

function ModelKit.taylor!(
    u::AbstractVector{<:TruncatedTaylorSeries},
    v::Val{K},
    F::LinearSystem,
    tx,
    p = nothing,
) where {K}
    m, n = size(F.A)
    u .= 0
    for j = 1:n, i = 1:m
        u[i] = u[i] + F.A[i, j] * tx[j]
    end
    for i = 1:m
        u[i, 1] -= F.b[i]
    end

    u
end

