"""
    totaldegree_solutions(F::Vector{<:MP.AbstractPolynomialLike})

Returns an iterator of the solutions of the total degree startsystem of `F`.
"""
function totaldegree_solutions(F::Vector{AP}) where {AP<:MP.AbstractPolynomialLike}
    totaldegree_solutions(MP.maxdegree.(F))
end
totaldegree_solutions(degrees::Vector{Int}) = TotalDegreeSolutionIterator(degrees)


"""
    TotalDegreeSolutionIterator(degrees)

Given the `Vector`s `degrees` `TotalDegreeSolutionIterator` enumerates all solutions
of the system
```math
\\begin{align*}
    z_1^{d_1} - 1 &= 0 \\\\
    z_1^{d_2} - 1 &= 0 \\\\
    &\\vdots \\\\
    z_n^{d_n} - 1 &= 0 \\\\
\\end{align*}
```
where ``d_i`` is `degrees[i]`.
"""
struct TotalDegreeSolutionIterator{Iter}
    degrees::Vector{Int}
    iterator::Iter
end
function TotalDegreeSolutionIterator(degrees::Vector{Int})
    iterator = Base.Iterators.product(map(d -> 0:d-1, degrees)...)
    TotalDegreeSolutionIterator(degrees, iterator)
end
function Base.show(io::IO, iter::TotalDegreeSolutionIterator)
    print(io, "TotalDegreeSolutionIterator(degrees=$(iter.degrees))")
end

function Base.iterate(iter::TotalDegreeSolutionIterator)
    indices, state = iterate(iter.iterator)
    _value(iter, indices), state
end
function Base.iterate(iter::TotalDegreeSolutionIterator, state)
    it = iterate(iter.iterator, state)
    it === nothing && return nothing
    _value(iter, it[1]), it[2]
end

function _value(iter::TotalDegreeSolutionIterator, indices)
    value = Vector{Complex{Float64}}()
    for (i, k) in enumerate(indices)
        push!(value, cis(2π * k / iter.degrees[i]))
    end
    value
end

Base.length(iter::TotalDegreeSolutionIterator) = length(iter.iterator)
Base.eltype(iter::Type{<:TotalDegreeSolutionIterator}) = Vector{Complex{Float64}}


##################
# Multihomogenous
##################

struct MultiBezoutIndicesIterator
    D::Matrix{Int} # Degree matrix
    k::Vector{Int} # size of (affine) variable groups
end
MultiBezoutIndicesIterator(D::Matrix, k) = MultiBezoutIndicesIterator(D, collect(k))

Base.IteratorSize(::Type{MultiBezoutIndicesIterator}) = Base.SizeUnknown()
Base.eltype(::Type{MultiBezoutIndicesIterator}) = Tuple{Int, Vector{Int}}
function Base.iterate(iter::MultiBezoutIndicesIterator, state=[i for i=1:length(iter.k) for j=1:iter.k[i]])
    d = 0
    n, m = size(iter.D)
    p = state
    while d == 0
        state[1] > n && return nothing
        p, state = nextpermutation(1:m, state)
        d = prod(i -> iter.D[i, p[i]], 1:n)
    end

    (d, p), state
end

# Adopted from Combinatorics.jl:
# https://github.com/JuliaMath/Combinatorics.jl/blob/8d9571402319799b29da2005a65b627e8771c1e4/src/permutations.jl#L47
function nextpermutation(m, state)
    n = length(state)
    perm = [m[state[i]] for i in 1:n]
    s = copy(state)
    i = n - 1
    while i>=1 && s[i] >= s[i+1]
        i -= 1
    end
    if i > 0
        j = n
        while j > i && s[i] >= s[j]
            j -= 1
        end
        s[i], s[j] = s[j], s[i]
        reverse!(s, i+1)
    else
        s[1] = n+1
    end
    return (perm, s)
end

"""
    bezout_number(multidegrees, groups::VariableGroups)

Compute the multi-homogenous bezout number associated to the given multidegrees and variable groups.
"""
bezout_number(D::Matrix, groups::VariableGroups) = bezout_number(D, projective_dims(groups))
bezout_number(D::Matrix, k) = sum(first, MultiBezoutIndicesIterator(D, k))


"""
    totaldegree_polysystem(multidegrees::Matrix, variables, variable_groups::VariableGroups)

The multi-homogenous totaldegree start system described in [Wampler, 93]. Returns a tuple, the system
and the coefficients``c_{i,j,l}`` described in [Wampler, 93] as a `Matrix{Vector{Float64}}`.

[Wampler, 93]: An efficient start system for multi-homogeneous polynomial continuation (https://link.springer.com/article/10.1007/BF01385710).
"""
function totaldegree_polysystem(multidegrees::Matrix, variables, variable_groups::VariableGroups)
    n, m = size(multidegrees)
    Z = groups(variable_groups, variables)
    c(i, j) = begin
        kⱼ = length(Z[j])
        map(1:kⱼ) do
            if (i == l && i ≤ kⱼ) || l == 0
                1
            elseif i ≠ l && i ≤ kⱼ
                0
            else
                randn()
            end
        end
    end
    C = [c(i,j) for i ∈ 1:n, j ∈ 1:m]
    G = map(1:n) do i
        prod(1:m) do j
            dᵢⱼ = multidegrees[i,j]
            dᵢⱼ == 0 && return 0
            kⱼ = length(Z[j])
            bᵢⱼ = sum(1:kⱼ) do l
                C[i,j][l] * Z[j][l]
            end
        end
    end
    G, C
end
