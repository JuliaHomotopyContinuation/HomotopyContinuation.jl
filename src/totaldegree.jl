export bezout_number

"""
    totaldegree_solutions(F::Vector{<:MP.AbstractPolynomialLike})

Returns an iterator of the solutions of the total degree startsystem of `F`.
"""
function totaldegree_solutions(F::Vector{AP}) where {AP<:MP.AbstractPolynomialLike}
    TotalDegreeSolutionIterator(MP.maxdegree.(F))
end
totaldegree_solutions(degrees::Vector{Int}) = TotalDegreeSolutionIterator(degrees)
totaldegree_solutions(F::TotalDegreeSystem, _) = TotalDegreeSolutionIterator(F.degrees)
function totaldegree_solutions(F::MultiHomTotalDegreeSystem, vargroups)
    MultiBezoutSolutionsIterator(F.D, F.C, vargroups)
end


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
# Multihomogeneous
##################

struct MultiBezoutIndicesIterator{N}
    D::Matrix{Int} # Degree matrix
    k::NTuple{N, Int} # (projective) dimensions

    function MultiBezoutIndicesIterator(D::Matrix{Int}, k::NTuple{M, Int}) where {M}
        @assert size(D, 1) == M
        new{M}(D, k)
    end
end
MultiBezoutIndicesIterator(D::Matrix, k) = MultiBezoutIndicesIterator(D, tuple(k...))
MultiBezoutIndicesIterator(D::Matrix, groups::VariableGroups) = MultiBezoutIndicesIterator(D, projective_dims(groups))

Base.IteratorSize(::Type{MultiBezoutIndicesIterator}) = Base.SizeUnknown()
Base.eltype(::Type{MultiBezoutIndicesIterator}) = Tuple{Int, Vector{Int}}
function Base.iterate(iter::MultiBezoutIndicesIterator, state=[i for i=1:length(iter.k) for j=1:iter.k[i]])
    d = 0
    m, n = size(iter.D)
    p = state
    while d == 0
        state[1] > n && return nothing
        p, state = nextpermutation(1:m, state)
        let p=p, D=iter.D
            d = prod(i -> D[p[i], i], 1:n)::Int
        end
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
    bezout_number(F::MPPolys; variable_groups=[variables(F)], homvars=nothing, parameters=nothing)
    bezout_number(multidegrees, groups::VariableGroups)

Compute the multi-homogeneous bezout number associated to the given system and variable groups.
"""
bezout_number(D::Matrix, groups::VariableGroups) = bezout_number(D, projective_dims(groups))
bezout_number(D::Matrix, k) = sum(first, MultiBezoutIndicesIterator(D, k))
function bezout_number(F::MPPolys; parameters=nothing, variable_groups=nothing, homvars=nothing)
    if variable_groups === nothing
        return prod(maxdegrees(F; parameters=parameters))
    end

    vars = variables(F; parameters=parameters)
    hominfo = HomogenizationInformation(variable_groups=variable_groups, homvars=homvars)
    D = multidegrees(F, variable_groups)
    if is_homogeneous(F, hominfo; parameters=parameters)
        bezout_number(D, length.(variable_groups) .- 1)
    else
        bezout_number(D, length.(variable_groups))
    end
end

"""
    multi_bezout_coefficients(D, projective_dims)

Compute coefficients for a multi-bezout start system.
`projective_dims` is the projective dimensions of the projective spaces.
"""
function multi_bezout_coefficients(D::Matrix, k::NTuple{M, Int}) where M
    @assert size(D, 1) == M
    @assert size(D, 2) == sum(k)

    n = size(D, 2)
    C = zeros(n+M, n)
    for i in 1:n
        j = 1
        for kⱼ in k
            for l in 1:kⱼ
                C[j, i] = begin
                    if (i == l && i ≤ kⱼ)
                        1.0
                    elseif i ≠ l && i ≤ kⱼ
                        0.0
                    else
                        randn()
                    end
                end
                j += 1
            end
            C[j, i] = 1.0
            j += 1
        end
    end
    C
end


struct MultiBezoutSolutionsIterator{N}
    indices::MultiBezoutIndicesIterator{N}
    coeffs::Matrix{Float64}
    # precomputed stuff and cache
    roots_of_unity::Matrix{ComplexF64} # Lower Triangular matrix storing the precomputed roots
    A::NTuple{N, Matrix{Float64}}
    b::NTuple{N, Vector{ComplexF64}}
    ranges::NTuple{N, UnitRange{Int}}
end

function MultiBezoutSolutionsIterator(indices::MultiBezoutIndicesIterator, coeffs::Matrix)
    d_max = maximum(indices.D)
    roots_of_unity = zeros(ComplexF64, d_max, d_max)
    for i in 1:d_max, j in 0:i-1
        roots_of_unity[i, j+1] = cis(2π * j / i)
    end
    A = map(kⱼ -> zeros(kⱼ, kⱼ), indices.k)
    b = map(kⱼ -> zeros(ComplexF64, kⱼ), indices.k)
    r_stop = Ref(1)
    ranges = map(indices.k) do kⱼ
        range = r_stop[]:(kⱼ + r_stop[] - 1)
        r_stop[] += kⱼ + 1
        range
    end
    MultiBezoutSolutionsIterator(indices, coeffs, roots_of_unity, A, b, ranges)
end
function MultiBezoutSolutionsIterator(D::Matrix, C::Matrix, vargroups::VariableGroups)
    MultiBezoutSolutionsIterator(MultiBezoutIndicesIterator(D, vargroups), C)
end

function Base.show(io::IO, iter::MultiBezoutSolutionsIterator)
    print(io, "Solutions iterator for a multi-homogeneous start system")
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::MultiBezoutSolutionsIterator) = x

Base.length(iter::MultiBezoutSolutionsIterator) = bezout_number(iter.indices.D, iter.indices.k)
Base.IteratorSize(::Type{<:MultiBezoutSolutionsIterator}) = Base.SizeUnknown()
Base.eltype(::Type{MultiBezoutSolutionsIterator{N}}) where N = ProjectiveVectors.PVector{ComplexF64, N}

function Base.iterate(iter::MultiBezoutSolutionsIterator)
    (_, perm), indices_state = iterate(iter.indices)
    dᵢ = ntuple(i -> iter.indices.D[perm[i],i], length(perm))
    d_iter = Iterators.product(map(dᵢ -> 1:dᵢ, dᵢ)...)
    q, d_state = iterate(d_iter)

    x = compute_solution(iter, perm, q, dᵢ)
    x, (perm, indices_state, dᵢ, d_iter, d_state)
end

function Base.iterate(iter::MultiBezoutSolutionsIterator, state)
    perm, indices_state, dᵢ, d_iter, d_state = state

    next_d_state = iterate(d_iter, d_state)
    if next_d_state !== nothing
        q, d_state = next_d_state
        x = compute_solution(iter, perm, q, dᵢ)
        return x, (perm, indices_state, dᵢ, d_iter, d_state)
    end

    # We exhausted our current batch. move indices forward
    indices_newstate = iterate(iter.indices, indices_state)
    # check whether we are completely done
    indices_newstate === nothing && return nothing

    (_, perm), indices_state = indices_newstate
    # Create new d_iter
    new_dᵢ::typeof(dᵢ) = let D=iter.indices.D, perm=perm
        ntuple(i -> D[perm[i],i], length(perm))
    end
    new_d_iter::typeof(d_iter) = Iterators.product(map(dᵢⱼ -> 1:dᵢⱼ, new_dᵢ)...)
    q, d_state = iterate(new_d_iter)
    x = compute_solution(iter, perm, q, new_dᵢ)
    return x, (perm, indices_state, new_dᵢ, new_d_iter, d_state)
end


function compute_solution(iter::MultiBezoutSolutionsIterator, perm, q, dᵢ)
    m, n = size(iter.indices.D)
    # for each variable group we have to solve a linear system
    t = 1
    data = zeros(ComplexF64, n + m)
    for j=1:m
        kⱼ = iter.indices.k[j]
        Aⱼ = iter.A[j]
        bⱼ = iter.b[j]
        s = 1
        for i=1:n
            perm[i] == j || continue
            range = iter.ranges[j]
            for (t, l)  in enumerate(range)
                Aⱼ[s, t] = iter.coeffs[l, i]
            end
            bⱼ[s] = iter.roots_of_unity[dᵢ[i], q[i]]
            s += 1
        end

        LinearAlgebra.ldiv!(LinearAlgebra.lu!(Aⱼ), bⱼ)

        data[t:(t+kⱼ-1)] .= bⱼ
        data[t+kⱼ] = 1
        t += kⱼ + 1
    end

    x = ProjectiveVectors.PVector(data, iter.indices.k)
end
