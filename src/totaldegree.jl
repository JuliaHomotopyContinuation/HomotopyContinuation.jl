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
