export total_degree, total_degree_start_solutions

struct TotalDegreeStartSolutionsIterator{Iter}
    degrees::Vector{Int}
    iterator::Iter
end
function TotalDegreeStartSolutionsIterator(degrees)
    iterator = Iterators.product(map(d -> 0:(d-1), degrees)...)
    TotalDegreeStartSolutionsIterator(degrees, iterator)
end
function Base.show(io::IO, iter::TotalDegreeStartSolutionsIterator)
    print(io, "$(length(iter)) total degree start solutions for degrees $(iter.degrees)")
end

function Base.iterate(iter::TotalDegreeStartSolutionsIterator)
    indices, state = iterate(iter.iterator)
    _value(iter, indices), state
end
function Base.iterate(iter::TotalDegreeStartSolutionsIterator, state)
    it = iterate(iter.iterator, state)
    it === nothing && return nothing
    _value(iter, first(it)), last(it)
end

_value(iter::TotalDegreeStartSolutionsIterator, indices) =
    map((k, d) -> cis(2π * k / d), indices, iter.degrees)


Base.length(iter::TotalDegreeStartSolutionsIterator) = length(iter.iterator)
Base.eltype(iter::Type{<:TotalDegreeStartSolutionsIterator}) = Vector{Complex{Float64}}

"""
    total_degree_start_solutions(degrees)

Returns an iterator of the start solutions of the system
```math
\\begin{array}{rl}
    z_1^{d_1} - 1 &= 0 \\\\
    z_2^{d_2} - 1 &= 0 \\\\
    &\\vdots \\\\
    z_n^{d_n} - 1 &= 0 \\\\
\\end{array}
```
where ``d_i`` is `degrees[i]`.

## Example
```julia-repl
julia> iter = total_degree_start_solutions([2, 2])
4 total degree start solutions for degrees [2, 2]

julia> collect(iter)
4-element Array{Array{Complex{Float64},1},1}:
 [1.0 + 0.0im, 1.0 + 0.0im]
 [-1.0 + 1.2246467991473532e-16im, 1.0 + 0.0im]
 [1.0 + 0.0im, -1.0 + 1.2246467991473532e-16im]
 [-1.0 + 1.2246467991473532e-16im, -1.0 + 1.2246467991473532e-16im]
```
"""
total_degree_start_solutions(degrees) = TotalDegreeStartSolutionsIterator(degrees)

"""
    total_degree(F::System; parameters = nothing, gamma = cis(2π * rand()))
    total_degree(F::AbstractSystem; gamma = cis(2π * rand()))

Construct a total degree homotopy ``H`` using the 'γ-trick' such that ``H(x,0) = F(x)``
and ``H(x,1)`` is the generated start system.
Returns a `StraightLineHomotopy` and the start solutions as an interator.
"""
total_degree(F::System; target_parameters = nothing, kwargs...) =
    total_degree(ModelKitSystem(F, target_parameters); kwargs...)

function total_degree(
    F::AbstractSystem;
    γ = cis(2π * rand()),
    gamma = γ,
    scaling::Union{Nothing,AbstractVector} = nothing,
    iterator::Bool = true,
)
    m, n = size(F)
    m == n ||
    throw(ArgumentError("Given system does not have the same number of polynomials as variables."))

    @unique_var x[1:n] t s[1:n]
    dicts = ModelKit.to_dict.(expand.(F(x)), Ref(x))
    D = maximum.(sum, keys.(dicts))

    if scaling === nothing
        scaling = map(C -> maximum(float ∘ abs ∘ to_number, C), values.(dicts))
    end
    G = ModelKitSystem(System(s .* (x .^ D .- 1), x, s), scaling)
    StraightLineHomotopy(G, F; γ = gamma), total_degree_start_solutions(D)
end
