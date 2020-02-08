export ParameterHomotopy, ModelKitHomotopy, total_degree_homotopy

abstract type AbstractHomotopy end

Base.size(H::AbstractHomotopy, i::Integer) = size(H)[i]

include("homotopies/model_kit_homotopy.jl")
include("homotopies/parameter_homotopy.jl")


##################
## Total Degree ##
##################
"""
    TotalDegreeStarts(degrees)

Given the `Vector`s `degrees` `TotalDegreeStarts` enumerates all solutions
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
struct TotalDegreeStarts{Iter}
    degrees::Vector{Int}
    iterator::Iter
end
function TotalDegreeStarts(degrees::Vector{Int})
    iterator = Iterators.product(map(d -> 0:(d-1), degrees)...)
    TotalDegreeStarts(degrees, iterator)
end
function Base.show(io::IO, iter::TotalDegreeStarts)
    print(io, "$(length(iter)) total degree start solutions for degrees $(iter.degrees)")
end

function Base.iterate(iter::TotalDegreeStarts)
    indices, state = iterate(iter.iterator)
    _value(iter, indices), state
end
function Base.iterate(iter::TotalDegreeStarts, state)
    it = iterate(iter.iterator, state)
    it === nothing && return nothing
    _value(iter, first(it)), last(it)
end

_value(iter::TotalDegreeStarts, indices) =
    map((k, d) -> cis(2π * k / d), indices, iter.degrees)


Base.length(iter::TotalDegreeStarts) = length(iter.iterator)
Base.eltype(iter::Type{<:TotalDegreeStarts}) = Vector{Complex{Float64}}

"""
    total_degree_homotopy(
        f::Vector{Expression},
        vars::Vector{Variable},
        params = Variable[];
        parameters = ComplexF64[],
        gamma = cis(2π * rand()))

Construct a total degree homotopy using the 'γ-trick'.
Returns a `ModelKitHomotopy` and the start solutions as an interator.
"""
function total_degree_homotopy(
    f::Vector{Expression},
    vars::Vector{Variable},
    params::Vector{Variable} = Variable[];
    parameters::AbstractVector = ComplexF64[],
    gamma = cis(2π * rand()),
)
    length(f) == length(vars) || throw(ArgumentError("Given system does not have the same number of polynomials as variables."))

    D = ModelKit.degrees(f, vars)

    t = ModelKit.unique_variable(:t, vars, params)
    γ = ModelKit.unique_variable(:γ, vars, params)
    h = γ .* t .* (vars .^ D .- 1) .+ (1 .- t) .* f

    H = ModelKitHomotopy(Homotopy(h, vars, t, [params; γ]), [parameters; gamma])
    S = TotalDegreeStarts(D)
    H, S
end
