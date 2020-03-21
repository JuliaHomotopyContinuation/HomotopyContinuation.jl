export ParameterHomotopy, ModelKitHomotopy, total_degree_homotopy

abstract type AbstractHomotopy end

Base.size(H::AbstractHomotopy, i::Integer) = size(H)[i]
on_chart!(x::Vector, H::AbstractHomotopy) = x

include("homotopies/differentiation.jl")

include("homotopies/model_kit_homotopy.jl")
include("homotopies/parameter_homotopy.jl")
include("homotopies/affine_chart_homotopy.jl")
include("homotopies/straight_line_homotopy.jl")
include("homotopies/coefficient_homotopy.jl")
include("homotopies/polyhedral_homotopy.jl")


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
    total_degree_homotopy(F::System; parameters = ComplexF64[], gamma = cis(2π * rand()))

Construct a total degree homotopy ``H`` using the 'γ-trick' such that ``H(x,0) = F(x)``
and ``H(x,1)`` is the generated start system.
Returns a `ModelKitHomotopy` and the start solutions as an interator.
"""
total_degree_homotopy(f::AbstractVector{Expression}, vars::AbstractVector{Variable}, params::AbstractVector{Variable} = Variable[]; kwargs...) =
    total_degree_homotopy(System(f, vars, params); kwargs...)

function total_degree_homotopy(
    F::System;
    target_parameters::AbstractVector = ComplexF64[],
    gamma = cis(2π * rand()),
    scaling = nothing,
)
    n = length(F)
    n == ModelKit.nvariables(F) ||
    throw(ArgumentError("Given system does not have the same number of polynomials as variables."))
    ModelKit.nparameters(F) == length(target_parameters) ||
    throw(ArgumentError("Given system does not have the same number of parameter values provided as parameters."))

    vars = variables(F)
    params = ModelKit.parameters(F)
    D = ModelKit.degree(F.expressions, vars)

    t = ModelKit.unique_variable(:t, vars, params)

    γ = ModelKit.unique_variable(:γ, vars, params)
    all_params = [params; γ]
    _s_ = Variable[]
    for i = 1:n
        _s_i = ModelKit.unique_variable(Symbol(:_s_, i), vars, all_params)
        push!(_s_, _s_i)
        push!(all_params, _s_i)
    end

    f = F.expressions
    h = γ .* t .* _s_ .* (vars .^ D .- 1) .+ (1 .- t) .* f
    if scaling == nothing
        s = ones(n)
        H = ModelKitHomotopy(
            Homotopy(h, vars, t, all_params),
            [target_parameters; gamma; s],
        )
        u = zeros(ComplexF64, n)
        for i in 1:3
            evaluate!(u, H, LA.normalize!(randn(ComplexF64, n)), 0.0)
            s .= max.(s, exp2.(last.(frexp.(fast_abs.(u)))))
        end
        H.parameters[end-n+1:end] .= s
    else
        H = ModelKitHomotopy(
            Homotopy(h, vars, t, all_params),
            [target_parameters; gamma; scaling],
        )
    end

    S = TotalDegreeStarts(D)
    (homotopy = H, start_solutions = S)
end
