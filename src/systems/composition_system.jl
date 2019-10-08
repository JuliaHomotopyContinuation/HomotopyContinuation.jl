export CompositionSystem

"""
    CompositionSystem(composition::Composition, systems_constructor) <: AbstractSystem

A system representing the composition of polynomial maps.
"""
struct CompositionSystem{S1<:AbstractSystem,S2<:AbstractSystem} <: AbstractSystem
    # The system is g ∘ f
    g::S2 # Never a composition system
    f::S1 # Can be a CompositionSystem again

    g_has_parameters::Bool
    f_has_parameters::Bool
end

function CompositionSystem(
    C::Composition,
    system_constructor;
    variables = nothing,
    parameters = nothing,
    homvars = nothing,
)
    n = length(C.polys)
    f = CompositionSystem(
        C.polys[n-1],
        C.polys[n],
        system_constructor;
        input_variables = variables,
        parameters = parameters,
        homvars = homvars,
    )
    for k = (n - 2):-1:1
        f = CompositionSystem(
            C.polys[k],
            f,
            system_constructor;
            parameters = parameters,
            homvars = homvars,
        )
    end
    f
end

function CompositionSystem(
    g::Vector{<:MP.AbstractPolynomialLike},
    f::Union{CompositionSystem,Vector{<:MP.AbstractPolynomialLike}},
    system_constructor;
    input_variables = nothing,
    parameters = nothing,
    homvars = nothing,
)

    vars = variables(g; parameters = parameters)
    homvars_to_end!(vars, homvars)
    G = system_constructor(g, variables = vars, parameters = parameters)

    if isa(f, CompositionSystem)
        CompositionSystem(G, f, hasparameters(g, parameters), hasparameters(f))
    else
        F = system_constructor(f, variables = input_variables, parameters = parameters)
        CompositionSystem(G, F, hasparameters(g, parameters), hasparameters(f, parameters))
    end
end


hasparameters(C::CompositionSystem) = C.g_has_parameters || C.f_has_parameters

homvars_to_end!(variables, ::Nothing) = variables
homvars_to_end!(vars, homvars::MP.AbstractVariable) = homvars_to_end!(vars, (homvars,))
function homvars_to_end!(vars, homvars)
    deleteat!(vars, findall(v -> v ∈ homvars, vars))
    append!(vars, homvars)
    vars
end

struct CompositionSystemCache{
    C1<:AbstractSystemCache,
    C2<:AbstractSystemCache,
    T,
} <: AbstractSystemCache
    cache_f::C1
    cache_g::C2

    eval_f::Vector{T}
    J_f::Matrix{T}
    J_g::Matrix{T}
    Jp_f::Union{Nothing,Matrix{T}}
    Jp_g::Union{Nothing,Matrix{T}}
end

function cache(C::CompositionSystem, x, p = nothing)
    f, g = C.f, C.g
    c_f = p === nothing ? cache(f, x) : cache(f, x, p)
    eval_f = p === nothing ? evaluate(f, x, c_f) : evaluate(f, x, p, c_f)
    c_g = p === nothing ? cache(g, eval_f) : cache(g, eval_f, p)
    J_f = p === nothing ? jacobian(f, x, c_f) : jacobian(f, x, p, c_f)
    J_g = p === nothing ? jacobian(g, eval_f, c_g) : jacobian(g, eval_f, p, c_g)
    # need to bring eval_f, J_f, J_g to the same element type
    T = promote_type(eltype(eval_f), eltype(J_f), eltype(J_g))
    eval_f_T = convert(Vector{T}, eval_f)
    J_f_T = convert(Matrix{T}, J_f)
    J_g_T = convert(Matrix{T}, J_g)
    if p === nothing
        Jp_f = Jp_g = nothing
    else
        Jp_f = similar(J_f_T, size(J_f_T, 1), length(p))
        Jp_g = similar(J_g_T, size(J_g_T, 1), length(p))
    end

    CompositionSystemCache(c_f, c_g, eval_f_T, J_f_T, J_g_T, Jp_f, Jp_g)
end

Base.size(C::CompositionSystem) = (length(C.g), size(C.f)[2])

function evaluate!(u, C::CompositionSystem, x, c::CompositionSystemCache)
    evaluate!(u, C.g, evaluate!(c.eval_f, C.f, x, c.cache_f), c.cache_g)
end
function evaluate!(u, C::CompositionSystem, x, p, c::CompositionSystemCache)
    evaluate!(u, C.g, evaluate!(c.eval_f, C.f, x, p, c.cache_f), p, c.cache_g)
end
function evaluate(C::CompositionSystem, x, c::CompositionSystemCache)
    evaluate(C.g, evaluate!(c.eval_f, C.f, x, c.cache_f), c.cache_g)
end
function evaluate(C::CompositionSystem, x, p, c::CompositionSystemCache)
    evaluate(C.g, evaluate!(c.eval_f, C.f, x, p, c.cache_f), p, c.cache_g)
end
function jacobian!(U, C::CompositionSystem, x, c::CompositionSystemCache)
    # chain rule
    evaluate_and_jacobian!(c.eval_f, c.J_f, C.f, x, c.cache_f)
    jacobian!(c.J_g, C.g, c.eval_f, c.cache_g)
    LinearAlgebra.mul!(U, c.J_g, c.J_f)
end
function jacobian!(U, C::CompositionSystem, x, p, c::CompositionSystemCache)
    # chain rule
    evaluate_and_jacobian!(c.eval_f, c.J_f, C.f, x, p, c.cache_f)
    jacobian!(c.J_g, C.g, c.eval_f, p, c.cache_g)
    LinearAlgebra.mul!(U, c.J_g, c.J_f)
end
function jacobian(C::CompositionSystem, x, c::CompositionSystemCache)
    jacobian!(similar(c.J_g, size(c.J_g, 1), size(c.J_f, 2)), C, x, c)
end
function jacobian(C::CompositionSystem, x, p, c::CompositionSystemCache)
    jacobian!(similar(c.J_g, size(c.J_g, 1), size(c.J_f, 2)), C, x, p, c)
end
function evaluate_and_jacobian!(u, U, C::CompositionSystem, x, c::CompositionSystemCache)
    evaluate_and_jacobian!(c.eval_f, c.J_f, C.f, x, c.cache_f)
    evaluate_and_jacobian!(u, c.J_g, C.g, c.eval_f, c.cache_g)
    LinearAlgebra.mul!(U, c.J_g, c.J_f)
    nothing
end
function evaluate_and_jacobian!(u, U, C::CompositionSystem, x, p, c::CompositionSystemCache)
    evaluate_and_jacobian!(c.eval_f, c.J_f, C.f, x, p, c.cache_f)
    evaluate_and_jacobian!(u, c.J_g, C.g, c.eval_f, p, c.cache_g)
    LinearAlgebra.mul!(U, c.J_g, c.J_f)
    nothing
end

function differentiate_parameters!(U, C::CompositionSystem, x, p, c::CompositionSystemCache)
    @assert c.Jp_f !== nothing && c.Jp_g !== nothing

    # d/(dp)(g(f(x, p), p)) = g_y(f(x, p), p)f_p(x, p) + g_p(f(x, p), p)
    evaluate!(c.eval_f, C.f, x, p, c.cache_f)
    # See whether we can simplify things
    if C.g_has_parameters && C.f_has_parameters
        differentiate_parameters!(c.Jp_g, C.g, c.eval_f, p, c.cache_g)
        differentiate_parameters!(c.Jp_f, C.f, x, p, c.cache_f)
        jacobian!(c.J_g, C.g, c.eval_f, p, c.cache_g)
        LinearAlgebra.mul!(U, c.J_g, c.Jp_f)
        U .= U .+ c.Jp_g
    elseif C.g_has_parameters # && !C.f_has_parameters
        differentiate_parameters!(c.Jp_g, C.g, c.eval_f, p, c.cache_g)
        U .= c.Jp_g
    elseif C.f_has_parameters  # !C.g_has_parameters &&
        differentiate_parameters!(c.Jp_f, C.f, x, p, c.cache_f)
        jacobian!(c.J_g, C.g, c.eval_f, p, c.cache_g)
        LinearAlgebra.mul!(U, c.J_g, c.Jp_f)
    else
        U .= zero(eltype(U))
    end

    U
end

function differentiate_parameters(F::CompositionSystem, x, p, c::CompositionSystemCache)
    U = similar(c.J_g, size(c.J_g, 1), size(c.Jp_f, 2))
    differentiate_parameters!(U, F, x, p, c)
end
