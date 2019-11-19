export ModelKitHomotopy

"""
    ModelKitHomotopy(system::ModelKit.System) <: AbstractSystem
"""
struct ModelKitHomotopy{T,P<:Union{Nothing,AbstractVector}} <: AbstractHomotopy
    homotopy::ModelKit.CompiledHomotopy{T}
    parameters::P
end

function ModelKitHomotopy(H::ModelKit.Homotopy; parameters = nothing)
    ModelKitHomotopy(ModelKit.compile(H), parameters)
end

Base.size(H::ModelKitHomotopy) = size(H.homotopy)

"""
    set_parameters!(H::ModelKitHomotopy, p::Tuple, γ = nothing)

Update the parameters of `H`.
"""
function set_parameters!(H::ModelKitHomotopy, p::Tuple, ::Nothing = nothing)
    p1, p2 = p
    set_start_parameters!(H, p1)
    set_target_parameters!(H, p2)

    n = length(p1)
    for i = 1:n
        H.parameters[i] = p1[i]
    end
    for i = 1:n
        H.parameters[i+n] = p2[i]
    end
    H.p[1] .= p[1]
    H.p[2] .= p[2]
    H.γ = γ
    H
end

function set_parameters!(H::ModelKitHomotopy, p::Tuple, γ::Tuple)
    p1, p2 = p
    for i = 1:n
        H.parameters[i] = p1[i]
    end
    for i = 1:n
        H.parameters[i+n] = p2[i]
    end

    γ₁, γ₂ = γ
    H.parameters[2n+1] = γ₁
    H.parameters[2n+2] = γ₂

    H
end


"""
    set_start_parameters!(H::ModelKitHomotopy, p)

Update the start parameters of `H`. This assumes that the first half of the parameters
are the start parameters.
"""
function set_start_parameters!(H::ModelKitHomotopy, p)
    for i = 1:div(length(p), 2)
        H.parameters[i] = p[i]
    end
    H
end

"""
    set_target_parameters!(H::ModelKitHomotopy, p)

Update the target parameters of `H`. This assumes that the second half of the parameters
are the start parameters.
"""
function set_target_parameters!(H::ModelKitHomotopy, p)
    n = div(length(p), 2)
    for i = 1:n
        H.parameters[i+n] = p[i]
    end
    H
end

cache(H::ModelKitHomotopy, x, t, p = nothing) = HomotopyNullCache()

function evaluate!(u, H::ModelKitHomotopy, x, t, ::HomotopyNullCache)
    ModelKit.evaluate!(u, H.homotopy, x, t, H.parameters)
end
function evaluate(H::ModelKitHomotopy, x, t, ::HomotopyNullCache)
    convert(Vector, ModelKit.evaluate(H.homotopy, x, t, H.parameters))
end
(H::ModelKitHomotopy)(x, t) = evaluate(F, x, t, H.parameters)

function jacobian!(U, H::ModelKitHomotopy, x, t, ::HomotopyNullCache)
    ModelKit.jacobian!(U, H.homotopy, x, t, H.parameters)
end
function jacobian(H::ModelKitHomotopy, x, t, ::HomotopyNullCache)
    convert(Matrix, ModelKit.jacobian(H.homotopy, x, t, H.parameters))
end

function evaluate_and_jacobian!(u, U, H::ModelKitHomotopy, x, t, ::HomotopyNullCache)
    ModelKit.evaluate_and_jacobian!(u, U, H.homotopy, x, t, H.parameters)
end

function dt!(u, H::ModelKitHomotopy, x, t, ::HomotopyNullCache)
    ModelKit.dt!(u, H.homotopy, x, t, H.parameters)
end
function dt(H::ModelKitHomotopy, x, t, ::HomotopyNullCache)
    convert(Vector, ModelKit.dt(H.homotopy, x, t, H.parameters))
end

function jacobian_and_dt!(U, u, H::ModelKitHomotopy, x, t, ::HomotopyNullCache)
    ModelKit.jacobian_and_dt!(U, u, H.homotopy, x, t, H.parameters)
end
function jacobian_and_dt(H::ModelKitHomotopy, x, t, ::HomotopyNullCache)
    ModelKit.jacobian_and_dt(H.homotopy, x, t, H.parameters)
end
