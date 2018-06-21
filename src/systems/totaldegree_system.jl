export TotalDegreeSystem

"""
    TotalDegreeSystem(polynomials, vars) <: AbstractSystem

Create a system using the `StaticPolynomials` package.
"""
struct TotalDegreeSystem <: AbstractSystem
    degrees::Vector{Int}
    degree_idxs::Vector{Int}
    hom_idx::Int

    function TotalDegreeSystem(degrees::Vector{Int}, degree_idxs::Vector{Int}, hom_idx::Int)
        @assert length(degrees) == length(degree_idxs)
        @assert 1 ≤ hom_idx ≤ length(degrees) + 1
        new(degrees, degree_idxs, hom_idx)
    end
end

function TotalDegreeSystem(F::Vector{<:MP.AbstractPolynomialLike}, vars, homvar=nothing)
    TotalDegreeSystem(MP.maxdegree.(F), vars, homvar)
end
function TotalDegreeSystem(degrees::Vector{Int}, vars::Vector{<:MP.AbstractVariable}, homvar::MP.AbstractVariable)
    degree_idxs = Int[]
    hom_idx = 0
    for (i, var) in enumerate(vars)
        if var == homvar
            hom_idx = i
        else
            push!(degree_idxs, i)
        end
    end
    TotalDegreeSystem(degrees, degree_idxs, hom_idx)
end

Base.size(F::TotalDegreeSystem) = (length(F.degrees), length(F.degrees) + 1)

function evaluate!(u, F::TotalDegreeSystem, x)
    for i=1:length(F.degrees)
        d = F.degrees[i]
        dix = F.degree_idxs[i]
        u[i] = x[dix]^d - x[F.hom_idx]^d
    end
    u
end
function evaluate(F::TotalDegreeSystem, x)
    u = similar(x, size(F, 1))
    evaluate!(u, F, x)
    u
end

function jacobian!(U, F::TotalDegreeSystem, x)
    U .= zero(eltype(x))
    hidx = F.hom_idx
    for i=1:length(F.degrees)
        d = F.degrees[i]
        didx = F.degree_idxs[i]
        if d == 1
            U[i, didx] = one(eltype(x))
            U[i, hidx] = -one(eltype(x))
        elseif d > 1
            U[i, didx] = d * x[didx]^(d-1)
            U[i, hidx] = -d * x[hidx]^(d-1)
        end
    end
    U
end
function jacobian(F::TotalDegreeSystem, x)
    U = similar(x, size(F))
    jacobian!(U, F, x)
    U
end

function evaluate_and_jacobian!(u, U, F::TotalDegreeSystem, x)
    U .= zero(x)
    hidx = F.hom_idx
    for i=1:length(F.degrees)
        d = F.degrees[i]
        didx = F.degree_idxs[i]
        if d == 1
            u[i] = x[didx] - x[hidx]
            U[i, didx] = one(x)
            U[i, hidx] = -one(x)
        elseif d > 1
            xd = x[didx]^(d-1)
            xh = x[hidx]^(d-1)
            u[i] = xd * x[didx] - xh * x[hidx]
            U[i, didx] = d * xd
            U[i, hidx] = -d * xh
        end
    end
    nothing
end

function evaluate_and_jacobian(F::TotalDegreeSystem, x)
    u = similar(x, size(F, 1))
    U = similar(x, size(F))
    evaluate_and_jacobian!(u, U, F, x)
    u, U
end
