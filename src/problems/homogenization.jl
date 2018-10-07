export AbstractHomogenization, NullHomogenization, Homogenization

# The HomogenizationStrategy will be useful if we have multi-homogenous stuff
"""
    AbstractHomogenization

Abstract type for describing the homogenization strategy.
"""
abstract type AbstractHomogenization end

"""
    NullHomogenization()

Homogenization was not necessary since it is already homogenous.
"""
struct NullHomogenization <: AbstractHomogenization end

"""
    Homogenization(homvaridx::Int)

Homogenization happend by adding a new variable with index `homvaridx`.
"""
struct Homogenization <: AbstractHomogenization
    homvaridx::Int
end
function Homogenization(homvar::MP.AbstractVariable, variables::Vector{<:MP.AbstractVariable})
    homvaridx = findfirst(var -> var == homvar, variables)
    if homvaridx === nothing
        throw(error("homvar $homvar doesn't occur in the variables!"))
    end
    Homogenization(homvaridx)
end

homogenization(::Nothing) = NullHomogenization()
function homogenization(homvar::MP.AbstractVariable, variables)
    Homogenization(homvar, variables)
end
homogenization(homvar::Integer) = Homogenization(homvar)

ishomogenized(::Homogenization) = true
ishomogenized(::NullHomogenization) = false
