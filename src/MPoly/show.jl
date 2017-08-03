function show{T,N}(io::IO, p::Poly{T,N})
    first = true
    for exp_ in sorted_exponents(p)
        coeff = p[exp_]
        if (!first && show_plus(coeff))
            Base.print(io, "+")
        end
        first = false

        if (coeff != one(T) && coeff != -one(T)) || exp_ == ntuple(_ -> 0, N)
            show_coeff(io, coeff)
        elseif coeff == -one(T)
            Base.print(io, "-")
        end

        for (var, power) in zip(p.vars, exp_)
            if power == 1
                Base.print(io, var)
            elseif power > 1
                Base.print(io, "$(var)^$(power)")
            end
        end
    end

    if first
        Base.print(io, zero(T))
    end
end

function show{T}(io::IO, p::FixedPoly{T})
    first = true
    exps = exponents(p)
    cfs = coeffs(p)

    m, n = size(exps)

    for i=1:n
        exp = exps[:, i]
        coeff = cfs[i]

        if (!first && show_plus(coeff))
            Base.print(io, "+")
        end
        first = false

        if (coeff != one(T) && coeff != -one(T)) || exp == zeros(Int, m)
            show_coeff(io, coeff)
        elseif coeff == -one(T)
            Base.print(io, "-")
        end

        for (var, power) in zip(vars(p), exp)
            if power == 1
                Base.print(io, var)
            elseif power > 1
                Base.print(io, "$(var)^$(power)")
            end
        end
    end

    if first
        Base.print(io, zero(T))
    end
end

# helpers
show_plus{T<:Real}(x::T) = x >= zero(T)
show_plus{T<:Complex}(x::T) = x == -one(T) ? false : true

show_coeff{T<:Real}(io::IO, x::T) = Base.print(io, x)
function show_coeff{T<:Complex}(io::IO, x::T)
    if imag(x) ≈ 0
        Base.print(io, convert(Float64, x))
    else
        Base.print(io, "($(x))")
    end
end

function total_deg_isless{N,T<:Number}(a::NTuple{N,T},b::NTuple{N,T})
    if sum(a) < sum(b)
        true
    else
        a < b
    end
end
sort_bytotal{T<:Number, N}(xs::Vector{NTuple{N,T}}) = sort(xs, lt=total_deg_isless, rev=true)
sort_bytotal!{T<:Number, N}(xs::Vector{NTuple{N,T}}) = sort!(xs, lt=total_deg_isless, rev=true)
sorted_exponents{T,N}(p::Poly{T,N})::Vector{NTuple{N,Int}} = sort_bytotal(collect(keys(p.terms)))
