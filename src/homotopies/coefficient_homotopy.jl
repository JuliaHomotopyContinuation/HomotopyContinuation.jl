export CoefficientHomotopy


struct CoefficientHomotopy{S} <: AbstractHomotopy
    system::ModelKit.CompiledSystem{S}
    start_coeffs::Vector{ComplexF64}
    target_coeffs::Vector{ComplexF64}

    taylor_coeffs::NTuple{2,Vector{ComplexF64}}
    t_coeffs::Base.RefValue{ComplexF64}
    # these are just here to avoid unnecessary allocations
    dc1::Tuple{Vector{ComplexF64}}
end

function CoefficientHomotopy(
    system::ModelKit.CompiledSystem,
    start_coefficients::AbstractVector{<:AbstractVector{<:Number}},
    target_coefficients::AbstractVector{<:AbstractVector{<:Number}}
)
    m = ModelKit.nparameters(system)
    m1 = sum(length, start_coefficients)
    m2 = sum(length, target_coefficients)
    m == m1 == m2 ||
    throw(ArgumentError("System parameters and coefficients do not have the same size, got $m and $m1"))

    n = nvariables(system)

    start_coeffs = ComplexF64[]
    for c in start_coefficients
        append!(start_coeffs, c)
    end
    target_coeffs = ComplexF64[]
    for c in target_coefficients
        append!(target_coeffs, c)
    end
    taylor_coeffs = tuple((zeros(ComplexF64, m) for i = 0:1)...)
    dc1 = (taylor_coeffs[2],)


    CoefficientHomotopy(
        system,
        start_coeffs,
        target_coeffs,
        taylor_coeffs,
        Ref(Complex(NaN,NaN)),
        dc1,
    )
end

Base.size(H::CoefficientHomotopy) = size(H.system)

function taylor_coeffs!(H::CoefficientHomotopy, t::Number)
    H.t_coeffs[] != t || return H.taylor_coeffs

    c, c¹ = H.taylor_coeffs
    t0 = (1 - t)

    c .= t .* H.start_coeffs .+ t0 .* H.target_coeffs
    c¹ .= H.start_coeffs .- H.target_coeffs

    H.t_coeffs[] = t

    H.taylor_coeffs
end

function evaluate!(u, H::CoefficientHomotopy, x::AbstractVector, t)
    c, _ = taylor_coeffs!(H::CoefficientHomotopy, t)
    ModelKit.evaluate!(u, H.system, x, c)
end

function evaluate_and_jacobian!(u, U, H::CoefficientHomotopy, x::AbstractVector, t)
    c, _ = taylor_coeffs!(H::CoefficientHomotopy, t)
    ModelKit.evaluate_and_jacobian!(u, U, H.system, x, c)
    nothing
end

function diff_t!(u, H::CoefficientHomotopy, x, t, dx::Tuple)
    c, _ = taylor_coeffs!(H::CoefficientHomotopy, t)
    ModelKit.diff_t!(u, H.system, x, dx, c, H.dc1)
    u
end
