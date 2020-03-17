import Base.GMP: MPZ
import MixedSubdivisions: MixedCell

struct BinomialSystemSolver
    A::Matrix{Int32} # binomial system lhs
    b::Vector{ComplexF64} # binomial system rhs
    X::ElasticArray{ComplexF64,2,1} # Solutions
    # HNF
    H::Matrix{Int64} # hermite normal form
    U::Matrix{Int64} # trafo matrix for hnf
    H_big::Matrix{BigInt} # hermite normal form
    U_big::Matrix{BigInt} # trafo matrix for hnf
    # Solve data structures
    γ::Vector{Float64} # angle of the elements in b
    μ::Vector{Float64} # γ^U
    αs::Vector{Float64} # partial angles in the triangular solve
    Aᵀ::Matrix{Float64} # transposed of A to solve coords
    unit_roots_table::ElasticArray{Int32,2,1}
end

function BinomialSystemSolver(A::Matrix, b::Vector)
    BSS = BinomialSystemSolver(size(A, 1))
    BSS.A .= A
    BSS.b .= b
    BSS
end

function BinomialSystemSolver(n::Int)
    A = zeros(Int32, n, n)
    b = zeros(ComplexF64, n)
    X = ElasticArray{ComplexF64}(undef, n, 0)
    H = zeros(Int64, n, n)
    U = zeros(Int64, n, n)
    H_big = Matrix{BigInt}(undef, n, n)
    H_big .= zero.(BigInt)
    U_big = Matrix{BigInt}(undef, n, n)
    U_big .= zero.(BigInt)
    γ = zeros(Float64, n)
    μ = zeros(Float64, n)
    αs = zeros(Float64, n)
    Aᵀ = zeros(Float64, n, n)
    S = ElasticArray{Int32}(undef, n, 0)
    BinomialSystemSolver(A, b, X, H, U, H_big, U_big, γ, μ, αs, Aᵀ, S)
end

function init!(BSS::BinomialSystemSolver, support, coeffs, cell::MixedCell)
    for (i, (aᵢ, bᵢ)) in enumerate(cell.indices)
        for j = 1:size(BSS.A, 1)
            BSS.A[j, i] = support[i][j, aᵢ] - support[i][j, bᵢ]
        end
        BSS.b[i] = -coeffs[i][bᵢ] / coeffs[i][aᵢ]
    end
    BSS
end

function fill_unit_roots_combinations!(S::ElasticArray, H, d̂)
    n = size(H, 1)
    resize!(S, n, d̂)
    d, e = d̂, 1
    @inbounds for i = 1:n
        dᵢ = Int64(H[i, i])
        d = d ÷ dᵢ
        k = 1
        for _ = 1:e, j = 0:(dᵢ-1), _ = 1:d
            S[i, k] = j
            k += 1
        end
        e *= dᵢ
    end
    S
end

function compute_angular_part!(
    BSS::BinomialSystemSolver,
    H::Matrix{Int64},
    U::Matrix{Int64},
    d̂::Int,
)
    @unpack X, γ, b, μ, unit_roots_table, αs = BSS
    n = size(H, 1)

    # apply coordinate change to b
    γ .= angle.(b) ./ 2π
    μ .= zero(eltype(μ))
    for j = 1:n
        for i = 1:n
            μ[j] += rem(U[i, j] * γ[i], 2.0, RoundNearest)
        end
        μ[j] = rem(μ[j], 2.0, RoundNearest)
    end
    # solve triangular system
    @inbounds for i = 1:d̂, j = n:-1:1
        α = μ[j] + unit_roots_table[j, i]
        α = rem((μ[j] + unit_roots_table[j, i]) / H[j, j], 2.0, RoundNearest)
        for k = n:-1:(j+1)
            α -= rem(αs[k] * H[k, j] / H[j, j], 2.0, RoundNearest)
            α = rem(α, 2.0, RoundNearest)
        end
        # α = rem(α, 2.0, RoundNearest)
        X[j, i] = cis(2π * α)
        αs[j] = α
    end
    BSS
end


function compute_angular_part!(
    BSS::BinomialSystemSolver,
    H::Matrix{BigInt},
    U::Matrix{BigInt},
    d̂::Int,
)
    # @unpack X, γ, b, μ, unit_roots_table = BSS
    @unpack X, unit_roots_table, b = BSS
    n = size(H, 1)

    # compute precision necessary
    p = max(maximum(x -> MPZ.sizeinbase(x, 2), H), maximum(x -> MPZ.sizeinbase(x, 2), U))
    prec = max((ceil(Int, (p + 53) / 32)) * 32, 64)

    γ = angle.(b) ./ 2π

    μ = [BigFloat(0.0; precision = prec) for i in 1:n]
    α = γᵢ = BigFloat(0.0; precision = prec)
    αk = γᵢⱼ = BigFloat(0.0; precision = prec)
    m = BigFloat(2.0; precision = prec)

    for i = 1:n
        set!(γᵢ, γ[i])
        for j = 1:n
            mul!(γᵢⱼ, γᵢ, U[i, j])
            if γᵢⱼ < -1 || γᵢⱼ > 1
                rem!(γᵢⱼ, γᵢⱼ, m, RoundNearest)
            end
            add!(μ[j], μ[j], γᵢⱼ)
        end
    end
    αs = [BigFloat(0.0; precision = prec) for i in 1:n]
    # solve triangular system
    for i = 1:d̂, j = n:-1:1
        add!(α, μ[j], Int64(unit_roots_table[j, i]))
        div!(α, α, H[j,j])
        for k = n:-1:(j+1)
            mul!(αk, αs[k], H[k, j])
            div!(αk, αk, H[j, j])
            sub!(α, α, αk)
        end
        if α < -1 || α > 1
            rem!(α, α, m, RoundNearest)
        end
        X[j, i] = cis(2π * Float64(α))
        set!(αs[j], α)
    end
    BSS
end

function solve(BSS::BinomialSystemSolver, support, coeffs, cell::MixedCell)
    init!(BSS, support, coeffs, cell)
    solve!(BSS)
end

function solve(BSS::BinomialSystemSolver, A, b)
    BSS.A .= A
    BSS.b .= b
    solve!(BSS)
end

function solve!(BSS::BinomialSystemSolver)
    try
        hnf!(BSS.H, BSS.U, BSS.A)
        solve!(BSS, BSS.H, BSS.U)
    catch e
        isa(e, OverflowError) || rethrow(e)
        hnf!(BSS.H_big, BSS.U_big, BSS.A)
        solve!(BSS, BSS.H_big, BSS.U_big)
    end
    BSS.X
end

function solve!(BSS::BinomialSystemSolver, H::Matrix{<:Integer}, U::Matrix{<:Integer})
    n = size(H, 1)
    d̂ = 1
    for i = 1:n
        d̂ *= Int64(H[i, i])
    end
    resize!(BSS.X, n, d̂)
    fill_unit_roots_combinations!(BSS.unit_roots_table, H, d̂)

    # We split the computation of the solution in 2 stages
    # 1) Compute solutions for the angular part -> d̂ many solutions
    compute_angular_part!(BSS, H, U, d̂)

    @unpack A, Aᵀ, b, μ, X = BSS
    # 2) Solve for the absolute value
    μ .= log.(fast_abs.(b))
    LA.transpose!(Aᵀ, A)
    LA.ldiv!(LA.lu!(Aᵀ), μ)
    μ .= exp.(μ)
    for j = 1:d̂, i = 1:n
        X[i, j] *= μ[i]
    end
    X
end

#########################
## Hermite Normal Form ##
#########################

"""
    hnf(A, T=elytpe(A))
Compute the hermite normal form `H` of `A` by overwriting `A` and the corresponding transformation
matrix `U` using precision T. This is `A*U == H` and `H` is a lower triangular matrix.
The implementation follows the algorithm described in [1].

[1] Kannan, Ravindran, and Achim Bachem. "Polynomial algorithms for computing the Smith and Hermite normal forms of an integer matrix." SIAM Journal on Computing 8.4 (1979): 499-507.
"""
function hnf(A, T = eltype(A))
    H = similar(A, T)
    H .= A
    U = similar(A, T)
    U .= zero.(T)
    hnf!(H, U)
    H, U
end

# use checked arithmethic
⊡(x, y) = Base.checked_mul(x, y)
⊞(x, y) = Base.checked_add(x, y)

"""
    hnf!(H, U, A)
Inplace version of [hnf](@ref).

    hnf!(A, U)
Inplace version of [hnf](@ref) overwriting `A` with `H`.
"""
hnf!(H, U, A) = hnf!(copyto!(H, A), U)
function hnf!(A, U)
    n = size(A, 1)
    U .= 0
    @inbounds for i = 1:n
        U[i, i] = one(eltype(U))
    end
    @inbounds for i = 1:(n-1)
        ii = i ⊞ 1
        for j = 1:i
            if !iszero(A[j, j]) || !iszero(A[j, ii])
                # 4.1
                r, p, q = gcdx(A[j, j], A[j, ii])
                # 4.2
                d_j = -A[j, ii] ÷ r
                d_ii = A[j, j] ÷ r
                for k = 1:n
                    a_kj, a_kii = A[k, j], A[k, ii]
                    A[k, j] = a_kj ⊡ p ⊞ a_kii ⊡ q
                    A[k, ii] = a_kj ⊡ d_j ⊞ a_kii ⊡ d_ii

                    u_kj, u_kii = U[k, j], U[k, ii]
                    U[k, j] = u_kj ⊡ p ⊞ u_kii ⊡ q
                    U[k, ii] = u_kj ⊡ d_j ⊞ u_kii ⊡ d_ii
                end
            end
            # 4.3
            if j > 1
                reduce_off_diagonal!(A, U, j)
            end
        end
        #5
        reduce_off_diagonal!(A, U, ii)
    end
    # Special case for 1 × 1 matrices to guarantee that the diagonal is positive
    # This comes up in polyhedral homotopy for univariate polynomials
    if n == 1
        U[1, 1] = flipsign(U[1, 1], A[1, 1])
        A[1, 1] = flipsign(A[1, 1], A[1, 1])
    end

    nothing
end

function cdiv(x, y)
    d = ceil(typeof(x), x / y)
    if 0 ≥ x - d * y > -y
        return d
    else
        d, r = divrem(x, y)
        d + (r > 0)
    end
end


@inline function reduce_off_diagonal!(A, U, k)
    n = size(A, 1)
    @inbounds if A[k, k] < 0
        for i = 1:n
            A[i, k] = -A[i, k]
            U[i, k] = -U[i, k]
        end
    end
    @inbounds for z = 1:(k-1)
        if !iszero(A[z, z])
            # r = -ceil(eltype(A), big(A[k, z]) / A[z, z])
            r = -cdiv(A[k, z], A[z, z])
            for i = 1:n
                U[i, z] = U[i, z] ⊞ r ⊡ U[i, k]
                A[i, z] = A[i, z] ⊞ r ⊡ A[i, k]
            end
        end
    end
    nothing
end

function hnf!(A::AbstractMatrix{BigInt}, U::AbstractMatrix{BigInt})
    n = size(A, 1)
    for j = 1:n, i = 1:n
        if i == j
            MPZ.set_si!(U[i, j], 1)
        else
            MPZ.set_si!(U[i, j], 0)
        end
    end
    r, p, q = zero(BigInt), zero(BigInt), zero(BigInt)
    d_j, d_ii = zero(BigInt), zero(BigInt)
    t1, t2, t3 = zero(BigInt), zero(BigInt), zero(BigInt)
    @inbounds for i = 1:(n-1)
        ii = i + 1
        for j = 1:i
            if !iszero(A[j, j]) || !iszero(A[j, ii])
                # 4.1
                MPZ.gcdext!(r, p, q, A[j, j], A[j, ii])
                # 4.2
                # d_j = -A[j, ii] ÷ r
                MPZ.mul_si!(d_j, A[j, ii], -1)
                MPZ.tdiv_q!(d_j, d_j, r)

                # d_ii = A[j, j] ÷ r
                MPZ.tdiv_q!(d_ii, A[j, j], r)
                for k = 1:n
                    # a_kj, a_kii = A[k, j], A[k, ii]
                    a_kj = MPZ.set!(t1, A[k, j])
                    a_kii = MPZ.set!(t2, A[k, ii])
                    # A[k, j] = a_kj ⊡ p ⊞ a_kii ⊡ q
                    MPZ.mul!(A[k, j], a_kj, p)
                    MPZ.mul!(t3, a_kii, q)
                    MPZ.add!(A[k, j], A[k, j], t3)
                    # A[k, ii] = a_kj ⊡ d_j ⊞ a_kii ⊡ d_ii
                    MPZ.mul!(A[k, ii], a_kj, d_j)
                    MPZ.mul!(t3, a_kii, d_ii)
                    MPZ.add!(A[k, ii], A[k, ii], t3)

                    # u_kj, u_kii = U[k, j], U[k, ii]
                    u_kj = MPZ.set!(t1, U[k, j])
                    u_kii = MPZ.set!(t2, U[k, ii])
                    # U[k, j] = u_kj ⊡ p ⊞ u_kii ⊡ q
                    MPZ.mul!(U[k, j], u_kj, p)
                    MPZ.mul!(t3, u_kii, q)
                    MPZ.add!(U[k, j], U[k, j], t3)

                    # U[k, ii] = u_kj ⊡ d_j ⊞ u_kii ⊡ d_ii
                    MPZ.mul!(U[k, ii], u_kj, d_j)
                    MPZ.mul!(t3, u_kii, d_ii)
                    MPZ.add!(U[k, ii], U[k, ii], t3)
                end
            end
            # 4.3
            if j > 1
                reduce_off_diagonal!(A, U, j, t1, t2)
            end
        end
        #5
        reduce_off_diagonal!(A, U, ii, t1, t2)
    end
    # Special case for 1 × 1 matrices to guarantee that the diagonal is positive
    # This comes up in polyhedral homotopy for univariate polynomials
    if n == 1
        MPZ.mul_si!(U[1, 1], A[1, 1] < 0 ? -1 : 1)
        MPZ.mul_si!(A[1, 1], A[1, 1] < 0 ? -1 : 1)
    end

    nothing
end

cdiv_q!(x::BigInt, a::BigInt, b::BigInt) = (
    ccall(
        (:__gmpz_cdiv_q, :libgmp),
        Cvoid,
        (MPZ.mpz_t, MPZ.mpz_t, MPZ.mpz_t),
        x,
        a,
        b,
    ); x
)

function reduce_off_diagonal!(
    A::AbstractMatrix{BigInt},
    U::AbstractMatrix{BigInt},
    k,
    t1,
    t2,
)
    n = size(A, 1)
    @inbounds if A[k, k] < 0
        for i = 1:n
            MPZ.mul_si!(A[i, k], -1)
            MPZ.mul_si!(U[i, k], -1)
        end
    end

    @inbounds for z = 1:(k-1)
        if !iszero(A[z, z])
            r = cdiv_q!(t1, A[k, z], A[z, z])
            r = MPZ.mul_si!(r, r, -1)
            # r = -divceil(A[k, z], A[z,z])
            for i = 1:n
                MPZ.mul!(t2, r, U[i, k])
                MPZ.add!(U[i, z], U[i, z], t2)
                # U[i, z] = U[i, z] ⊞ r ⊡ U[i, k]
                MPZ.mul!(t2, r, A[i, k])
                MPZ.add!(A[i, z], A[i, z], t2)
                # A[i, z] = A[i, z] ⊞ r ⊡ A[i, k]
            end
        end
    end
    nothing
end
