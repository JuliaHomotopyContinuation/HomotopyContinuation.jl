import LinearAlgebra
using LinearAlgebra: BLAS, LAPACK, BlasInt
using LinearAlgebra.BLAS: @blasfunc, libblastrampoline
import RecursiveFactorization


Base.@kwdef mutable struct LinearSolveWorkspaceStatistics
    factorizations::Int = 0
    solves::Int = 0
end
function init!(s::LinearSolveWorkspaceStatistics)
    s.factorizations = 0
    s.solves = 0
    s
end

Base.@kwdef struct LinearSolveWorkspace
    norm::WeightedNorm{InfNorm}
    # We solve A * u = b
    A::Matrix{ComplexF64}
    LU::Matrix{ComplexF64}
    b::Vector{ComplexF64}
    u::Vector{ComplexF64}
    x::Vector{ComplexF64}
    ipiv::Vector{Int}
    rcond::RefValue{Float64}
    ferr::RefValue{Float64}
    berr::RefValue{Float64}
    # row and column scalings
    R::Vector{Float64}
    C::Vector{Float64}
    rowcond::RefValue{Float64}
    colcond::RefValue{Float64}
    amax::RefValue{Float64}
    work_n::Vector{ComplexF64}
    work_2n::Vector{ComplexF64}
    rwork_n::Vector{Float64}
    rwork_2n::Vector{Float64}
    work_n_extended::Vector{ComplexDF64}
    weighted_norm::WeightedNorm{InfNorm}
    # state
    factorized::RefValue{Bool}
    solved::RefValue{Bool}
    equilibriate::RefValue{Bool}
    equilibriated::RefValue{Bool}
    statistics::LinearSolveWorkspaceStatistics = LinearSolveWorkspaceStatistics()
end

function LinearSolveWorkspace(
    A::Matrix{ComplexF64},
    b::Vector{ComplexF64},
    norm::WeightedNorm{InfNorm} = WeightedNorm(InfNorm(), size(A, 2));
    equilibriate::Bool = true,
)
    m, n = size(A)
    m == n || throw(DimensionMismatch("A must be square"))
    m == length(b) || throw(DimensionMismatch("A and b must have the same number of rows"))
    ws = LinearSolveWorkspace(
        norm = norm,
        A = copy(A),
        LU = Matrix{ComplexF64}(undef, m, n),
        b = copy(b),
        u = Vector{ComplexF64}(undef, n),
        x = Vector{ComplexF64}(undef, n),
        ipiv = Vector{Int}(undef, m),
        rcond = Ref{Float64}(),
        ferr = Ref{Float64}(),
        berr = Ref{Float64}(),
        R = Vector{Float64}(undef, m),
        C = Vector{Float64}(undef, n),
        rowcond = Ref{Float64}(),
        colcond = Ref{Float64}(),
        amax = Ref{Float64}(),
        work_n = Vector{ComplexF64}(undef, n),
        work_2n = Vector{ComplexF64}(undef, 2 * n),
        rwork_n = Vector{Float64}(undef, n),
        rwork_2n = Vector{Float64}(undef, 2 * n),
        work_n_extended = Vector{ComplexDF64}(undef, n),
        weighted_norm = WeightedNorm(Vector{Float64}(undef, n), InfNorm()),
        factorized = Ref{Bool}(false),
        solved = Ref{Bool}(false),
        equilibriate = Ref{Bool}(equilibriate),
        equilibriated = Ref{Bool}(false),
    )
    set!(ws, A, b)
end
LinearSolveWorkspace(
    m::Int,
    n::Int,
    norm::WeightedNorm{InfNorm} = WeightedNorm(InfNorm, n),
) = LinearSolveWorkspace(
    Matrix{ComplexF64}(undef, m, n),
    Vector{ComplexF64}(undef, n),
    norm,
)

solution(ws::LinearSolveWorkspace) = ws.x
Base.size(ws::LinearSolveWorkspace) = size(ws.A)
get_A(ws::LinearSolveWorkspace) = ws.A
get_b(ws::LinearSolveWorkspace) = ws.b

function Base.show(io::IO, ws::LinearSolveWorkspace)
    m, n = size(ws.A)
    print(io, "LinearSolveWorkspace(A,b) with size(A) = ($(m),$(n))")
end

function init!(ws::LinearSolveWorkspace)
    init!(ws.statistics)
    ws
end

function set!(ws::LinearSolveWorkspace, A::Matrix{ComplexF64}, b::Vector{ComplexF64})
    set_A!(ws, A)
    set_b!(ws, b)
    ws
end

function set_A!(ws::LinearSolveWorkspace, A::AbstractMatrix)
    A === ws.A || copy!(ws.A, A)
    ws.R .= 1
    ws.rowcond[] = 1
    ws.colcond[] = 1
    ws.amax[] = 1
    ws.rcond[] = NaN
    ws.factorized[] = false
    ws.equilibriated[] = false
    ws
end
function set_b!(ws::LinearSolveWorkspace, b::AbstractVector)
    b === ws.b || copy!(ws.b, b)
    ws.solved[] = false
    ws.ferr[] = 0
    ws.berr[] = 0
    if ws.factorized[] && ws.equilibriated[]
        m = length(ws.b)
        @inbounds for i = 1:m
            ws.b[i] *= ws.R[i]
        end
    end
    ws
end

function equilibriate!(ws::LinearSolveWorkspace)
    if ws.equilibriate[] && !ws.equilibriated[]
        m, n = size(ws.A)
        # (_, _, rowcond, colcond, amax) = equilibriate!(ws.R, ws.C, ws.A)
        ws.C .= weights(ws.norm)
        skeel_row_scaling!(ws.R, ws.A, ws.C)

        @inbounds for i = 1:n
            ws.b[i] *= ws.R[i]
        end

        @inbounds for j = 1:n
            cj = ws.C[j]
            for i = 1:m
                ws.A[i, j] = cj * ws.R[i] * ws.A[i, j]
            end
        end
        ws.equilibriated[] = true
    end
    ws
end

function factorize!(ws::LinearSolveWorkspace)
    ws.factorized[] && return true

    equilibriate!(ws)
    @inbounds copy!(ws.LU, ws.A)
    @inbounds RecursiveFactorization.lu!(ws.LU, ws.ipiv)

    ws.statistics.factorizations += 1
    ws.factorized[] = true
end

function rcond!(ws::LinearSolveWorkspace)
    try
        factorize!(ws)
        if isnan(ws.rcond[])
            ws.rcond[] = gecon!('I', ws.LU, norm(ws.A, Inf), ws.work_2n, ws.rwork_2n)
        end
    catch e
        if e isa LA.SingularException
            ws.rcond[] = 0.0
        else
            rethrow(e)
        end
    end

    ws.rcond[]
end

function transform_u_to_x!(ws::LinearSolveWorkspace)
    if ws.equilibriated[]
        @inbounds for i = 1:length(ws.u)
            ws.x[i] = ws.C[i] * ws.u[i]
        end
    else
        ws.x .= ws.u
    end
    ws.x
end

function solve!(ws::LinearSolveWorkspace)
    factorize!(ws)
    if (!ws.solved[])
        lu_ldiv!(ws.u, ws.LU, ws.ipiv, ws.b)
        transform_u_to_x!(ws)
        ws.statistics.solves += 1
        ws.solved[] = true
    end

    ws.x
end


Base.@kwdef struct RefinementResult
    ferr::Float64
    initial_ferr::Float64
    ω::Float64
    iters::Int
end
Base.show(io::IO, r::RefinementResult) = print_fieldnames(io, r)

"""
refine!(ws::LinearSolveWorkspace; max_iters=3, tol=0.0, extended_precision::Bool=false)

Perform one step of mixed precision iterative refinement.
If `norm` is an `AbstractNorm` then the normwise relative error before
the refinement step is returned otherwise `x` is returned.
"""
function refine!(
    ws::LinearSolveWorkspace;
    max_iters::Int = 3,
    tol::Float64 = 0.0,
    extended_precision::Bool = false,
)
    ferr = Inf
    iter = 0
    initial_ferr = NaN
    ω = NaN
    norm_u = norm(ws.u, InfNorm())
    ū = ws.work_n
    while iter <= max_iters
        iter += 1
        ferr_prev = ferr

        if extended_precision
            ferr = extended_precision_iterative_refinement_step!(ū, ws) / norm_u
        else
            ferr = fixed_precision_iterative_refinement_step!(ū, ws) / norm_u
        end

        if iter == 1
            initial_ferr = ferr
        end
        if iter == 2
            ω = 2 * ferr / ferr_prev^2
        end

        if 2 * ferr > ferr_prev
            # If the error is increasing, then we have probably reached the
            # machine precision limit. In this case, we return the previous
            # error estimate.
            ferr = ferr_prev
            break
        elseif ferr < tol
            ws.u .= ū
            break
        else
            ws.u .= ū
        end
    end

    # We don't need to scale ferr for the transformation from u to x
    # since the column scaling makes it such that we already compute
    # in the currect norm
    ws.ferr[] = ferr

    transform_u_to_x!(ws)

    return RefinementResult(ferr, initial_ferr, ω, iter)
end

"""
  fixed_precision_iterative_refinement_step!(ws::LinearSolveWorkspace

Perform one step of mixed precision iterative refinement.
"""
function fixed_precision_iterative_refinement_step!(u, ws::LinearSolveWorkspace)
    Δu = r = ws.work_n
    residual!(r, ws)
    lu_ldiv!(Δu, ws.LU, ws.ipiv, r)
    norm_Δu = norm(Δu, InfNorm())
    u .= ws.u .- Δu

    return norm_Δu
end

"""
  extended_precision_iterative_refinement_step!(ws::LinearSolveWorkspace)

Perform one step of mixed precision iterative refinement.
"""
function extended_precision_iterative_refinement_step!(u, ws::LinearSolveWorkspace)
    r̄ = ws.work_n_extended
    Δu = ws.work_n
    residual!(r̄, ws)
    lu_ldiv!(Δu, ws.LU, ws.ipiv, r̄)
    norm_Δu = norm(Δu, InfNorm())
    u .= ws.u .- Δu

    return norm_Δu
end

"""
    residual!(r, A, x, b)

Compute the residual `Ax-b` and store it in `r`
"""
function residual!(r::AbstractVector, ws::LinearSolveWorkspace)
    m, n = size(ws.A)
    r .= 0

    @inbounds for j = 1:n
        u_j = convert(eltype(r), ws.u[j])
        for i = 1:m
            r[i] += ws.A[i, j] * u_j
        end
    end

    @inbounds for i = 1:m
        r[i] -= ws.b[i]
    end

    r
end

"""
  gecon!(normtype, A, anorm)
Finds the reciprocal condition number of matrix `A`. If `normtype = I`,
the condition number is found in the infinity norm. If `normtype = O` or
`1`, the condition number is found in the one norm. `A` must be the
result of `getrf!` and `anorm` is the norm of `A` in the relevant norm.
"""
gecon!(
    normtype::AbstractChar,
    A::AbstractMatrix,
    anorm,
    work::AbstractVector,
    rwork::AbstractVector,
)


for (gecon, elty, relty) in
    ((:zgecon_, :ComplexF64, :Float64), (:cgecon_, :ComplexF32, :Float32))
    @eval begin
        #       SUBROUTINE ZGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK,
        #      $                   INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          NORM
        #       INTEGER            INFO, LDA, N
        #       DOUBLE PRECISION   ANORM, RCOND
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   RWORK( * )
        #       COMPLEX*16         A( LDA, * ), WORK( * )
        function gecon!(
            normtype::AbstractChar,
            A::AbstractMatrix{$elty},
            anorm::$relty,
            work::AbstractVector{$elty},
            rwork::AbstractVector{$relty},
        )
            LAPACK.chkstride1(A)
            n = LAPACK.checksquare(A)
            lda = max(1, stride(A, 2))
            rcond = Ref{$relty}()
            info = Ref{BlasInt}()
            ccall(
                (@blasfunc($gecon), libblastrampoline),
                Cvoid,
                (
                    Ref{UInt8},
                    Ref{BlasInt},
                    Ptr{$elty},
                    Ref{BlasInt},
                    Ref{$relty},
                    Ref{$relty},
                    Ptr{$elty},
                    Ptr{$relty},
                    Ptr{BlasInt},
                    Clong,
                ),
                normtype,
                n,
                A,
                lda,
                anorm,
                rcond,
                work,
                rwork,
                info,
                1,
            )
            LAPACK.chklapackerror(info[])
            rcond[]
        end


    end
end
