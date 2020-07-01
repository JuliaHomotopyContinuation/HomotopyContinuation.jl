# TODOs
#
# Modify certify(F,S) such that
#   F can be <: DP
#   S can be <: Result or MonodromyResult
#
# If there are zero certified solutions, this should also be displayed
#
# Write refine function (using Newton from HC.jl)
#
# Combine certify with the refine function
#
# Progress bar


export certify,
    CertifiedSolution,
    CertificationResult,
    ncertified,
    nreal_certified,
    ndistinct_certified,
    ndistinct_real_certified

Base.@kwdef struct SolutionCertificate
    initial_solution::Vector{ComplexF64}
    certified_solution::Union{Nothing,Vector{IComplexF64}} = nothing
    certified_extended_solution::Union{Nothing,Vector{IComplex{DoubleF64}}} = nothing
    is_real::Bool = false
    index::Union{Nothing,Int} = nothing
end

is_certified(C::SolutionCertificate) =
    !isnothing(C.certified_solution) || !isnothing(C.certified_extended_solution)
is_real(C::SolutionCertificate) = C.is_real

struct CertificationResult
    certified::Vector{SolutionCertificate}
    duplicates::Vector{Vector{Int}}
end
ncertified(R::CertificationResult) = count(is_certified, R.certified)
nreal_certified(R::CertificationResult) =
    count(r -> is_certified(r) && is_real(r), R.certified)
function ndistinct_certified(R::CertificationResult)
    ncert = ncertified(R)
    if isempty(R.duplicates)
        return ncert
    else
        return ncert - sum(length, R.duplicates) + length(R.duplicates)
    end
end
function ndistinct_real_certified(R::CertificationResult)
    ncert = nreal_certified(R)
    if isempty(R.duplicates)
        return ncert
    else
        ncert - sum(R.duplicates) do dup
            is_real(R.certified[dup[1]]) ? length(dup) - 1 : 0
        end
    end
end

function Base.show(io::IO, R::CertificationResult)
    println(io, "CertificationResult")
    println(io, "===================")
    println(io, "• $(length(R.certified)) solutions given")
    print(io, "• $(ncertified(R)) certified solutions")
    # if nreal_certified(R) > 0
    print(io, " ($(nreal_certified(R)) real)")
    # end
    println(io)
    print(io, "• $(ndistinct_certified(R)) distinct certified solutions")
    # if ndistinct_real_certified(R) > 0
    print(io, " ($(ndistinct_real_certified(R)) real)")
    # end
    println(io)
end

struct CertifyCache{T,M}
    jac_interpreter::ModelKit.Interpreter{T,2}
    newton_cache::NewtonCache{M}
    # norm
    norm::WeightedNorm{InfNorm}

    # data for krawczyc_step
    # complex case
    C::Matrix{ComplexF64}
    IJ::Matrix{IComplexF64}
    IJ_cache::ModelKit.InterpreterCache{IComplexF64}
    A::Matrix{IComplexF64}
    δx::Vector{IComplexF64}
    # complex extended precision case
    # TODO: Replace wih arb
    C_ext::Matrix{ComplexDF64}
    IJ_ext::Matrix{IComplex{DoubleF64}}
    IJ_ext_cache::ModelKit.InterpreterCache{IComplex{DoubleF64}}
end
function CertifyCache(F::AbstractSystem)
    m, n = size(F)
    jac_interpreter = ModelKit.jacobian_interpreter(System(F))
    IJ_cache = ModelKit.InterpreterCache(IComplexF64, jac_interpreter)
    CertifyCache(
        jac_interpreter,
        NewtonCache(F),
        WeightedNorm(InfNorm(), n),
        zeros(ComplexF64, m, n),
        zeros(IComplexF64, m, n),
        IJ_cache,
        zeros(IComplexF64, m, n),
        zeros(IComplexF64, m),
        zeros(ComplexDF64, m, n),
        zeros(IComplex{DoubleF64}, m, n),
        ModelKit.InterpreterCache(IComplex{DoubleF64}, jac_interpreter),
    )
end

"""
    certify(F::System, solutions; check_real, compile = $(COMPILE_DEFAULT[]))


"""
function certify(
    F::AbstractVector{Expression},
    X;
    parameters = Variable[],
    variables = setdiff(variables(F), parameters),
    variable_ordering = variables,
    variable_groups = nothing,
    kwargs...,
)
    sys = System(
        F,
        variables = variable_ordering,
        parameters = parameters,
        variable_groups = variable_groups,
    )
    certify(sys, X; kwargs...)
end
function certify(
    F::AbstractVector{<:MP.AbstractPolynomial},
    X;
    parameters = similar(MP.variables(F), 0),
    variables = setdiff(MP.variables(F), parameters),
    variable_ordering = variables,
    variable_groups = nothing,
    target_parameters = nothing,
    kwargs...,
)
    # as in solver_startsolutions():
    # handle special case that we have no parameters
    # to shift the coefficients of the polynomials to the parameters
    # this was the behaviour of HC.jl v1
    if isnothing(target_parameters) && isempty(parameters)
        sys, target_parameters = ModelKit.system_with_coefficents_as_params(
            F,
            variables = variable_ordering,
            variable_groups = variable_groups,
        )
    else
        sys = System(
            F,
            variables = variable_ordering,
            parameters = parameters,
            variable_groups = variable_groups,
        )
    end
    certify(sys, X; target_parameters = target_parameters, kwargs...)
end
function certify(F::System, args...; compile::Bool = COMPILE_DEFAULT[], kwargs...)
    certify(fixed(F; compile = compile), args...; kwargs...)
end
function certify(F::AbstractSystem, X, cache = CertifyCache(F); target_parameters = nothing, check_real = true, show_progress = true)
    _certify(F, X, target_parameters, cache; check_real = check_real, show_progress=show_progress)
end
function certify(
    F::AbstractSystem,
    X::MonodromyResult,
    cache = CertifyCache(F);
    target_parameters = parameters(X),
    check_real = true, show_progress = true
)
    _certify(F, solutions(X), target_parameters, cache; check_real = check_real, show_progress=show_progress)
end

function _certify(F::AbstractSystem, X::Result, p, cache::CertifyCache; check_real::Bool, show_progress::Bool)
    _certify(F, results(X; only_nonsingular = true), p, cache; check_real = check_real, show_progress=show_progress)
end
function _certify(F::AbstractSystem, X, p, cache::CertifyCache; check_real::Bool, show_progress::Bool)


    if show_progress
        n = length(X)
        desc = "Certifying $n solutions... "
        barlen = min(ProgressMeter.tty_width(desc), 40)
        progress = ProgressMeter.Progress(n; dt = 0.2, desc = desc, barlen = barlen, color=:green)
    end

    certified = map(enumerate(X)) do (i, x)
        if show_progress
            ProgressMeter.next!(progress)
        end
        _certify(F, x, p, cache, i; check_real = check_real)
    end
    duplicates = find_duplicates(certified)
    CertificationResult(certified, duplicates)
end
function _certify(
    F::AbstractSystem,
    r::PathResult,
    p,
    cache::CertifyCache,
    index = nothing;
    check_real::Bool,
)
    _certify(F, solution(r), p, cache, index; check_real = check_real)
end
function _certify(
    F::AbstractSystem,
    x₀::AbstractVector{<:Number},
    p::Union{Nothing,AbstractVector},
    cache::CertifyCache,
    index = nothing;
    check_real::Bool,
)
    @unpack C, IJ, IJ_cache, A, δx, norm, newton_cache = cache
    @unpack Δx, J, r = newton_cache

    @show typeof(F)

    init!(norm, x₀)
    res = newton(
        F,
        x₀,
        p,
        cache.norm,
        cache.newton_cache;
        atol = 4 * eps(),
        rtol = 4 * eps(),
        extended_precision = true,
    )
    #TODO: Don't assume res is of acc < 4eps()

    # The Newton method actually the accuracy of x + Δx, so use this, i.e., x̃ = x + Δx
    # This has the benefit that we don't need to compute F(x̃) or J(x̃) or an LU factorization
    x̃ = solution(res)
    x̃ .+= Δx

    # compute the inverse
    LA.inv!(C, J)
    # compute condition number
    κ = LA.opnorm(C, Inf) * LA.opnorm(J, Inf)
    # We have to make a first guess on the interval expansion factor with 1024, i.e.,
    # 1/(1024κ) is our guess for a sufficient ε.
    # But we want to have at least 8eps() and upper bound it by sqrt(eps())
    ε = clamp(inv(1024κ), 4 * eps(), 1e-12)
    # perform inflation
    δ = complex(Interval(-ε, ε), Interval(-ε, ε))
    # multiply by weights to account for relative accuracies
    x = x̃ .+ weights(norm) .* δ

    ModelKit.execute!(IJ, cache.jac_interpreter, x, p, IJ_cache)

    # Perform Krawczyk step: x′ = x̃ - Δx + (I - C * IJ) * (x - x̃)
    x′ = krawczyk_step(x̃, Δx, x, C, IJ, A, δx)
    certified = false
    if all2(isinterior, x′, x)
        certified = true
    else
        # # no success, check actual expansion factor
        K = maximum(diam, A) / (eps() * κ)
        # check if it makes sense to try again with new K
        ε = inv(16 * K * κ)
        if ε ≥ 4 * eps() # don't take eps() to have a buffer
            δ = complex(Interval(-ε, ε), Interval(-ε, ε))
            x .= x̃ .+ weights(norm) .* δ
            ModelKit.execute!(IJ, cache.jac_interpreter, x, p, IJ_cache)
            # Perform Krawczyk step: x′ = x̃ - Δx + (I - C * IJ) * (x - x̃)
            x′ = krawczyk_step(x̃, Δx, x, C, IJ, A, δx)

            certified = all2(isinterior, x′, x)
        else
            # TODO: extended precision necessary
            ε /= 16

            # 1) Perform newton steps to refine
            # 2) run krawczyc_step again
        end
    end

    if !certified
        return CertifiedSolution(initial_solution = x₀, index = index)
    end

    # We have a certified solution x.
    # Now we want to check if it is a real solution.
    if check_real
        is_real = all2((xᵢ′, xᵢ) -> isinterior(conj(xᵢ′), xᵢ), x′, x)
    else
        is_real = false
    end
    return SolutionCertificate(
        initial_solution = x₀,
        certified_solution = x,
        is_real = is_real,
        index = index,
    )
end

function krawczyk_step(x̃, Δx, x, C, IJ, A, δx)
    sqr_mul!(A, C, IJ)
    for i = 1:size(A, 1)
        A[i, i] -= 1
    end
    δx .= x .- x̃
    x′ = sqr_mul!(copy(δx), A, δx)
    @. x′ += x̃ - Δx
end

"Simple matrix multiplcation using the muladd method which is much faster four our case."
function sqr_mul!(C, A, B)
    n = size(A, 1)
    C .= zero(eltype(C))
    for j = 1:size(B, 2)
        for k = 1:n
            @inbounds bkj = B[k, j]
            for i = 1:n
                @inbounds C[i, j] = muladd(A[i, k], bkj, C[i, j])
            end
        end
    end
    C
end

function find_duplicates(R::AbstractVector{SolutionCertificate})
    # TODO: Do better than n^2 check
    duplicates = Dict{Int,Vector{Int}}()
    duplicated_indices = Set{Int}()
    for i = 1:length(R)
        is_certified(R[i]) && i ∉ duplicated_indices || continue
        rᵢ = R[i]
        for j = i+1:length(R)
            is_certified(R[j]) && j ∉ duplicated_indices || continue
            if all2(!isdisjoint, rᵢ.certified_solution, R[j].certified_solution)
                if haskey(duplicates, i)
                    push!(duplicated_indices, j)
                    push!(duplicates[i], j)
                else
                    push!(duplicated_indices, i)
                    push!(duplicated_indices, j)
                    duplicates[i] = [i, j]
                end
            end
        end
    end

    isempty(duplicates) ? Vector{Int}[] : collect(values(duplicates))
end
