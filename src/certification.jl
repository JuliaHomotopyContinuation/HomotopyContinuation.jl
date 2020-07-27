export certify,
    SolutionCertificate,
    CertificationResult,
    CertifyCache,
    is_certified,
    is_real,
    is_positive,
    certified_solution,
    initial_solution,
    certificates,
    ncertified,
    nreal_certified,
    ndistinct_certified,
    ndistinct_real_certified,
    save

"""
    SolutionCertificate

Result of [`certify`](@ref) for a single solution. Contains the initial solutions
and if the certification was successfull a vector of complex intervals where the true
solution is contained in.
"""
Base.@kwdef struct SolutionCertificate
    initial_solution::Vector{ComplexF64}
    certified_solution::Union{Nothing,Vector{IComplexF64}} = nothing
    is_real::Bool = false
    index::Union{Nothing,Int} = nothing
end

"""
    initial_solution(C::SolutionCertificate)

Returns the given initial solution.
"""
initial_solution(C::SolutionCertificate) = C.initial_solution

"""
    is_certified(C::SolutionCertificate)

Returns `true` if `C` certifies that the given initial solution corresponds to a true
solution.
"""
is_certified(C::SolutionCertificate) = !isnothing(C.certified_solution)

"""
    is_real(C::SolutionCertificate)

Returns `true` if `C` certifies that the given initial solution corresponds to a true
real solution of the system.
"""
is_real(C::SolutionCertificate) = C.is_real


"""
    is_positive(C::SolutionCertificate)

Returns `true` if `C` is certifiably a real, positive solution.
"""
function is_positive(C::SolutionCertificate)
    if isnothing(C.certified_solution) || !is_real(C)
        return false
    else
        all(zᵢ -> real(zᵢ).lo > 0, C.certified_solution)
    end
end


"""
    certified_solution(C::SolutionCertificate)

Returns a vector of complex intervals where the true solution is contained in.
"""
certified_solution(C::SolutionCertificate) = C.certified_solution


function Base.show(f::IO, cert::SolutionCertificate)
    if !isnothing(cert.index)
        println(f, "index = ", cert.index)
    end
    println(f, "initial_solution = [")
    for z in initial_solution(cert)
        print(f, "  ")
        print(f, z)
        println(f, ",")
    end
    println(f, "]")
    println(f, "certified_solution = [")
    for z in certified_solution(cert)
        print(f, "  ")
        print(f, z)
        println(f, ",")
    end
    println(f, "]")
    println(f, "is_real = ", is_real(cert))
end

"""
    CertificationResult

The result of [`certify`](@ref) for multiple solutions.
Contains a vector of [`SolutionCertificate`](@ref) as well as a list of certificates
which correspond to the same true solution.
"""
struct CertificationResult
    certificates::Vector{SolutionCertificate}
    duplicates::Vector{Vector{Int}}
end

"""
    certificates(R::CertificationResult)

Obtain the stored [`SolutionCertificate`](@ref)s.
"""
certificates(R::CertificationResult) = R.certificates

"""
    ncertified(R::CertificationResult)

Returns the number of certified solutions.
"""
ncertified(R::CertificationResult) = count(is_certified, R.certificates)

"""
    nreal_certified(R::CertificationResult)

Returns the number of certified real solutions.
"""
nreal_certified(R::CertificationResult) =
    count(r -> is_certified(r) && is_real(r), R.certificates)

"""
    ndistinct_certified(R::CertificationResult)

Returns the number of distinct certified solutions.
"""
function ndistinct_certified(R::CertificationResult)
    ncert = ncertified(R)
    if isempty(R.duplicates)
        return ncert
    else
        return ncert - sum(length, R.duplicates) + length(R.duplicates)
    end
end

"""
    ndistinct_real_certified(R::CertificationResult)

Returns the number of distinct certified real solutions.
"""
function ndistinct_real_certified(R::CertificationResult)
    ncert = nreal_certified(R)
    if isempty(R.duplicates)
        return ncert
    else
        ncert - sum(R.duplicates) do dup
            is_real(R.certificates[dup[1]]) ? length(dup) - 1 : 0
        end
    end
end

function Base.show(io::IO, R::CertificationResult)
    println(io, "CertificationResult")
    println(io, "===================")
    println(io, "• $(length(R.certificates)) solutions given")
    print(io, "• $(ncertified(R)) certified solutions")
    print(io, " ($(nreal_certified(R)) real)")
    println(io)
    print(io, "• $(ndistinct_certified(R)) distinct certified solutions")
    print(io, " ($(ndistinct_real_certified(R)) real)")
end

"""
    save(filename, C::CertificationResult)

Store a text representation of the certification result `C` on disk.
"""
function save(filename, R::CertificationResult)
    open(filename, "w") do f
        println(f, "## Summary")
        show(f, R)
        println(f, "\n\n## Certificates")
        for cert in certificates(R)
            show(f, cert)
            println(f)
        end
    end
    filename
end
"""
    CertifyCache(F::AbstractSystem)

Contains the necessary data structures for [`certify`](@ref).
"""
struct CertifyCache{T,M}
    eval_interpreter::ModelKit.Interpreter{T,1}
    jac_interpreter::ModelKit.Interpreter{T,2}
    newton_cache::NewtonCache{M}
    # norm
    norm::WeightedNorm{InfNorm}
    # data for krawczyc_step
    ir::Vector{IComplexF64}
    iΔx::Vector{IComplexF64}
    ix̃::Vector{IComplexF64}
    dix̃::Vector{IComplex{DoubleF64}}
    C::Matrix{ComplexF64}
    IJ::Matrix{IComplexF64}
    IJ_cache::ModelKit.InterpreterCache{IComplexF64}
    dIJ_cache::ModelKit.InterpreterCache{IComplex{DoubleF64}}
    A::Matrix{IComplexF64}
    δx::Vector{IComplexF64}
    # TODO: complex extended precision case with arb
end
function CertifyCache(F::AbstractSystem)
    m, n = size(F)
    m == n || error("Can only certify square polynomial systems.")
    eval_interpreter, jac_interpreter =
        ModelKit.promote_common_constants(ModelKit.evaluate_jacobian_interpreter(System(
            F,
        ))...)
    IJ_cache = ModelKit.InterpreterCache(Vector{IComplexF64}(
        undef,
        max(length(eval_interpreter.instructions), length(jac_interpreter.instructions)),
    ))
    dIJ_cache = ModelKit.InterpreterCache(Vector{IComplex{DoubleF64}}(
        undef,
        max(length(eval_interpreter.instructions), length(jac_interpreter.instructions)),
    ))
    CertifyCache(
        eval_interpreter,
        jac_interpreter,
        NewtonCache(F),
        WeightedNorm(InfNorm(), n),
        zeros(IComplexF64, m),
        zeros(IComplexF64, m),
        zeros(IComplexF64, m),
        zeros(IComplex{DoubleF64}, m),
        zeros(ComplexF64, m, n),
        zeros(IComplexF64, m, n),
        IJ_cache,
        dIJ_cache,
        zeros(IComplexF64, m, n),
        zeros(IComplexF64, m),
    )
end

"""
    certify(F, solutions, [p, certify_cache]; options...)
    certify(F, result, [p, certify_cache]; options...)

Attempt to certify that the given approximate `solutions` correspond to true solutions
of the polynomial system ``F(x;p)``. The system ``F`` has to be an (affine) square
polynomial system. Also attemps to certify for each solutions whether it
approximates a real solution. The certification is done using interval arithmetic and the
Krawczyk method[^Moo77]. Returns a [`CertificationResult`](@ref) which additionall returns
the number of distinct solutions. For more details of the implementation see [^BRT20].

## Options

* `show_progress = true`: If `true` shows a progress bar of the certification process.
* `compile = $(COMPILE_DEFAULT[])`: See the [`solve`](@ref) documentation.


## Example
We take the [first example](https://www.juliahomotopycontinuation.org/guides/introduction/#a-first-example) from our
introduction guide.
```julia
@var x y
# define the polynomials
f₁ = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
f₂ = x^2+2x*y^2 - 2y^2 - 1/2
F = System([f₁, f₂], variables = [x,y])
result = solve(F)
```
```
Result with 18 solutions
========================
• 18 paths tracked
• 18 non-singular solutions (4 real)
• random seed: 0xcaa483cd
• start_system: :polyhedral
```
We see that we obtain 18 solutions and it seems that 4 solutions are real. However,
this is based on heuristics. To be absolute certain we can certify the result

```julia
certify(F, result)
```
```
CertificationResult
===================
• 18 solutions given
• 18 certified solutions (4 real)
• 18 distinct certified solutions (4 real)
```

and see that there are indeed 18 solutions and that they are all distinct.


[^Moo77]: Moore, Ramon E. "A test for existence of solutions to nonlinear systems." SIAM Journal on Numerical Analysis 14.4 (1977): 611-615.
[^BRT20]: Breiding, P., Rose, K. and Timme, S. "Certifying roots of polynomial systems using interval arithmetic." In preparation (2020).
"""
function certify(
    F::AbstractVector{Expression},
    X,
    p = nothing;
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
    certify(sys, X, p; kwargs...)
end
function certify(
    F::AbstractVector{<:MP.AbstractPolynomial},
    X,
    p = nothing;
    parameters = similar(MP.variables(F), 0),
    variables = setdiff(MP.variables(F), parameters),
    variable_ordering = variables,
    variable_groups = nothing,
    target_parameters = nothing,
    kwargs...,
)
    if !isnothing(p)
        target_parameters = p
    end
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
    certify(sys, X, p; target_parameters = target_parameters, kwargs...)
end
function certify(F::System, args...; compile::Union{Bool,Symbol} = COMPILE_DEFAULT[], kwargs...)
    certify(fixed(F; compile = compile), args...; kwargs...)
end
function certify(
    F::AbstractSystem,
    X::MonodromyResult,
    cache::CertifyCache = CertifyCache(F);
    show_progress = true,
)
    certify(F, solutions(X), parameters(X), cache; show_progress = show_progress)
end
function certify(
    F::AbstractSystem,
    X,
    p = nothing,
    cache::CertifyCache = CertifyCache(F);
    target_parameters = nothing,
    show_progress = true,
)
    if !isnothing(p)
        target_parameters = p
    end
    if !isnothing(target_parameters)
        # convert target parameters to something certifyable
        if isa(target_parameters, AbstractVector{<:Rational}) ||
           isa(target_parameters, AbstractVector{<:Complex{<:Rational}})
            ip = target_parameters
            dip = target_parameters
        elseif isa(target_parameters, AbstractVector{<:Complex})
            ip = IComplexF64.(target_parameters)
            dip = IComplex{DoubleF64}.(target_parameters)
        else
            ip = Interval{Float64}.(target_parameters)
            dip = Interval{DoubleF64}.(target_parameters)
        end
    else
        ip = dip = nothing
    end

    if X isa PathResult || X isa AbstractVector{<:Number}
        X = [X]
    end
    _certify(F, X, target_parameters, ip, dip, cache; show_progress = show_progress)
end

function _certify(
    F::AbstractSystem,
    X::Result,
    p,
    ip,
    dip,
    cache::CertifyCache;
    show_progress::Bool,
)
    _certify(
        F,
        results(X; only_nonsingular = true),
        p,
        ip,
        dip,
        cache;
        show_progress = show_progress,
    )
end
function _certify(
    F::AbstractSystem,
    X,
    p,
    ip,
    dip,
    cache::CertifyCache;
    show_progress::Bool,
)
    certificates = SolutionCertificate[]
    if show_progress
        n = length(X)
        desc = "Certifying $n solutions... "
        barlen = min(ProgressMeter.tty_width(desc), 40)
        progress = ProgressMeter.Progress(
            n;
            dt = 0.2,
            desc = desc,
            barlen = barlen,
            color = :green,
            output = stdout,
        )
        progress.tlast += progress.dt
        ncertified = 0
        nreal_certified = 0
        for (i, x) in enumerate(X)
            r = _certify(F, x, p, ip, dip, cache, i)
            push!(certificates, r)
            ncertified += is_certified(r)
            nreal_certified += is_certified(r) && is_real(r)
            if show_progress
                update_certify_progress!(progress, i, ncertified, nreal_certified)
            end
        end
    else
        for (i, x) in enumerate(X)
            push!(certificates, _certify(F, x, p, ip, dip, cache, i))
        end
    end
    duplicates = find_duplicates(certificates; show_progress = show_progress)
    CertificationResult(certificates, duplicates)
end

function update_certify_progress!(progress, k, ncertified, nreal_certified)
    t = time()
    if k == progress.n || t > progress.tlast + progress.dt
        showvalues = make_certifiy_showvalues(k, ncertified, nreal_certified)
        ProgressMeter.update!(progress, k; showvalues = showvalues)
    end
    nothing
end
@noinline function make_certifiy_showvalues(k, ncertified, nreal_certified)
    (
        ("# solutions considered", k),
        ("# certified solutions (real)", "$(ncertified) ($nreal_certified)"),
    )
end

function _certify(
    F::AbstractSystem,
    r::PathResult,
    p,
    ip,
    dip,
    cache::CertifyCache,
    index = nothing,
)
    _certify(F, solution(r), p, ip, dip, cache, index)
end
function _certify(
    F::AbstractSystem,
    x₀::AbstractVector{<:Number},
    p::Union{Nothing,AbstractVector},
    ip,
    dip,
    cache::CertifyCache,
    index = nothing,
)
    _certify(F, convert(Vector{ComplexF64}, x₀), p, ip, dip, cache, index)
end
function _certify(
    F::AbstractSystem,
    x₀::Vector{ComplexF64},
    p::Union{Nothing,AbstractVector},
    ip,
    dip,
    cache::CertifyCache,
    index = nothing,
)
    @unpack C, IJ, IJ_cache, dIJ_cache, A, δx, ir, iΔx, ix̃, dix̃, norm, newton_cache =
        cache
    @unpack Δx, J, r = newton_cache

    init!(norm, x₀)
    res = newton(
        F,
        x₀,
        p,
        norm,
        newton_cache;
        atol = 4 * eps(),
        rtol = 4 * eps(),
        extended_precision = true,
        max_iters = 4,
    )
    x̃ = solution(res)
    LA.inv!(C, J)

    ix̃ .= IComplexF64.(x̃)
    ModelKit.execute!(ir, cache.eval_interpreter, ix̃, ip, IJ_cache)
    sqr_mul!(iΔx, C, ir)

    x = map(x̃, iΔx) do x̃ᵢ, iΔxᵢ
        Δᵢ = abs(mid(iΔxᵢ)) + diam(iΔxᵢ)
        εᵢ = max(4Δᵢ, 8eps() * abs(x̃ᵢ))
        complex(
            Interval(real(x̃ᵢ) - εᵢ, real(x̃ᵢ) + εᵢ),
            Interval(imag(x̃ᵢ) - εᵢ, imag(x̃ᵢ) + εᵢ),
        )
    end
    ModelKit.execute!(IJ, cache.jac_interpreter, x, ip, IJ_cache)
    x′ = krawczyk_step(x̃, iΔx, x, C, IJ, A, δx)
    certified = maximum(mag, A) < 0.7 && all2(isinterior, x′, x)

    if !certified
        dix̃ .= IComplex{DoubleF64}.(x̃)
        ModelKit.execute!(ir, cache.eval_interpreter, dix̃, dip, dIJ_cache)
        # multiplication in standard precision
        sqr_mul!(iΔx, C, ir)
        # ε-inflation
        x = map(x̃, iΔx) do x̃ᵢ, iΔxᵢ
            Δᵢ = mag(iΔxᵢ)
            εᵢ = max(8Δᵢ, 8eps() * abs(x̃ᵢ))
            a, b = reim(ComplexDF64(x̃ᵢ))
            complex(Interval(a - εᵢ, a + εᵢ), Interval(b - εᵢ, b + εᵢ))
        end

        # Compute interval jacobian in extended precision
        ModelKit.execute!(IJ, cache.jac_interpreter, x, dip, dIJ_cache)
        x′ = krawczyk_step(x̃, iΔx, x, C, IJ, A, δx)
        certified = maximum(mag, A) < 0.7 && all2(isinterior, x′, x)
    end

    if !certified
        return SolutionCertificate(initial_solution = x₀, index = index)
    end

    # We have a certified solution x.
    # Now we want to check if it is a real solution.
    is_real = certified && all2((xᵢ′, xᵢ) -> isinterior(conj(xᵢ′), xᵢ), x′, x)

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

"Simple matrix multiplication using the muladd method which is much faster four our case."
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

function find_duplicates(R::AbstractVector{SolutionCertificate}; show_progress::Bool = true)
    # TODO: Do better than n^2 check
    duplicates = Dict{Int,Vector{Int}}()
    duplicated_indices = Set{Int}()
    t0 = time()
    displayed_info = !show_progress
    for i = 1:length(R)
        is_certified(R[i]) && i ∉ duplicated_indices || continue
        rᵢ = R[i]
        if !displayed_info && (time() - t0 > 0.2)
            @info "Checking for duplicates"
            displayed_info = true
        end
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
