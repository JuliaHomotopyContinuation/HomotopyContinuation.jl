using Arblib: Acb, AcbRef, AcbVector, AcbRefVector, AcbMatrix, AcbRefMatrix

export certify,
    SolutionCertificate,
    CertificationResult,
    CertificationCache,
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
    solution_candidate::AbstractVector
    certified::Bool
    real::Bool = false
    index::Union{Nothing,Int} = nothing
    prec::Int = 53
    # x₀, x₁ ∈ 𝕀ℂⁿ and are certified solution has x₁ ⊊ x₀
    x₀::Union{Nothing,AcbMatrix} = nothing
    x₁::Union{Nothing,AcbMatrix} = nothing
    # We also store a double precision representation of the midpoint of x₁
    # as the best available double precision estimate of the solution
    solution::Union{Nothing,Vector{ComplexF64}} = nothing
end


"""
    solution_candidate(C::SolutionCertificate)

Returns the given provided solution candidate.
"""
solution_candidate(C::SolutionCertificate) = C.solution_candidate
@deprecate initial_solution(cert::SolutionCertificate) solution_candidate(cert)

"""
    is_certified(C::SolutionCertificate)

Returns `true` if `C` certifies that the given initial solution corresponds to a true
solution.
"""
is_certified(C::SolutionCertificate) = C.certified

"""
    is_real(C::SolutionCertificate)

Returns `true` if `C` certifies that the given initial solution corresponds to a true
real solution of the system.
"""
is_real(C::SolutionCertificate) = C.real


"""
    is_positive(C::SolutionCertificate)

Returns `true` if `C` is certifiably a real, positive solution.
"""
function is_positive(C::SolutionCertificate)
    if !certified(C) || !is_real(C)
        return false
    else
        all(i -> Arblib.is_positive(real(Arblib.ref(x₀, i))), 1:length(C.x₀))
    end
end


"""
    certified_solution(C::SolutionCertificate)

Returns a vector of complex intervals where the true solution is contained in.
"""
certified_solution(C::SolutionCertificate) = C.x₀


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
    CertificationCache(F::AbstractSystem)

Contains the necessary data structures for [`certify`](@ref).
"""
Base.@kwdef mutable struct CertificationCache{T₁,T₂}
    eval_interpreter::ModelKit.Interpreter{T₁}
    jac_interpreter::ModelKit.Interpreter{T₂}
    newton_cache::NewtonCache{MatrixWorkspace{Matrix{ComplexF64}}}
    # norm
    norm::WeightedNorm{InfNorm}
    # data for krawczyc_step
    C::Matrix{ComplexF64}
    r₀::Vector{IComplexF64}
    Δx₀::Vector{IComplexF64}
    ix̃₀::Vector{IComplexF64}
    Jx₀::Matrix{IComplexF64}
    M::Matrix{IComplexF64}
    δx::Vector{IComplexF64}
    interpreter_cache::ModelKit.InterpreterCache{Vector{IComplexF64},Nothing}

    arb_prec::Int
    arb_C::AcbRefMatrix
    arb_r₀::AcbRefMatrix # m × 1
    arb_Δx₀::AcbRefMatrix # m × 1
    arb_x̃₀::AcbRefMatrix # m × 1
    arb_J_x₀::AcbRefMatrix
    arb_M::AcbRefMatrix
    arb_δx::AcbRefMatrix # m × 1
    arb_mag::Arblib.Mag
    arb_interpreter_cache::ModelKit.InterpreterCache{AcbRefVector,Acb}
end
function CertificationCache(F::AbstractSystem)
    m, n = size(F)
    m == n || error("Can only certify square polynomial systems.")
    f = System(F)
    interpreter = ModelKit.interpreter(f)
    jac_interpreter = ModelKit.jacobian_interpreter(f)
    tape_length = max(
        ModelKit.cache_min_length(interpreter),
        ModelKit.cache_min_length(jac_interpreter),
    )
    arb_prec = 128
    CertificationCache(;
        eval_interpreter = interpreter,
        jac_interpreter = jac_interpreter,
        newton_cache = NewtonCache(F; optimize_data_structure = false),
        norm = WeightedNorm(InfNorm(), m),
        C = zeros(ComplexF64, m, m),
        r₀ = zeros(IComplexF64, m),
        Δx₀ = zeros(IComplexF64, m),
        ix̃₀ = zeros(IComplexF64, m),
        Jx₀ = zeros(IComplexF64, m, m),
        M = zeros(IComplexF64, m, m),
        δx = zeros(IComplexF64, m),
        interpreter_cache = ModelKit.InterpreterCache(IComplexF64, tape_length),
        arb_prec = arb_prec,
        arb_C = AcbRefMatrix(m, m; prec = arb_prec),
        arb_r₀ = AcbRefMatrix(m, 1; prec = arb_prec),
        arb_Δx₀ = AcbRefMatrix(m, 1; prec = arb_prec),
        arb_x̃₀ = AcbRefMatrix(m, 1; prec = arb_prec),
        arb_J_x₀ = AcbRefMatrix(m, m; prec = arb_prec),
        arb_M = AcbRefMatrix(m, m; prec = arb_prec),
        arb_δx = AcbRefMatrix(m, 1; prec = arb_prec),
        arb_mag = Arblib.Mag(),
        arb_interpreter_cache = ModelKit.InterpreterCache(
            Acb,
            tape_length;
            prec = arb_prec,
        ),
    )
end

Base.setprecision(M::AcbRefMatrix, p::Int) = AcbRefMatrix(M.acb_mat, p)
function set_arb_precision!(cache::CertificationCache, p::Int)
    cache.arb_prec == p && return cache
    cache.arb_prec = p
    cache.arb_r₀ = setprecision(cache.arb_r₀, p)
    cache.arb_Δx₀ = setprecision(cache.arb_Δx₀, p)
    cache.arb_x̃₀ = setprecision(cache.arb_x̃₀, p)
    cache.arb_J_x₀ = setprecision(cache.arb_J_x₀, p)
    cache.arb_M = setprecision(cache.arb_M, p)
    cache.arb_δx = setprecision(cache.arb_δx, p)
    set_arb_precision!(arb_interpreter_cache)

    cache
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
function certify end


struct CertificationParameters
    params::Vector{ComplexF64}
    interval_params::Vector{IComplexF64}
    arb_interval_params::AcbRefVector
end
#
function CertificationParameters(p::AbstractVector; prec::Int = 128)
    arb_ip = AcbRefVector(length(p); prec = prec)
    for (i, p_i) in enumerate(p)
        x = arb_ip[i]
        x[] = p_i
    end
    CertificationParameters(
        convert(Vector{ComplexF64}, p),
        convert(Vector{IComplexF64}, p),
        arb_ip,
    )
end

certification_parameters(p::AbstractVector) = CertificationParameters(p)
certification_parameters(::Nothing) = nothing


complexF64_params(C::CertificationParameters) = C.params
complexF64_interval_params(C::CertificationParameters) = C.interval_params
arb_interval_params(C::CertificationParameters) = C.arb_interval_params
complexF64_params(::Nothing) = nothing
complexF64_interval_params(::Nothing) = nothing
arb_interval_params(::Nothing) = nothing

function certify(
    F::AbstractSystem,
    solution_candidates::AbstractArray{<:AbstractArray{<:Number}},
    p::Union{Nothing,CertificationParameters} = nothing,
    cache::CertificationCache = CertificationCache(F);
    show_progress::Bool = true,
    max_precision::Int = 256,
)
    certificates = SolutionCertificate[]

    if !show_progress
        for (i, s) in enumerate(solution_candidates)
            push!(certificates, certify_solution(F, x, s, cache, i))
        end
    else
        # Create progress meter
        n = length(solution_candidates)
        desc = "Certifying $n solutions... "
        barlen = min(ProgressMeter.tty_width(desc, stdout), 40)
        progress = ProgressMeter.Progress(
            n;
            dt = 0.2,
            desc = desc,
            barlen = barlen,
            color = :green,
            output = stdout,
        )
        progress.tlast += progress.dt

        # Go over all solution
        ncertified = 0
        nreal_certified = 0
        for (i, s) in enumerate(solution_candidates)
            r = certify_solution(F, s, p, cache, i)
            push!(certificates, r)
            ncertified += is_certified(r)
            nreal_certified += is_certified(r) && is_real(r)
            update_certify_progress!(progress, i, ncertified, nreal_certified)
        end
    end
    # duplicates = find_duplicates_fast(certificates; show_progress = show_progress)
    @warn "Duplicate check disabled"
    duplicates = Vector{Vector{Int}}()
    CertificationResult(certificates, duplicates)
end

function update_certify_progress!(progress, k, ncertified, nreal_certified)
    t = time()
    if k == progress.n || t > progress.tlast + progress.dt
        showvalues = make_certify_showvalues(k, ncertified, nreal_certified)
        ProgressMeter.update!(progress, k; showvalues = showvalues)
    end
    nothing
end
@noinline function make_certify_showvalues(k, ncertified, nreal_certified)
    (
        ("# solutions considered", k),
        ("# certified solutions (real)", "$(ncertified) ($nreal_certified)"),
    )
end


function certify_solution(
    F::AbstractSystem,
    solution_candidate::Vector{ComplexF64},
    cert_params::Union{Nothing,CertificationParameters},
    cert_cache::CertificationCache,
    index::Int,
    max_precision::Int = 256,
    refine_solution::Bool = true,
)
    @unpack C, arb_C, arb_x̃₀ = cert_cache

    # refine solution to machine precicision
    init!(cert_cache.norm, solution_candidate)
    res = newton(
        F,
        solution_candidate,
        complexF64_params(cert_params),
        cert_cache.norm,
        cert_cache.newton_cache;
        # we have a weighted norm, so atol is fine
        atol = 8 * eps(),
        # this should already be an approximate zero
        contraction_factor = 0.5,
        extended_precision = true,
        max_iters = 8,
    )

    # Abort if we cannot refine to this accuracy
    if !is_success(res)
        return SolutionCertificate(
            solution_candidate = solution_candidate,
            certified = false,
            index = index,
        )
    end

    x̃₀ = solution(res)
    LA.inv!(C, cert_cache.newton_cache.J)

    certified, x₁, x₀, is_real = ε_inflation_krawczyk(x̃₀, cert_params, C, cert_cache)

    if certified
        return SolutionCertificate(
            solution_candidate = solution_candidate,
            certified = true,
            real = is_real,
            index = index,
            prec = 53,
            x₀ = AcbMatrix(x₀; prec = 53),
            x₁ = AcbMatrix(x₁; prec = 53),
            solution = mid.(x₁),
        )
    end

    # If not certified we do the computation in higher precision using Arb
    prec = 128
    set_arb_precision!(cert_cache, prec)

    # We keep the same C matrix for now.
    for j = 1:size(C, 2), i = 1:size(C, 1)
        arb_C[i, j][] = C[i, j]
    end
    for (i, x̃₀_i) in enumerate(x̃₀)
        arb_x̃₀[i][] = x̃₀_i
    end
    while (prec <= max_precision)
        certified, arb_x₁, arb_x₀, is_real =
            arb_ε_inflation_krawczyk(arb_x̃₀, cert_params, arb_C, cert_cache; prec = prec)

        if certified
            return SolutionCertificate(
                solution_candidate = solution_candidate,
                certified = true,
                real = is_real,
                index = index,
                prec = prec,
                x₀ = arb_x₀,
                x₁ = arb_x₁,
                solution = [ComplexF64(arb_x₁[i]) for i = 1:size(C, 1)],
            )
        end

        prec += 64
    end

    # certification failed
    return SolutionCertificate(solution_candidate = solution_candidate, index = index)
end



function Base.Vector{ComplexF64}(A::Arblib.AcbMatrixLike)
    @assert size(A, 2) == 1
    [ComplexF64(Arblib.ref(A, i, 1).acb_ptr) for i = 1:size(A, 1)]
end

function ε_inflation_krawczyk(x̃₀, p::Union{Nothing,CertificationParameters}, C, cert_cache)
    @unpack C, r₀, Δx₀, ix̃₀, Jx₀, M, δx, interpreter_cache = cert_cache
    J_x₀ = Jx₀

    ix̃₀ .= IComplexF64.(x̃₀)
    ModelKit.execute!(
        r₀,
        cert_cache.eval_interpreter,
        ix̃₀,
        complexF64_interval_params(p),
        interpreter_cache,
    )
    # iΔx = C * r₀
    sqr_mul!(Δx₀, C, r₀)

    # Perform ε-inflation. We choose a different εᵢ per coordinate. This matches our strategy
    # to use a weighted norm.
    x₀ = map(x̃₀, Δx₀) do x̃₀_i, Δx₀_i
        Δᵢ = mag(Δx₀_i)
        εᵢ = max(4Δᵢ, 8eps() * abs(x̃₀_i))
        complex(
            Interval(real(x̃₀_i) - εᵢ, real(x̃₀_i) + εᵢ),
            Interval(imag(x̃₀_i) - εᵢ, imag(x̃₀_i) + εᵢ),
        )
    end

    ModelKit.execute!(
        nothing,
        J_x₀,
        cert_cache.jac_interpreter,
        x₀,
        complexF64_interval_params(p),
        interpreter_cache,
    )

    x₁ = similar(x₀)
    # x₁ = (x̃₀ - C * F(x₀)) + (I - C * J(x₀)) * (x₀ - x̃₀)
    #    = (x̃₀ - C * F(x₀)) - (C * J(x₀) - I) * (x₀ - x̃₀)
    #    = (x̃₀ - Δx₀) - (C * J(x₀) - I) * (x₀ - x̃₀)

    # Define M = C * J(x₀) - I
    sqr_mul!(M, C, J_x₀)
    for i = 1:size(M, 1)
        M[i, i] -= 1
    end
    # Necessary condition is ||M|| < 1 / √2
    # We lower bound 1 / √2 by 0.7071
    if IntervalArithmetic.inf_norm_bound(M) < 0.7071
        δx .= x₀ .- x̃₀
        sqr_mul!(x₁, M, δx)
        x₁ .= (x̃₀ .- Δx₀) .- x₁
        certified = all2(isinterior, x₁, x₀)
        is_real =
            certified ? all2((x₁_i, x₀_i) -> isinterior(conj(x₁_i), x₀_i), x₁, x₀) : false
    else
        certified = false
        is_real = false
    end
    certified, x₁, x₀, is_real
end

function arb_ε_inflation_krawczyk(
    x̃₀::AcbRefMatrix,
    p::Union{Nothing,CertificationParameters},
    C::Arblib.AcbMatrixLike,
    cert_cache;
    prec::Int,
)
    @unpack arb_r₀, arb_Δx₀, arb_x̃₀, arb_J_x₀, arb_M, arb_δx, arb_mag = cert_cache
    r₀, Δx₀, x̃₀, J_x₀, M, δx, m =
        arb_r₀, arb_Δx₀, arb_x̃₀, arb_J_x₀, arb_M, arb_δx, arb_mag
    @unpack arb_interpreter_cache = cert_cache

    m = arb_mag
    x̃₀′ = x̄₁ = Δx₀
    # Perform a simple newton refinement using arb_C until we cannot improve the
    # accuracy anymore
    acc = max_iters = 10
    for i = 1:max_iters
        ModelKit.execute!(
            arb_r₀,
            cert_cache.eval_interpreter,
            x̃₀,
            arb_interval_params(p),
            arb_interpreter_cache,
        )
        Arblib.get_mid!(r₀, r₀)
        Arblib.mul!(Δx₀, C, r₀)

        acc_new = Arblib.get(Arblib.bound_inf_norm!(m, Δx₀))
        if acc_new > acc || i == max_iters
            break
        end
        acc = acc_new
        Arblib.sub!(x̃₀, x̃₀, Δx₀)
        Arblib.get_mid!(x̃₀, x̃₀)
    end

    # Perform ε-inflation
    n = length(x̃₀)
    x₀ = AcbRefMatrix(n, 1; prec = prec)
    Arblib.get_mid!(x₀, x̃₀)
    for i = 1:n
        m[] = 128magF64(Δx₀[i], m)
        Arblib.add_error!(x₀[i], m)
    end

    ModelKit.execute!(
        nothing,
        J_x₀,
        cert_cache.jac_interpreter,
        x₀,
        arb_interval_params(p),
        arb_interpreter_cache,
    )

    # x₁ = (x̃₀ - C * F(x₀)) + (I - C * J(x₀)) * (x₀ - x̃₀)
    #    = (x̃₀ - C * F(x₀)) - (C * J(x₀) - I) * (x₀ - x̃₀)
    #    = (x̃₀ - Δx₀) - (C * J(x₀) - I) * (x₀ - x̃₀)

    # Define M = C * J(x₀) - I
    Arblib.mul!(M, C, J_x₀)
    for i = 1:n
        Arblib.sub!(M[i, i], M[i, i], 1)
    end

    x₁ = similar(x₀)

    # Necessary condition is ||M|| < 1 / √2
    # We lower bound 1 / √2 by 0.7071
    if Arblib.get(Arblib.bound_inf_norm!(m, M)) < 0.7071
        # x₀ - x̃₀
        Arblib.sub!(δx, x₀, x̃₀)
        # (C * J(x₀) - I) * (x₀ - x̃₀)
        Arblib.mul!(x₁, M, δx)
        # (x̃₀ - Δx₀)
        Arblib.sub!(x̃₀′, x̃₀, Δx₀)
        # (x̃₀ - Δx₀) - (C * J(x₀) - I) * (x₀ - x̃₀)
        Arblib.sub!(x₁, x̃₀′, x₁)

        certified = Bool(Arblib.contains(x₀, x₁))

        # check realness
        Arblib.conjugate!(x̄₁, x₁)
        is_real = Bool(Arblib.contains(x₀, x̄₁))
    else
        certified = false
        is_real = false
    end

    if !certified
        # Update the approximation of the inverse
        Arblib.get_mid!(J_x₀, J_x₀)
        Arblib.inv!(C, J_x₀)
    end

    certified, AcbMatrix(x₁), AcbMatrix(x₀), is_real
end

magF64(x, mag) = Arblib.get(Arblib.get!(mag, x))

"Simple matrix multiplication using the muladd method which is much faster four our case."
function sqr_mul!(C, A, B)
    n = size(A, 1)
    C .= zero(eltype(C))
    @inbounds for j = 1:size(B, 2)
        for k = 1:n
            bkj = B[k, j]
            for i = 1:n
                C[i, j] = muladd(A[i, k], bkj, C[i, j])
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
            if Bool(Arblib.overlaps(rᵢ.certified_solution, R[j].certified_solution))
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

radF64(x::Arblib.AcbLike) = max(
    Float64(Arblib.radref(Arblib.realref(x))),
    Float64(Arblib.radref(Arblib.imagref(x))),
)

function find_duplicates_fast(
    R::AbstractVector{SolutionCertificate};
    show_progress::Bool = true,
)
    # TODO: Do better than n^2 check
    duplicates = Dict{Int,Vector{Int}}()
    duplicated_indices = Set{Int}()
    tree =
        VoronoiTree{Float64,Int}(2 * length(first(R).initial_solution); metric = InfNorm())
    t0 = time()
    displayed_info = !show_progress
    max_norm = 0.0
    max_diam = 0.0
    inf_norm = InfNorm()
    m = Arblib.Mag()

    add_duplicate = (i, j) -> begin
        if haskey(duplicates, i)
            push!(duplicated_indices, j)
            push!(duplicates[i], j)
        else
            push!(duplicated_indices, i)
            push!(duplicated_indices, j)
            duplicates[i] = [i, j]
        end
    end

    for i = 1:length(R)
        is_certified(R[i]) && i ∉ duplicated_indices || continue
        rᵢ = R[i]
        if !displayed_info && (time() - t0 > 0.2)
            @info "Checking for duplicates"
            displayed_info = true
        end
        mid = reinterpret(Float64, rᵢ.certified_midpoint::Vector{ComplexF64})
        norm_mid = inf_norm(mid)
        max_norm = max(max_norm, norm_mid)
        rad = maximum(radF64, rᵢ.certified_solution)
        max_diam = max(max_diam, 2 * rad) # this is exact
        tol = nextfloat(max_diam + nextfloat(eps() * max_norm))

        (id, added) = add!(tree, mid, i, tol)
        if !added
            if Bool(Arblib.overlaps(rᵢ.certified_solution, R[id].certified_solution))
                add_duplicate(i, id)
            else
                for j = 1:i-1
                    is_certified(R[j]) && j ∉ duplicated_indices || continue
                    if Bool(Arblib.overlaps(rᵢ.certified_solution, R[j].certified_solution))
                        add_duplicate(i, j)
                    end
                end
            end
        end

    end

    isempty(duplicates) ? Vector{Int}[] : collect(values(duplicates))
end


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


function certify(F::System, args...; compile::Union{Bool,Symbol} = false, kwargs...)
    certify(fixed(F; compile = compile), args...; kwargs...)
end

function certify(
    F::AbstractSystem,
    X::MonodromyResult,
    cache::CertificationCache = CertificationCache(F);
    kwargs...,
)
    certify(F, solutions(X), parameters(X), cache; kwargs...)
end

function certify(
    F::AbstractSystem,
    X,
    p::AbstractArray,
    cache::CertificationCache = CertificationCache(F);
    target_parameters = nothing,
    kwargs...,
)
    if !isnothing(p)
        target_parameters = p
    end

    cert_params = certification_parameters(target_parameters)
    certify(F, X, cert_params, cache; kwargs...)
end

function certify(
    F::AbstractSystem,
    r::PathResult,
    p::Union{Nothing,CertificationParameters} = nothing,
    cache::CertificationCache = CertificationCache(F);
    kwargs...,
)
    certify(F, [solution(r)], p, cache, index; kwargs...)
end


function certify(
    F::AbstractSystem,
    r::AbstractVector{PathResult},
    p::Union{Nothing,CertificationParameters} = nothing,
    cache::CertificationCache = CertificationCache(F);
    kwargs...,
)
    certify(F, solution.(r), p, cache, index; kwargs...)
end

function certify(
    F::AbstractSystem,
    x::AbstractVector{<:Number},
    p::Union{Nothing,CertificationParameters} = nothing,
    cache::CertificationCache = CertificationCache(F);
    kwargs...,
)
    certify(F, [x], p, cache, index; kwargs...)
end

function certify(
    F::AbstractSystem,
    X::Result,
    p::Union{Nothing,CertificationParameters} = nothing,
    cache::CertificationCache = CertificationCache(F);
    kwargs...,
)
    certify(F, solutions(X; only_nonsingular = true), p, cache; kwargs...)
end
