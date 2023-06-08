using Arblib: Acb, AcbRef, AcbVector, AcbRefVector, AcbMatrix, AcbRefMatrix

export certify,
    SolutionCertificate,
    CertificationResult,
    CertificationCache,
    is_certified,
    is_real,
    is_complex,
    is_positive,
    solution_candidate,
    certified_solution_interval,
    certified_solution_interval_after_krawczyk,
    certificate_index,
    solution_approximation,
    certificates,
    distinct_certificates,
    distinct_solutions,
    ncertified,
    nreal_certified,
    ncomplex_certified,
    ndistinct_certified,
    ndistinct_real_certified,
    ndistinct_complex_certified,
    show_straight_line_program,
    save,
    DistinctCertifiedSolutions,
    add_solution!,
    solutions,
    certificates,
    distinct_certified_solutions,
    # deprecated
    initial_solution,
    certified_solution


abstract type AbstractSolutionCertificate end


"""
    SolutionCertificate

Result of [`certify`](@ref) for a single solution. Contains the initial solutions
and if the certification was successfull a vector of complex intervals where the true
solution is contained in.
"""
Base.@kwdef struct SolutionCertificate <: AbstractSolutionCertificate
    solution_candidate::AbstractVector
    certified::Bool
    real::Bool = false
    complex::Bool = false
    index::Union{Nothing,Int} = nothing
    prec::Int = 53
    I::Union{Nothing,AcbMatrix} = nothing

    # We also store a double precision representation of the midpoint of x₁
    # as the best available double precision estimate of the solution
    solution::Union{Nothing,Vector{ComplexF64}} = nothing
end

Base.@kwdef struct ExtendedSolutionCertificate <: AbstractSolutionCertificate
    solution_candidate::AbstractVector
    certified::Bool
    real::Bool = false
    complex::Bool = false
    index::Union{Nothing,Int} = nothing
    prec::Int = 53
    # I, I′ ∈ 𝕀ℂⁿ and a certified solution has I′ ⊊ I
    I::Union{Nothing,AcbMatrix} = nothing
    I′::Union{Nothing,AcbMatrix} = nothing
    x̃::Union{Nothing,AcbMatrix} = nothing
    Y::Union{Nothing,AcbMatrix} = nothing

    # We also store a double precision representation of the midpoint of x₁
    # as the best available double precision estimate of the solution
    solution::Union{Nothing,Vector{ComplexF64}} = nothing
end


"""
    solution_candidate(certificate::AbstractSolutionCertificate)

Returns the given provided solution candidate.
"""
solution_candidate(C::AbstractSolutionCertificate) = C.solution_candidate
@deprecate initial_solution(cert::AbstractSolutionCertificate) solution_candidate(cert)

"""
    is_certified(certificate::AbstractSolutionCertificate)

Returns `true` if `certificate` is a certificate that `certified_solution_interval(certificate)`
contains a unique zero.
"""
is_certified(C::AbstractSolutionCertificate) = C.certified

"""
    is_real(certificate::AbstractSolutionCertificate)

Returns `true` if `certificate` certifies that the certified solution interval
contains a true real zero of the system.
If `false` is returned then this does not necessarily mean that the true solution is not real.
"""
is_real(C::AbstractSolutionCertificate) = C.real

"""
    is_complex(certificate::AbstractSolutionCertificate)

Returns `true` if `certificate` certifies that the certified solution interval
contains a true complex zero of the system.
"""
is_complex(C::AbstractSolutionCertificate) = C.complex

"""
    is_positive(certificate::AbstractSolutionCertificate)

Returns `true` if `is_certified(certificate)` is `true` and the unique zero contained
in `certified_solution_interval(certificate)` is real and positive.
"""
function is_positive(C::AbstractSolutionCertificate)
    if !is_certified(C) || !is_real(C)
        return false
    else
        all(i -> Arblib.is_positive(real(Arblib.ref(C.I, i, 1))), 1:length(C.I))
    end
end


"""
    is_positive(certificate::AbstractSolutionCertificate, i)

Returns `true` if `is_certified(certificate)` is `true` and `i`-th coordinate of the
unique zero contained in `certified_solution_interval(certificate)` is real and positive.
"""
function is_positive(C::AbstractSolutionCertificate, i::Integer)
    if !is_certified(C) || !is_real(C)
        return false
    else
        Arblib.is_positive(real(Arblib.ref(C.I, i, 1)))
    end
end


"""
    certified_solution_interval(certificate::AbstractSolutionCertificate)

Returns an `Arblib.AcbMatrix` representing a vector of complex intervals where a unique
zero of the system is contained in.
Returns `nothing` if `is_certified(certificate)` is `false`.
"""
certified_solution_interval(certificate::AbstractSolutionCertificate) = certificate.I
@deprecate certified_solution(certificate) certified_solution_interval(certificate)

"""
    certified_solution_interval_after_krawczyk(certificate::ExtendedSolutionCertificate)

Returns an `Arblib.AcbMatrix` representing a vector of complex intervals where a unique
zero of the system is contained in.
This is the result of applying the Krawczyk operator to `certified_solution_interval(certificate)`.
Returns `nothing` if `is_certified(certificate)` is `false`.
"""
certified_solution_interval_after_krawczyk(certificate::ExtendedSolutionCertificate) =
    certificate.I′

"""
    certificate_index(certificate::AbstractSolutionCertificate)

Return the index of the solution certificate. Here the index refers to the index of the
provided solution candidates.
"""
certificate_index(C::AbstractSolutionCertificate) = C.index

"""
    precision(certificate::AbstractSolutionCertificate)

Returns the maximal precision used to produce the given `certificate`.
"""
Base.precision(certificate::AbstractSolutionCertificate) = certificate.prec

"""
    solution_approximation(certificate::AbstractSolutionCertificate)

If `is_certified(certificate)` is `true` this returns the midpoint of the
[`certified_solution_interval`](@ref) of the given `certificate` as a `Vector{ComplexF64}`.
Returns `nothing` if `is_certified(certificate)` is `false`.
"""
solution_approximation(certificate::AbstractSolutionCertificate) = certificate.solution

"""
    krawczyk_operator_parameters(cert::ExtendedSolutionCertificate)

Returns a `NamedTuple` `(x, Y)` with the parameters of the Krawczyk operator following [BRT20].
"""
krawczyk_operator_parameters(cert::ExtendedSolutionCertificate) = (x = cert.x̃, Y = cert.Y)

function Base.show(f::IO, cert::AbstractSolutionCertificate)
    println(f, "SolutionCertificate:")
    println(f, "solution_candidate = [")
    for z in solution_candidate(cert)
        print(f, "  ")
        print(f, z)
        println(f, ",")
    end
    println(f, "]")
    print(f, "is_certified = ", is_certified(cert))
    if !isnothing(certified_solution_interval(cert))
        println(f)
        println(f, "certified_solution_interval = [")
        for z in certified_solution_interval(cert)
            print(f, "  ")
            print(f, z)
            println(f, ",")
        end
        println(f, "]")
        println(f, "precision = ", cert.prec)
        print(f, "is_real = ", is_real(cert))
    end
    if !isnothing(cert.index)
        println(f)
        print(f, "index = ", cert.index)
    end
end



"""
    squared_distance_interval(cert, reference_point::Vector{ComplexF64})

Calculate the squared distance between a solution certificate and a reference point.
"""
function squared_distance_interval(
    cert::AbstractSolutionCertificate,
    reference_point::Vector{ComplexF64},
)
    a, b = Arblib.Arf(prec = 53), Arblib.Arf(prec = 53)
    n = length(solution_candidate(cert))
    d = zero(IntervalArithmetic.Interval{Float64})
    for i = 1:n
        yᵢ = IComplexF64(cert.I[i], a, b)
        d +=
            IntervalArithmetic.sqr(real(yᵢ) - real(reference_point[i])) +
            IntervalArithmetic.sqr(imag(yᵢ) - imag(reference_point[i]))
    end
    return IntervalTrees.Interval(d.lo, d.hi)
end
"""
    squared_distance_interval(solution_candidate, reference_point::Vector{ComplexF64})

Calculate the squared distance between a solution certificate and a reference point.
"""
function squared_distance_interval(
    solution_candidate::AbstractVector{ComplexF64},
    reference_point::Vector{ComplexF64},
)
    n = length(solution_candidate)
    d = zero(IntervalArithmetic.Interval{Float64})
    for i = 1:n
        yᵢ = IComplexF64(solution_candidate[i])
        d +=
            IntervalArithmetic.sqr(real(yᵢ) - real(reference_point[i])) +
            IntervalArithmetic.sqr(imag(yᵢ) - imag(reference_point[i]))
    end
    return IntervalTrees.Interval(d.lo, d.hi)
end


"""
    struct DistinctSolutionCertificates

A struct that holds a reference point and an interval tree of distinct solution certificates.
"""
struct DistinctSolutionCertificates{S<:AbstractSolutionCertificate}
    reference_point::Vector{ComplexF64}
    distinct_tree::IntervalTrees.IntervalMap{Float64,S}
    acb_solution_candidate::AcbMatrix
end

"""
    DistinctSolutionCertificates(dim::Integer)

Create a DistinctSolutionCertificates object with a random reference point of the given dimension.
"""
DistinctSolutionCertificates(dim::Integer; kwargs...) =
    DistinctSolutionCertificates(randn(ComplexF64, dim); kwargs...)

"""
    DistinctSolutionCertificates(reference_point::Vector{ComplexF64})

Create a DistinctSolutionCertificates object with the given reference point.
"""
DistinctSolutionCertificates(
    reference_point::Vector{ComplexF64};
    extended_certificate::Bool = false,
) = DistinctSolutionCertificates(
    reference_point,
    IntervalTrees.IntervalMap{
        Float64,
        extended_certificate ? ExtendedSolutionCertificate : SolutionCertificate,
    }(),
    AcbMatrix(length(reference_point), length(reference_point)),
)

"""
    length(d::DistinctSolutionCertificates)

Return the number of distinct solution certificates in the interval tree.
"""
Base.length(d::DistinctSolutionCertificates) = length(d.distinct_tree)

Base.show(io::IO, d::DistinctSolutionCertificates) =
    print(io, "DistinctSolutionCertificates with $(length(d)) certificates")

"""
    add_certificate!(distinct_sols::DistinctSolutionCertificates, cert::SolutionCertificate)

Add a solution certificate to the interval tree if it is not a duplicate.
"""
function add_certificate!(
    distinct_sols::DistinctSolutionCertificates,
    cert::AbstractSolutionCertificate,
)
    d = squared_distance_interval(cert, distinct_sols.reference_point)

    for match in IntervalTrees.intersect(distinct_sols.distinct_tree, d)
        certᵢ = IntervalTrees.value(match)
        if Bool(Arblib.overlaps(cert.I, certᵢ.I))
            return (false, certᵢ)
            break
        end

    end

    distinct_sols.distinct_tree[d] = cert
    return (true, cert)
end

"""
    is_solution_candidate_guaranteed_duplicate(distinct_sols::DistinctSolutionCertificates,  s::Vector{ComplexF64})

Check if a solution candidate is a duplicate. This is the case if the point is within some certificate's interval.
"""
function is_solution_candidate_guaranteed_duplicate(
    distinct_sols::DistinctSolutionCertificates,
    s::Vector{ComplexF64},
)
    d = squared_distance_interval(s, distinct_sols.reference_point)

    assigned = false
    for match in IntervalTrees.intersect(distinct_sols.distinct_tree, d)
        certᵢ = IntervalTrees.value(match)
        if !assigned
            for (i, xᵢ) in enumerate(s)
                distinct_sols.acb_solution_candidate[i][] = xᵢ
            end
        end

        # Check if s is contained in certᵢ.I
        if Bool(Arblib.contains(certᵢ.I, distinct_sols.acb_solution_candidate))
            return true
            break
        end

    end

    return false
end

"""
    build_distinct_solution_certificates(certs::AbstractVector{SolutionCertificate})

Create a DistinctSolutionCertificates object from a vector of solution certificates.
"""
function build_distinct_solution_certificates(
    certs::AbstractVector{<:AbstractSolutionCertificate},
)
    distinct_sols = DistinctSolutionCertificates(length(solution_candidate(first(certs))))
    for cert in certs
        add_certificate!(distinct_sols, cert)
    end
    return distinct_sols
end


"""
    CertificationResult

The result of [`certify`](@ref) for multiple solutions.
Contains a vector of [`SolutionCertificate`](@ref) as well as a list of certificates
which correspond to the same true solution.
"""
struct CertificationResult{C<:AbstractSolutionCertificate}
    certificates::Vector{C}
    duplicates::Vector{Vector{Int}}
    slp::ModelKit.Interpreter{Vector{IComplexF64}}
    slp_jacobian::ModelKit.Interpreter{Vector{IComplexF64}}
end

"""
    certificates(R::CertificationResult)

Obtain the stored [`SolutionCertificate`](@ref)s.
"""
certificates(R::CertificationResult) = R.certificates


"""
    distinct_certificates(R::CertificationResult)

Obtain the certificates corresponding to the determined distinct solution intervals.
"""
function distinct_certificates(C::CertificationResult)
    cs = certificates(C)
    isempty(C.duplicates) && return cs
    indices = trues(length(cs))
    for d in C.duplicates, k = 2:length(d)
        indices[d[k]] = false
    end
    cs[indices]
end


"""
    distinct_solutions(R::CertificationResult)

Obtain the certificates corresponding to the determined distinct solution intervals.
"""
function distinct_solutions(C::CertificationResult)
    solution_approximation.(distinct_certificates(C))
end

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
    ncomplex_certified(R::CertificationResult)

Returns the number of certified complex solutions.
"""
ncomplex_certified(R::CertificationResult) =
    count(r -> is_certified(r) && is_complex(r), R.certificates)

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

"""
    ndistinct_complex_certified(R::CertificationResult)

Returns the number of distinct certified complex solutions.
"""
function ndistinct_complex_certified(R::CertificationResult)
    ncert = ncomplex_certified(R)
    if isempty(R.duplicates)
        return ncert
    else
        ncert - sum(R.duplicates) do dup
            is_complex(R.certificates[dup[1]]) ? length(dup) - 1 : 0
        end
    end
end

"""
    show_straight_line_program(R::CertificationResult)
    show_straight_line_program(io::IO, R::CertificationResult)

Print a representation of the used straight line program.
"""
show_straight_line_program(R::CertificationResult) = show_straight_line_program(stdout, R)
show_straight_line_program(io::IO, R::CertificationResult) =
    ModelKit.show_instructions(io, R.slp)

function Base.show(io::IO, R::CertificationResult)
    println(io, "CertificationResult")
    println(io, "===================")
    println(io, "• $(length(R.certificates)) solution candidates given")
    ncert = ncertified(R)
    print(io, "• $ncert certified solution intervals")
    nreal = nreal_certified(R)
    ncomplex = ncomplex_certified(R)
    print(io, " ($nreal real, $ncomplex complex")
    if nreal + ncomplex < ncert
        println(io, ", $(ncert - (nreal + ncomplex)) undecided)")
    else
        println(io, ")")
    end
    ndist_cert = ndistinct_certified(R)
    print(io, "• $ndist_cert distinct certified solution intervals")
    ndist_real = ndistinct_real_certified(R)
    ndist_complex = ndistinct_complex_certified(R)
    print(io, " ($ndist_real real, $ndist_complex complex")
    if ndist_real + ndist_complex < ndist_cert
        print(io, ", $(ndist_cert - (ndist_real + ndist_complex)) undecided)")
    else
        print(io, ")")
    end
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
Base.@kwdef mutable struct CertificationCache
    eval_interpreter_F64::ModelKit.Interpreter{Vector{IComplexF64}}
    jac_interpreter_F64::ModelKit.Interpreter{Vector{IComplexF64}}
    eval_interpreter_acb::ModelKit.Interpreter{AcbRefVector}
    jac_interpreter_acb::ModelKit.Interpreter{AcbRefVector}
    newton_cache::NewtonCache{MatrixWorkspace{Matrix{ComplexF64}}}
    # data for krawczyc_step
    C::Matrix{ComplexF64}
    r₀::Vector{IComplexF64}
    Δx₀::Vector{IComplexF64}
    ix̃₀::Vector{IComplexF64}
    Jx₀::Matrix{IComplexF64}
    M::Matrix{IComplexF64}
    δx::Vector{IComplexF64}

    arb_prec::Int
    arb_C::AcbRefMatrix
    arb_r₀::AcbRefMatrix # m × 1
    arb_Δx₀::AcbRefMatrix # m × 1
    arb_x̃₀::AcbRefMatrix # m × 1
    arb_x₀::AcbRefMatrix # m × 1
    arb_x₁::AcbRefMatrix # m × 1
    arb_J_x₀::AcbRefMatrix
    arb_M::AcbRefMatrix
    arb_δx::AcbRefMatrix # m × 1
    arb_mag::Arblib.Mag
end
function CertificationCache(F::AbstractSystem)
    m, n = size(F)
    m == n || error("We can only certify square systems.")
    f = System(F)
    eval_interpreter_F64 = ModelKit.interpreter(IComplexF64, f)
    jac_interpreter_F64 = ModelKit.jacobian_interpreter(IComplexF64, f)
    eval_interpreter_acb = ModelKit.interpreter(AcbRefVector, eval_interpreter_F64)
    jac_interpreter_acb = ModelKit.interpreter(AcbRefVector, jac_interpreter_F64)

    arb_prec = 128
    CertificationCache(;
        eval_interpreter_F64 = eval_interpreter_F64,
        jac_interpreter_F64 = jac_interpreter_F64,
        eval_interpreter_acb = eval_interpreter_acb,
        jac_interpreter_acb = jac_interpreter_acb,
        newton_cache = NewtonCache(F; optimize_data_structure = false),
        C = zeros(ComplexF64, m, m),
        r₀ = zeros(IComplexF64, m),
        Δx₀ = zeros(IComplexF64, m),
        ix̃₀ = zeros(IComplexF64, m),
        Jx₀ = zeros(IComplexF64, m, m),
        M = zeros(IComplexF64, m, m),
        δx = zeros(IComplexF64, m),
        arb_prec = arb_prec,
        arb_C = AcbRefMatrix(m, m; prec = arb_prec),
        arb_r₀ = AcbRefMatrix(m, 1; prec = arb_prec),
        arb_Δx₀ = AcbRefMatrix(m, 1; prec = arb_prec),
        arb_x̃₀ = AcbRefMatrix(m, 1; prec = arb_prec),
        arb_x₀ = AcbRefMatrix(m, 1; prec = arb_prec),
        arb_x₁ = AcbRefMatrix(m, 1; prec = arb_prec),
        arb_J_x₀ = AcbRefMatrix(m, m; prec = arb_prec),
        arb_M = AcbRefMatrix(m, m; prec = arb_prec),
        arb_δx = AcbRefMatrix(m, 1; prec = arb_prec),
        arb_mag = Arblib.Mag(),
    )
end

Base.setprecision(M::AcbRefMatrix, p::Int) = AcbRefMatrix(M.acb_mat, p)
function set_arb_precision!(cache::CertificationCache, p::Int)
    cache.arb_prec == p && return cache
    cache.arb_prec = p
    cache.arb_r₀ = setprecision(cache.arb_r₀, p)
    cache.arb_Δx₀ = setprecision(cache.arb_Δx₀, p)
    cache.arb_x̃₀ = setprecision(cache.arb_x̃₀, p)
    cache.arb_x₀ = setprecision(cache.arb_x₀, p)
    cache.arb_x₁ = setprecision(cache.arb_x₁, p)
    cache.arb_J_x₀ = setprecision(cache.arb_J_x₀, p)
    cache.arb_M = setprecision(cache.arb_M, p)
    cache.arb_δx = setprecision(cache.arb_δx, p)
    ModelKit.setprecision!(cache.eval_interpreter_acb, p)
    ModelKit.setprecision!(cache.jac_interpreter_acb, p)

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
* `max_precision = 256`: The maximal accuracy (in bits) that is used in the certification process.
* `compile = false`: See the [`solve`](@ref) documentation.


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
• 18 solution candidates given
• 18 certified solution intervals (4 real, 14 complex)
• 18 distinct certified solution intervals (4 real, 14 complex)
```

and see that there are indeed 18 solutions and that they are all distinct.

[^Moo77]: Moore, Ramon E. "A test for existence of solutions to nonlinear systems." SIAM Journal on Numerical Analysis 14.4 (1977): 611-615.
[^BRT20]: Breiding, P., Rose, K. and Timme, S. "Certifying zeros of polynomial systems using interval arithmetic." arXiv:2011.05000.
"""
function certify end


struct CertificationParameters
    params::Vector{ComplexF64}
    interval_params::Vector{IComplexF64}
    arb_interval_params::AcbRefVector
end

function CertificationParameters(p::AbstractVector; prec::Int = 256)
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

certification_parameters(p::AbstractVector; prec::Int = 256) = CertificationParameters(p)
certification_parameters(::Nothing; prec::Int = 256) = nothing


complexF64_params(C::CertificationParameters) = C.params
complexF64_interval_params(C::CertificationParameters) = C.interval_params
arb_interval_params(C::CertificationParameters) = C.arb_interval_params
complexF64_params(::Nothing) = nothing
complexF64_interval_params(::Nothing) = nothing
arb_interval_params(::Nothing) = nothing

function _certify(
    F::AbstractSystem,
    solution_candidates::AbstractArray{<:AbstractArray{<:Number}},
    p::Union{Nothing,CertificationParameters},
    cache::CertificationCache;
    show_progress::Bool = true,
    max_precision::Int = 256,
    check_distinct::Bool = true,
    extended_certificate::Bool = false,
    threading::Bool = true,
)
    N = length(solution_candidates)
    certificates =
        extended_certificate ? Vector{ExtendedSolutionCertificate}(undef, N) :
        Vector{SolutionCertificate}(undef, N)

    m, n = size(F)
    m == n || throw(ArgumentError("We can only certify solutions to square systems."))

    if isnothing(p) && nparameters(System(F)) > 0
        throw(ArgumentError("The given system expects parameters but none are given."))
    end

    distinct_sols = DistinctSolutionCertificates(n)
    ncertified = Threads.Atomic{Int}(0)
    nreal_certified = Threads.Atomic{Int}(0)
    nconsidered = Threads.Atomic{Int}(0)
    ndistinct = Threads.Atomic{Int}(0)
    ndistinct_real = Threads.Atomic{Int}(0)

    desc = "Certifying $N solutions... "
    barlen = min(ProgressMeter.tty_width(desc, stdout, false), 40)
    progress = nothing
    if show_progress
        progress = ProgressMeter.Progress(
            N;
            dt = 0.2,
            desc = desc,
            barlen = barlen,
            color = :green,
            output = stdout,
        )
        progress.tlast += progress.dt
    end

    duplicates_dict = Dict{Int,Vector{Int}}()
    certificates = Vector{SolutionCertificate}(undef, N)
    if threading
        distinct_lock = ReentrantLock()
        nthreads = Threads.nthreads()

        Tf = [F; [deepcopy(F) for _ = 2:nthreads]]
        Tcache = [cache; [deepcopy(cache) for _ = 2:nthreads]]
        Tp = [p; [deepcopy(p) for _ = 2:nthreads]]

        Threads.@threads for i = 1:N
            tid = Threads.threadid()

            s = solution_candidates[i]
            cert = certify_solution(Tf[tid], s, Tp[tid], Tcache[tid], i)

            certificates[i] = cert
            Threads.atomic_add!(nconsidered, 1)
            if is_certified(cert)
                Threads.atomic_add!(ncertified, 1)
                Threads.atomic_add!(nreal_certified, Int(is_real(cert)))

                @lock distinct_lock begin
                    (is_distinct, distinct_cert) = add_certificate!(distinct_sols, cert)
                    if is_distinct
                        Threads.atomic_add!(ndistinct, 1)
                        Threads.atomic_add!(ndistinct_real, Int(is_real(cert)))
                    else
                        if !haskey(duplicates_dict, distinct_cert.index)
                            duplicates_dict[distinct_cert.index] =
                                [distinct_cert.index, cert.index]
                        else
                            push!(duplicates_dict[distinct_cert.index], cert.index)
                        end
                    end

                    if !isnothing(progress)
                        update_certify_progress!(
                            progress,
                            nconsidered[],
                            ncertified[],
                            nreal_certified[],
                            ndistinct[],
                            ndistinct_real[],
                        )
                    end
                end
            else
                @lock distinct_lock begin
                    if !isnothing(progress)
                        update_certify_progress!(
                            progress,
                            nconsidered[],
                            ncertified[],
                            nreal_certified[],
                            ndistinct[],
                            ndistinct_real[],
                        )
                    end
                end
            end
        end

    else
        for i = 1:N
            s = solution_candidates[i]
            cert = certify_solution(F, s, p, cache, i)

            certificates[i] = cert
            Threads.atomic_add!(nconsidered, 1)
            if is_certified(cert)
                Threads.atomic_add!(ncertified, 1)
                Threads.atomic_add!(nreal_certified, Int(is_real(cert)))

                (is_distinct, distinct_cert) = add_certificate!(distinct_sols, cert)
                if is_distinct
                    Threads.atomic_add!(ndistinct, 1)
                    Threads.atomic_add!(ndistinct_real, Int(is_real(cert)))
                else
                    if !haskey(duplicates_dict, distinct_cert.index)
                        duplicates_dict[distinct_cert.index] =
                            [distinct_cert.index, cert.index]
                    else
                        push!(duplicates_dict[distinct_cert.index], cert.index)
                    end
                end
            end

            if !isnothing(progress)
                update_certify_progress!(
                    progress,
                    nconsidered[],
                    ncertified[],
                    nreal_certified[],
                    ndistinct[],
                    ndistinct_real[],
                )
            end
        end
    end

    duplicates = isempty(duplicates_dict) ? Vector{Int}[] : collect(values(duplicates_dict))

    CertificationResult(
        certificates,
        duplicates,
        cache.eval_interpreter_F64,
        cache.jac_interpreter_F64,
    )
end


function update_certify_progress!(
    progress,
    k,
    ncertified,
    nreal_certified,
    ndistinct,
    ndistinct_real,
)
    showvalues =
        make_certify_showvalues(k, ncertified, nreal_certified, ndistinct, ndistinct_real)
    ProgressMeter.update!(progress, k; showvalues = showvalues)

    nothing
end
@noinline function make_certify_showvalues(
    k,
    ncertified,
    nreal_certified,
    ndistinct,
    ndistinct_real,
)
    (
        ("# processed", "$k"),
        ("# certified (real)", "$(ncertified) ($nreal_certified)"),
        ("# distinct (real)", "$(ndistinct) ($ndistinct_real)"),
    )
end


function certify_solution(
    F::AbstractSystem,
    solution_candidate::AbstractVector,
    cert_params::Union{Nothing,CertificationParameters},
    cert_cache::CertificationCache,
    index::Int;
    max_precision::Int = 256,
    refine_solution::Bool = true,
    extended_certificate::Bool = false,
)
    @unpack C, arb_C, arb_x̃₀ = cert_cache

    # refine solution to machine precicision
    res = newton(
        F,
        solution_candidate,
        complexF64_params(cert_params),
        InfNorm(),
        cert_cache.newton_cache;
        rtol = 8 * eps(),
        atol = 0.0,
        # this should already be an approximate zero
        extended_precision = true,
        max_iters = 8,
    )

    x̃₀ = solution(res)
    LA.inv!(C, cert_cache.newton_cache.J)

    certified, x₁, x₀, is_real = ε_inflation_krawczyk(x̃₀, cert_params, C, cert_cache)

    if certified
        if extended_certificate
            return ExtendedSolutionCertificate(
                solution_candidate = solution_candidate,
                certified = true,
                real = is_real,
                complex = any(xi -> !(0.0 in imag(xi)), x₁),
                index = index,
                prec = 53,
                I = AcbMatrix(x₀; prec = 53),
                I′ = AcbMatrix(x₁; prec = 53),
                x̃ = AcbMatrix(x̃₀; prec = 53),
                Y = AcbMatrix(C; prec = 53),
                solution = mid.(x₁),
            )
        else
            return SolutionCertificate(
                solution_candidate = solution_candidate,
                certified = true,
                real = is_real,
                complex = any(xi -> !(0.0 in imag(xi)), x₁),
                index = index,
                prec = 53,
                I = AcbMatrix(x₁; prec = 53),
                solution = mid.(x₁),
            )
        end
    end

    return extended_prec_certify_solution(
        F,
        solution_candidate,
        x̃₀,
        cert_params,
        cert_cache,
        index;
        max_precision = max_precision,
        extended_certificate = extended_certificate,
    )
end

function extended_prec_certify_solution(
    F::AbstractSystem,
    solution_candidate::AbstractVector,
    x̃₀::AbstractVector,
    cert_params::Union{Nothing,CertificationParameters},
    cert_cache::CertificationCache,
    index::Int;
    max_precision::Int = 256,
    extended_certificate::Bool = false,
)
    @unpack C, arb_C, arb_x̃₀ = cert_cache

    n = size(C, 1)

    # If not certified we do the computation in higher precision using Arb
    prec = 128
    set_arb_precision!(cert_cache, prec)

    # We keep the same C matrix for now.
    for j = 1:n, i = 1:n
        arb_C[i, j][] = C[i, j]
    end
    for (i, x̃₀_i) in enumerate(x̃₀)
        arb_x̃₀[i][] = x̃₀_i
    end
    while (prec <= max_precision)
        certified, arb_x₁, arb_x₀, is_real =
            arb_ε_inflation_krawczyk(arb_x̃₀, cert_params, arb_C, cert_cache; prec = prec)

        if certified
            I = AcbMatrix(n, 1; prec = prec)
            Arblib.set!(I, arb_x₀)

            if extended_certificate
                I′ = AcbMatrix(n, 1; prec = prec)
                x̃ = AcbMatrix(n, 1; prec = prec)
                Y = AcbMatrix(n, n; prec = prec)
                Arblib.set!(I′, arb_x₁)
                Arblib.set!(x̃, arb_x̃₀)
                Arblib.set!(Y, arb_C)
                cert = ExtendedSolutionCertificate(
                    solution_candidate = solution_candidate,
                    certified = true,
                    real = is_real,
                    complex = any(xi -> !Arblib.contains_zero(Arblib.imagref(xi)), arb_x₁),
                    index = index,
                    prec = prec,
                    I = I,
                    I′ = I′,
                    x̃ = x̃,
                    Y = Y,
                    solution = [ComplexF64(arb_x₁[i]) for i = 1:n],
                )
                return cert
            else
                cert = SolutionCertificate(
                    solution_candidate = solution_candidate,
                    certified = true,
                    real = is_real,
                    complex = any(xi -> !Arblib.contains_zero(Arblib.imagref(xi)), arb_x₁),
                    index = index,
                    prec = prec,
                    I = I,
                    solution = [ComplexF64(arb_x₁[i]) for i = 1:n],
                )
                return cert
            end
        end

        prec += 64
        set_arb_precision!(cert_cache, prec)
    end

    # certification failed
    if extended_certificate
        return ExtendedSolutionCertificate(
            solution_candidate = solution_candidate,
            certified = false,
            index = index,
        )
    else
        return SolutionCertificate(
            solution_candidate = solution_candidate,
            certified = false,
            index = index,
        )
    end
end


function Base.Vector{ComplexF64}(A::Arblib.AcbMatrixLike)
    @assert size(A, 2) == 1
    [ComplexF64(Arblib.ref(A, i, 1).acb_ptr) for i = 1:size(A, 1)]
end

function ε_inflation_krawczyk(x̃₀, p::Union{Nothing,CertificationParameters}, C, cert_cache)
    @unpack C, r₀, Δx₀, ix̃₀, Jx₀, M, δx = cert_cache
    J_x₀ = Jx₀

    ix̃₀ .= IComplexF64.(x̃₀)
    # r₀ = F(ix̃₀)
    ModelKit.execute!(
        r₀,
        cert_cache.eval_interpreter_F64,
        ix̃₀,
        complexF64_interval_params(p),
    )
    # iΔx = C * r₀
    sqr_mul!(Δx₀, C, r₀)

    # Perform ε-inflation. We choose a different εᵢ per coordinate. This matches our strategy
    # to use a weighted norm.
    x₀ = map(x̃₀, Δx₀) do x̃₀_i, Δx₀_i
        εᵢ = 512 * max(mag(Δx₀_i), eps())
        complex(
            Interval(real(x̃₀_i) - εᵢ, real(x̃₀_i) + εᵢ),
            Interval(imag(x̃₀_i) - εᵢ, imag(x̃₀_i) + εᵢ),
        )
    end

    ModelKit.execute!(
        nothing,
        J_x₀,
        cert_cache.jac_interpreter_F64,
        x₀,
        complexF64_interval_params(p),
    )

    x₁ = similar(x₀)
    # x₁ = (x̃₀ - C * F([x̃₀])) + (I - C * J(x₀)) * (x₀ - x̃₀)
    #    = (x̃₀ - C * F([x̃₀])) - (C * J(x₀) - I) * (x₀ - x̃₀)
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
    @unpack arb_r₀, arb_Δx₀, arb_x̃₀, arb_J_x₀, arb_M, arb_δx, arb_mag, arb_x₀, arb_x₁ =
        cert_cache
    r₀, Δx₀, x̃₀, x₀, x₁, J_x₀, M, δx, m =
        arb_r₀, arb_Δx₀, arb_x̃₀, arb_x₀, arb_x₁, arb_J_x₀, arb_M, arb_δx, arb_mag

    m = arb_mag
    x̃₀′ = x̄₁ = Δx₀
    # Perform a simple newton refinement using arb_C until we cannot improve the
    # accuracy anymore
    acc = Inf
    max_iters = 10
    for i = 1:max_iters
        ModelKit.execute!(r₀, cert_cache.eval_interpreter_acb, x̃₀, arb_interval_params(p))
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
    Arblib.get_mid!(x₀, x̃₀)
    # We increase the radius by 2^(prec/4) to avoid hitting the precision limit.
    # We choose a dynamic increase to account for bad situations where any fixed choice
    # would be insufficient.
    incr_factor = exp2(div(prec, 4))
    mach_eps = exp2(-prec)
    for i = 1:n
        m[] = max(magF64(Δx₀[i], m), mach_eps) * incr_factor
        Arblib.add_error!(x₀[i], m)
    end

    ModelKit.execute!(
        nothing,
        J_x₀,
        cert_cache.jac_interpreter_acb,
        x₀,
        arb_interval_params(p),
    )

    # x₁ = (x̃₀ - C * F([x̃₀])) + (I - C * J(x₀)) * (x₀ - x̃₀)
    #    = (x̃₀ - C * F([x̃₀])) - (C * J(x₀) - I) * (x₀ - x̃₀)
    #    = (x̃₀ - Δx₀) - (C * J(x₀) - I) * (x₀ - x̃₀)

    # Define M = C * J(x₀) - I
    Arblib.mul!(M, C, J_x₀)
    for i = 1:n
        Arblib.sub!(M[i, i], M[i, i], 1)
    end

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
    max_precision::Int = 256,
    kwargs...,
)
    _certify(
        F,
        solutions(X),
        certification_parameters(parameters(X)),
        cache;
        max_precision = max_precision,
        kwargs...,
    )
end

function certify(
    F::AbstractSystem,
    r::PathResult,
    p::Union{Nothing,AbstractArray} = nothing,
    cache::CertificationCache = CertificationCache(F);
    target_parameters = nothing,
    max_precision::Int = 256,
    kwargs...,
)
    cert_params =
        certification_parameters(isnothing(p) ? target_parameters : p; prec = max_precision)
    _certify(F, [solution(r)], cert_params, cache; max_precision = max_precision, kwargs...)
end

function certify(
    F::AbstractSystem,
    r::AbstractVector{PathResult},
    p::Union{Nothing,AbstractArray} = nothing,
    cache::CertificationCache = CertificationCache(F);
    target_parameters = nothing,
    max_precision::Int = 256,
    kwargs...,
)
    cert_params =
        certification_parameters(isnothing(p) ? target_parameters : p; prec = max_precision)
    _certify(F, solution.(r), cert_params, cache; max_precision = max_precision, kwargs...)
end

function certify(
    F::AbstractSystem,
    x::AbstractVector{<:Number},
    p::Union{Nothing,AbstractArray} = nothing,
    cache::CertificationCache = CertificationCache(F);
    target_parameters = nothing,
    max_precision::Int = 256,
    kwargs...,
)
    cert_params =
        certification_parameters(isnothing(p) ? target_parameters : p; prec = max_precision)
    _certify(F, [x], cert_params, cache; max_precision = max_precision, kwargs...)
end

function certify(
    F::AbstractSystem,
    X::Result,
    p::Union{Nothing,AbstractArray} = nothing,
    cache::CertificationCache = CertificationCache(F);
    target_parameters = nothing,
    max_precision::Int = 256,
    kwargs...,
)
    params = isnothing(p) ? target_parameters : p
    cert_params = certification_parameters(params; prec = max_precision)
    _certify(
        F,
        solutions(X; only_nonsingular = true),
        cert_params,
        cache;
        max_precision = max_precision,
        kwargs...,
    )
end

function certify(
    F::AbstractSystem,
    X,
    p::Union{Nothing,AbstractArray} = nothing,
    cache::CertificationCache = CertificationCache(F);
    target_parameters = nothing,
    max_precision::Int = 256,
    kwargs...,
)
    cert_params =
        certification_parameters(isnothing(p) ? target_parameters : p; prec = max_precision)
    _certify(F, X, cert_params, cache; max_precision = max_precision, kwargs...)
end

"""
    mutable struct DistinctCertifiedSolutions{S<:AbstractSystem}

A struct to create on the fly a list of distinct certified solutions.
"""
Base.@kwdef mutable struct DistinctCertifiedSolutions{S<:AbstractSystem}
    system::Vector{S}
    parameters::Union{Vector{Nothing},Vector{CertificationParameters}}
    cache::Vector{CertificationCache}
    access_lock::ReentrantLock
    distinct_solution_certificates::DistinctSolutionCertificates{SolutionCertificate}
end

"""
    DistinctCertifiedSolutions(F::System, params; thread_safe::Bool = false)

Create a DistinctCertifiedSolutions object from a system and certification parameters.
If `thread_safe` is `true`, the object will be thread-safe.
"""
DistinctCertifiedSolutions(F::System, params; thread_safe::Bool = false) =
    DistinctCertifiedSolutions(fixed(F), params; thread_safe = thread_safe)

"""
    DistinctCertifiedSolutions(F::AbstractSystem, params; thread_safe::Bool = false)

Create a DistinctCertifiedSolutions object from an abstract system and certification parameters.
"""
function DistinctCertifiedSolutions(F::AbstractSystem, params; thread_safe::Bool = false)
    cache = CertificationCache(F)
    distinct_solution_certificates = DistinctSolutionCertificates(length(F))
    parameters = certification_parameters(params)
    access_lock = ReentrantLock()
    nthreads = thread_safe ? Threads.nthreads() : 1
    return DistinctCertifiedSolutions(
        [F; [deepcopy(F) for _ = 2:nthreads]],
        [parameters; [deepcopy(parameters) for _ = 2:nthreads]],
        [cache; [deepcopy(cache) for _ = 2:nthreads]],
        access_lock,
        distinct_solution_certificates,
    )
end

"""
    length(d::DistinctCertifiedSolutions)

Return the number of distinct solution certificates in the DistinctSolutionCertificates object.
"""
Base.length(d::DistinctCertifiedSolutions) = length(d.distinct_solution_certificates)

Base.show(io::IO, d::DistinctCertifiedSolutions) =
    print(io, "DistinctCertifiedSolutions with $(length(d)) distinct solutions")

"""
    add_solution!(d::DistinctCertifiedSolutions, sol::Vector{ComplexF64}, i::Integer = 0; certify_solution_kwargs...)

Add a solution to the DistinctSolutionCertificates object if it is not a duplicate.
"""
function add_solution!(
    d::DistinctCertifiedSolutions,
    sol::Vector{ComplexF64},
    i::Integer = 0;
    kwargs...,
)
    @lock d.access_lock begin
        if is_solution_candidate_guaranteed_duplicate(d.distinct_solution_certificates, sol)
            return (false, :duplicate)
        end
    end
    tid = length(d.system) === 1 ? 1 : Threads.threadid()
    cert =
        certify_solution(d.system[tid], sol, d.parameters[tid], d.cache[tid], i; kwargs...)
    if is_certified(cert)
        @lock d.access_lock begin
            added, _cert = add_certificate!(d.distinct_solution_certificates, cert)
            if added
                return (true, :certified_distinct)
            end
            return (false, :duplicate)
        end
    end
    return (false, :not_certified)
end

"""
    certificates(d::DistinctCertifiedSolutions)

Return a vector of solution certificates in the DistinctSolutionCertificates object.
"""
function certificates(d::DistinctCertifiedSolutions)
    return IntervalTrees.values(d.distinct_solution_certificates.distinct_tree)
end

"""
    solutions(d::DistinctCertifiedSolutions)

Return a vector of solutions in the DistinctSolutionCertificates object.
"""
function solutions(d::DistinctCertifiedSolutions)
    return map(solution_approximation, certificates(d))
end

"""
    distinct_certified_solutions(F, S, p = nothing; threading::Bool = true, show_progress::Bool = true, certify_solution_kwargs...)

Return a `DistinctCertifiedSolutions` struct containing distinct certified solutions obtained from a vector of solutions `S` of a system `F`.
Compared to `certify` this only keeps the distinct certified solutions and not all certificates. This in in particular
useful if you want to merge multiple large solution sets into one set of distinct certified solutions.
"""
function distinct_certified_solutions(
    F,
    S,
    p = nothing;
    threading::Bool = true,
    show_progress::Bool = true,
    kwargs...,
)
    d = DistinctCertifiedSolutions(F, p; thread_safe = threading)
    progress = nothing
    if show_progress
        progress =
            ProgressMeter.Progress(length(S); desc = "Processing $(length(S)) solutions")
    end
    progress_lock = ReentrantLock()
    if threading
        ndistinct = Threads.Atomic{Int}(0)
        Threads.@threads for i = 1:length(S)
            sol = S[i]
            added, = add_solution!(d, sol, i; kwargs...)
            if added
                Threads.atomic_add!(ndistinct, 1)
            end
            if show_progress
                @lock progress_lock begin
                    ProgressMeter.next!(
                        progress;
                        showvalues = [(:distinct_certified, ndistinct[])],
                    )
                end
            end
        end
    else
        ndistinct = 0
        for sol in S
            added, = add_solution!(d, sol; kwargs...)
            ndistinct += added
            if show_progress
                ProgressMeter.next!(
                    progress;
                    showvalues = [(:distinct_certified, ndistinct)],
                )
            end
        end
    end

    return d
end

