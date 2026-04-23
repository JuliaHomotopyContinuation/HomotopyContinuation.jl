# Certification

We provide support for certifying non-singular solutions to polynomial systems.
The details of the implementation described in the article

Breiding, P., Rose, K. and Timme, S. "Certifying zeros of polynomial systems using interval arithmetic." arXiv:2011.05000

## Certify
```@docs
certify
```

## CertificationResult
The result of [`certify`](@ref) is a [`CertificationResult`](@ref):

```@docs
CertificationResult
```

## IteratorCertificationResult
Certifying a [`ResultIterator`](@ref) produces a [`IteratorCertificationResult`](@ref):
```@docs
IteratorCertificationResult
bsp
```

## SolutionCertificate
A [`CertificationResult`](@ref) contains in particular all [`SolutionCertificate`](@ref)s:
```@docs
SolutionCertificate
certificates
distinct_certificates
is_certified
is_real(::SolutionCertificate)
is_complex(::SolutionCertificate)
is_positive(::SolutionCertificate)
solution_candidate
certified_solution_interval
certified_solution_interval_after_krawczyk
certificate_index
solution_approximation
```

## Incremental Certified Sets
For iteratively building a certified set of distinct solutions, use [`DistinctCertifiedSolutions`](@ref)
and [`add_solution!`](@ref). This API can be updated batch-by-batch and merged in memory with
[`merge!`](https://docs.julialang.org/en/v1/base/collections/#Base.merge!).

```@docs
DistinctCertifiedSolutions
add_solution!
distinct_certified_solutions
distinct_certified_solutions!
stats
```

## Helper functions

```@docs
ncertified
nreal_certified
ncomplex_certified
ndistinct_certified
ndistinct_real_certified
ndistinct_complex_certified
nprocessed
nduplicates
nnotcertified
nresults
nfinite
nparts
max_bucket_size
oversized_buckets
unsplittable_buckets
refinement_rounds
save(filename, R::CertificationResult)
show_straight_line_program
```