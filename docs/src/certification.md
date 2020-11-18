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
certificates
ncertified
nreal_certified
ncomplex_certified
ndistinct_certified
ndistinct_real_certified
ndistinct_complex_certified
save(filename, R::CertificationResult)
show_straight_line_program
```

## SolutionCertificate
A [`CertificationResult`](@ref) contains in particular all [`SolutionCertificate`](@ref)s:
```@docs
SolutionCertificate
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
