# Certification

We provide support for certifying non-singular solutions to polynomial systems.
The details of the implementation described in the article

Breiding, P., Rose, K. and Timme, S. "Certifying roots of polynomial systems using interval arithmetic." In preparation (2020)

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
ndistinct_certified
ndistinct_real_certified
save(filename, R::CertificationResult)
```

## SolutionCertificate
A [`CertificationResult`](@ref) contains in particular all [`SolutionCertificate`](@ref)s:
```@docs
SolutionCertificate
is_certified
is_real(::SolutionCertificate)
is_positive(::SolutionCertificate)
solution_candidate
certified_solution_interval
certificate_index
solution_approximation
```
