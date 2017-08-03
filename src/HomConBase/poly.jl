import ..MPoly

Poly{T} = MPoly.FixedPoly{T}

@inline evaluate(f, x) = MPoly.evaluate(f,x)
@inline gradient(f) = MPoly.gradient(f)
@inline is_homogenous(f) = MPoly.is_homogenous(f)
@inline homogenize(f) = MPoly.homogenize(f)
@inline deg(f) = MPoly.deg(f)
@inline nvars(f) = MPoly.nvars(f)
@inline weyl_dot(f,g) = MPoly.weyl_dot(f,g)