module MPoly
    
    import Base: +, - , * , ^, one, zero, zeros, getindex, setindex!, start, next, done, length, deepcopy, show, eltype

    include("poly.jl")
    include("poly_operators.jl")
    include("fixedpoly.jl")

    include("show.jl")

    export Poly, system, zeropoly, zeropolys, onepoly, vars, numvars, deg
    export generator, generators, eval_var, is_homogenous, homogenize, dehomogenize, weyl_dot
end