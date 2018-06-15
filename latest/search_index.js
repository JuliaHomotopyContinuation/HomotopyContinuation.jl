var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "Homotopy.Continuation.jl is a package for solving systems of polynomials with only finitely many solutions using numerical homotopy continuation. If this is your first time reading this documentation, we recommend you start with the getting started guide."
},

{
    "location": "index.html#References-1",
    "page": "Introduction",
    "title": "References",
    "category": "section",
    "text": "Pages = [\"solving.md\", \"systems.md\", \"homotopies.md\", \"predictors.md\"]"
},

{
    "location": "solving.html#",
    "page": "Solving",
    "title": "Solving",
    "category": "page",
    "text": ""
},

{
    "location": "solving.html#Solving-systems-1",
    "page": "Solving",
    "title": "Solving systems",
    "category": "section",
    "text": ""
},

{
    "location": "solving.html#HomotopyContinuation.solve",
    "page": "Solving",
    "title": "HomotopyContinuation.solve",
    "category": "function",
    "text": "solve(F::Vector{<:MP.AbstractPolynomial}; options...)\n\nSolve the system F using a total degree homotopy.\n\nsolve(G, F, start_solutions; options...)\n\nSolve the system F by tracking the each solution of G (as provided by start_solutions).\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.AffineResult",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.AffineResult",
    "category": "type",
    "text": "AffineResult <: Result\n\nThe result of an (non-homogenous) system of polynomials.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.ProjectiveResult",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.ProjectiveResult",
    "category": "type",
    "text": "ProjectiveResult <: Result\n\nThe result of a homogenous system of polynomials.\n\n\n\n"
},

{
    "location": "solving.html#solve-1",
    "page": "Solving",
    "title": "solve",
    "category": "section",
    "text": "solveDepending on the input solve returns one of the following typesAffineResult\nProjectiveResult"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.finite",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.finite",
    "category": "function",
    "text": "finite(result::AffineResult)\n\nReturn all PathResults for which the result is successfull and the contained solution is indeed a solution of the system.\n\nfinite(f::Function, result)\n\nAdditionally you can apply a transformation f on each result.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.results",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.results",
    "category": "function",
    "text": "results(result; only_real=false, realtol=1e-6, onlynonsingular=true, singulartol=1e10, onlyfinite=true)\n\nReturn all PathResults for which the given conditions apply.\n\nresults(f::Function, result; kwargs...)\n\nAdditionally you can apply a transformation f on each result.\n\nExample\n\n```julia R = solve(F)\n\nThis gives us all solutions considered real (but still as a complex vector).\n\nrealsolutions = results(solution, R, only_real=true)\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.failed",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.failed",
    "category": "function",
    "text": "failed(result)\n\nGet all results where the path tracking failed.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.atinfinity",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.atinfinity",
    "category": "function",
    "text": "atinfinity(result::AffineResult)\n\nGet all results where the solutions is at infinity.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.singular",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.singular",
    "category": "function",
    "text": "singular(result; tol=1e10)\n\nGet all singular solutions. A solution is considered singular if its windingnumber is larger than 1 or the condition number is larger than tol.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.nonsingular",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.nonsingular",
    "category": "function",
    "text": "nonsingular(result::AffineResult)\n\nReturn all PathResults for which the solution is non-singular.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.seed",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.seed",
    "category": "function",
    "text": "seed(result)\n\nThe random seed used in the computation.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.nresults",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.nresults",
    "category": "function",
    "text": "nresults(result)\n\nThe number of proper solutions.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.nfinite",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.nfinite",
    "category": "function",
    "text": "nresults(affineresult)\n\nThe number of finite solutions.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.nsingular",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.nsingular",
    "category": "function",
    "text": "nsingular(result; tol=1e10)\n\nThe number of singular solutions. A solution is considered singular if its windingnumber is larger than 1 or the condition number is larger than tol.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.nnonsingular",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.nnonsingular",
    "category": "function",
    "text": "nnonsingular(result)\n\nThe number of non-singular solutions.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.natinfinity",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.natinfinity",
    "category": "function",
    "text": "natinfinity(affineresult)\n\nThe number of solutions at infinity.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.nfailed",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.nfailed",
    "category": "function",
    "text": "nafailed(result)\n\nThe number of failed paths.\n\n\n\n"
},

{
    "location": "solving.html#Result-1",
    "page": "Solving",
    "title": "Result",
    "category": "section",
    "text": "finite\nresults\nfailed\natinfinity\nsingular\nnonsingular\nseed\nnresults\nnfinite\nnsingular\nnnonsingular\nnatinfinity\nnfailed\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.solution",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.solution",
    "category": "function",
    "text": "solution(pathresult)\n\nGet the solution of the path.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.residual",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.residual",
    "category": "function",
    "text": "residual(pathresult)\n\nGet the residual of the solution x of the path, i.e., H(x t)_infty.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.start_solution",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.start_solution",
    "category": "function",
    "text": "residual(pathresult)\n\nGet the residual of the solution x of the path, i.e., H(x t)_infty.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.issuccess",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.issuccess",
    "category": "function",
    "text": "issuccess(pathresult)\n\nChecks whether the path is successfull.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.isfailed",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.isfailed",
    "category": "function",
    "text": "isfailed(pathresult)\n\nChecks whether the path failed.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.isaffine",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.isaffine",
    "category": "function",
    "text": "isaffine(pathresult)\n\nChecks whether the path result is affine.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.isprojective",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.isprojective",
    "category": "function",
    "text": "isprojective(pathresult)\n\nChecks whether the path result is projective.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.isatinfinity",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.isatinfinity",
    "category": "function",
    "text": "isatinfinity(pathresult)\n\nChecks whether the path goes to infinity.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.issingular",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.issingular",
    "category": "function",
    "text": "issingular(pathresult; tol=1e10)\nissingular(pathresult, tol)\n\nChecks whether the path result is singular. This is true if the winding number > 1 or if the condition number of the Jacobian is larger than tol.\n\n\n\n"
},

{
    "location": "solving.html#HomotopyContinuation.Solving.isnonsingular",
    "page": "Solving",
    "title": "HomotopyContinuation.Solving.isnonsingular",
    "category": "function",
    "text": "isnonsingular(pathresult; tol=1e10)\n\nChecks whether the path result is non-singular. This is true if it is not singular.\n\n\n\n"
},

{
    "location": "solving.html#PathResult-1",
    "page": "Solving",
    "title": "PathResult",
    "category": "section",
    "text": "solution\nresidual\nstart_solution\nissuccess\nisfailed\nisaffine\nisprojective\nisatinfinity\nissingular\nisnonsingular"
},

{
    "location": "systems.html#",
    "page": "Systems",
    "title": "Systems",
    "category": "page",
    "text": ""
},

{
    "location": "systems.html#Polynomial-systems-1",
    "page": "Systems",
    "title": "Polynomial systems",
    "category": "section",
    "text": ""
},

{
    "location": "systems.html#HomotopyContinuation.Systems.FPSystem",
    "page": "Systems",
    "title": "HomotopyContinuation.Systems.FPSystem",
    "category": "type",
    "text": "FPSystem(polynomials, vars) <: AbstractSystem\n\nCreate a system using the FixedPolynomials package.\n\n\n\n"
},

{
    "location": "systems.html#HomotopyContinuation.Systems.SPSystem",
    "page": "Systems",
    "title": "HomotopyContinuation.Systems.SPSystem",
    "category": "type",
    "text": "SPSystem(polynomials, vars) <: AbstractSystem\n\nCreate a system using the StaticPolynomials package.\n\n\n\n"
},

{
    "location": "systems.html#Pre-defined-systems-1",
    "page": "Systems",
    "title": "Pre-defined systems",
    "category": "section",
    "text": "FPSystem\nSPSystem"
},

{
    "location": "systems.html#Interface-for-custom-systems-1",
    "page": "Systems",
    "title": "Interface for custom systems",
    "category": "section",
    "text": ""
},

{
    "location": "systems.html#HomotopyContinuation.SystemsBase.AbstractSystem",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.AbstractSystem",
    "category": "type",
    "text": "AbstractSystem\n\nRepresenting a system of polynomials.\n\n\n\n"
},

{
    "location": "systems.html#HomotopyContinuation.SystemsBase.AbstractSystemCache",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.AbstractSystemCache",
    "category": "type",
    "text": "AbstractSystemCache\n\nA cache to avoid allocations for the evaluation of an AbstractSystem.\n\n\n\n"
},

{
    "location": "systems.html#Abstract-types-1",
    "page": "Systems",
    "title": "Abstract types",
    "category": "section",
    "text": "Systems.AbstractSystem\nSystems.AbstractSystemCache"
},

{
    "location": "systems.html#HomotopyContinuation.SystemsBase.evaluate!",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.evaluate!",
    "category": "function",
    "text": "evaluate!(u, F::AbstractSystem, x [, cache::AbstractSystemCache])\n\nEvaluate the system F at x and store the result in u.\n\n\n\n"
},

{
    "location": "systems.html#HomotopyContinuation.SystemsBase.evaluate",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.evaluate",
    "category": "function",
    "text": "evaluate(F::AbstractSystem, x::AbstractVector [, cache::AbstractSystemCache])\n\nEvaluate the system F at x.\n\n\n\n"
},

{
    "location": "systems.html#HomotopyContinuation.SystemsBase.jacobian!",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.jacobian!",
    "category": "function",
    "text": "jacobian!(u, F::AbstractSystem, x [, cache::AbstractSystemCache])\n\nEvaluate the Jacobian of the system F at x and store the result in u.\n\n\n\n"
},

{
    "location": "systems.html#HomotopyContinuation.SystemsBase.jacobian",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.jacobian",
    "category": "function",
    "text": "jacobian(F::AbstractSystem, x [, cache::AbstractSystemCache])\n\nEvaluate the Jacobian of the system F at x.\n\n\n\n"
},

{
    "location": "systems.html#Base.size-Tuple{HomotopyContinuation.SystemsBase.AbstractSystem}",
    "page": "Systems",
    "title": "Base.size",
    "category": "method",
    "text": "Base.size(F::AbstractSystem)\n\nReturns a tuple (m, n) indicating that F is a system of m polynomials m in n variables.\n\n\n\n"
},

{
    "location": "systems.html#Mandatory-1",
    "page": "Systems",
    "title": "Mandatory",
    "category": "section",
    "text": "Systems.evaluate!\nSystems.evaluate\nSystems.jacobian!\nSystems.jacobian\nBase.size(::Systems.AbstractSystem)"
},

{
    "location": "systems.html#HomotopyContinuation.SystemsBase.cache",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.cache",
    "category": "function",
    "text": "cache(F::AbstractSystem, x)\n\nCreate a cache for the evaluation (incl. Jacobian) of F with elements of the type of x. The default implementation returns a NullCache.\n\n\n\n"
},

{
    "location": "systems.html#HomotopyContinuation.SystemsBase.evaluate_and_jacobian!",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.evaluate_and_jacobian!",
    "category": "function",
    "text": "evaluate_and_jacobian!(u, U, F, x [, cache::AbstractSystemCache])\n\nEvaluate the system F and its Jacobian at x and store the results in u (evalution) and U (Jacobian).\n\n\n\n"
},

{
    "location": "systems.html#HomotopyContinuation.SystemsBase.evaluate_and_jacobian",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemsBase.evaluate_and_jacobian",
    "category": "function",
    "text": "evaluate_and_jacobian(F::AbstractSystem, x [, cache::AbstractSystemCache])\n\nEvaluate the system F and its Jacobian at x.\n\n\n\n"
},

{
    "location": "systems.html#Optional-1",
    "page": "Systems",
    "title": "Optional",
    "category": "section",
    "text": "Systems.cache\nSystems.evaluate_and_jacobian!\nSystems.evaluate_and_jacobian"
},

{
    "location": "homotopies.html#",
    "page": "Homotopies",
    "title": "Homotopies",
    "category": "page",
    "text": ""
},

{
    "location": "homotopies.html#Homotopies-1",
    "page": "Homotopies",
    "title": "Homotopies",
    "category": "section",
    "text": ""
},

{
    "location": "homotopies.html#HomotopyContinuation.Homotopies.StraightLineHomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.Homotopies.StraightLineHomotopy",
    "category": "type",
    "text": "StraightLineHomotopy(G, F; gamma=exp(i * 2π*rand()))\n\nConstruct the homotopy H(x t) = tG(x) + (1-t)F(x).\n\n\n\n"
},

{
    "location": "homotopies.html#HomotopyContinuation.Homotopies.FixedPointHomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.Homotopies.FixedPointHomotopy",
    "category": "type",
    "text": "FixedPointHomotopy(F, x₀; gamma=exp(i * 2π*rand()))\n\nConstruct the homotopy H(x t) = (1-t)F(x) + t(x-x).\n\n\n\n"
},

{
    "location": "homotopies.html#Pre-defined-homotopies-1",
    "page": "Homotopies",
    "title": "Pre-defined homotopies",
    "category": "section",
    "text": "StraightLineHomotopy\nFixedPointHomotopy"
},

{
    "location": "homotopies.html#Interface-for-custom-homotopies-1",
    "page": "Homotopies",
    "title": "Interface for custom homotopies",
    "category": "section",
    "text": ""
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.AbstractHomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.AbstractHomotopy",
    "category": "type",
    "text": "AbstractHomotopy\n\nRepresenting a homotopy.\n\n\n\n"
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.AbstractHomotopyCache",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.AbstractHomotopyCache",
    "category": "type",
    "text": "AbstractHomotopyCache\n\nA cache to avoid allocations for the evaluation of an AbstractHomotopy.\n\n\n\n"
},

{
    "location": "homotopies.html#Abstract-types-1",
    "page": "Homotopies",
    "title": "Abstract types",
    "category": "section",
    "text": "Homotopies.AbstractHomotopy\nHomotopies.AbstractHomotopyCache"
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.evaluate!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.evaluate!",
    "category": "function",
    "text": "evaluate!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H at (x, t) and store the result in u.\n\n\n\n"
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.jacobian!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.jacobian!",
    "category": "function",
    "text": "jacobian!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the Jacobian of the homotopy H at (x, t) and store the result in u.\n\n\n\n"
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.dt!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.dt!",
    "category": "function",
    "text": "dt!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H at (x, t) and store the result in u.\n\n\n\n"
},

{
    "location": "homotopies.html#Base.size-Tuple{HomotopyContinuation.HomotopiesBase.AbstractHomotopy}",
    "page": "Homotopies",
    "title": "Base.size",
    "category": "method",
    "text": "Base.size(H::AbstractHomotopy)\n\nReturns a tuple (m, n) indicating that H is a homotopy of m polynomials m in n variables.\n\n\n\n"
},

{
    "location": "homotopies.html#Mandatory-1",
    "page": "Homotopies",
    "title": "Mandatory",
    "category": "section",
    "text": "Homotopies.evaluate!\nHomotopies.jacobian!\nHomotopies.dt!\nBase.size(::Homotopies.AbstractHomotopy)"
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.cache",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.cache",
    "category": "function",
    "text": "cache(H::AbstractHomotopy, x, t)\n\nCreate a cache for the evaluation (incl. Jacobian) of F with elements of the type of x.\n\n\n\n"
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.evaluate_and_jacobian!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.evaluate_and_jacobian!",
    "category": "function",
    "text": "evaluate_and_jacobian!(u, U, F, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H and its Jacobian at (x, t) and store the results in u (evalution) and U (Jacobian).\n\n\n\n"
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.evaluate_and_jacobian",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.evaluate_and_jacobian",
    "category": "function",
    "text": "evaluate_and_jacobian(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H and its Jacobian at (x, t).\n\n\n\n"
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.jacobian_and_dt!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.jacobian_and_dt!",
    "category": "function",
    "text": "jacobian_and_dt!(U, u, H, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H and its derivative w.r.t. t at (x, t) and store the results in u (evalution) and v (∂t).\n\n\n\n"
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.evaluate",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.evaluate",
    "category": "function",
    "text": "evaluate(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H at (x, t).\n\n\n\n"
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.jacobian",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.jacobian",
    "category": "function",
    "text": "jacobian(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the Jacobian of the homotopy H at (x, t).\n\n\n\n"
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.dt",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.dt",
    "category": "function",
    "text": "dt(H::AbstractHomotopy, x::AbstractVector, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H at (x, t).\n\n\n\n"
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.precondition!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.precondition!",
    "category": "function",
    "text": "precondition!(H::AbstractHomotopy, x, t, cache)\n\nPrepare a homotopy for things like pathtracking starting at x and t. This can modify x as well as H and anything in cache. By default this is a no-op. If H wraps another homotopy this should call precondition! on this as well.\n\n\n\n"
},

{
    "location": "homotopies.html#HomotopyContinuation.HomotopiesBase.update!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopiesBase.update!",
    "category": "function",
    "text": "update!(H::AbstractHomotopy, x, t, cache)\n\nUpdate a homotopy for new values of x and x, i.e., update an affine patch. This can modify x as well as H and anything in cache. By default this is a no-op. If H wraps another homotopy this should call update! on this as well.\n\n\n\n"
},

{
    "location": "homotopies.html#Optional-1",
    "page": "Homotopies",
    "title": "Optional",
    "category": "section",
    "text": "Homotopies.cache\nHomotopies.evaluate_and_jacobian!\nHomotopies.evaluate_and_jacobian\nHomotopies.jacobian_and_dt!\nHomotopies.evaluate\nHomotopies.jacobian\nHomotopies.dt\nHomotopies.precondition!\nHomotopies.update!"
},

{
    "location": "predictors.html#",
    "page": "Predictors",
    "title": "Predictors",
    "category": "page",
    "text": ""
},

{
    "location": "predictors.html#HomotopyContinuation.Predictors.RK4",
    "page": "Predictors",
    "title": "HomotopyContinuation.Predictors.RK4",
    "category": "type",
    "text": "RK4()\n\nThe classical Runge-Kutta predictor of order 4.\n\n\n\n"
},

{
    "location": "predictors.html#HomotopyContinuation.Predictors.Euler",
    "page": "Predictors",
    "title": "HomotopyContinuation.Predictors.Euler",
    "category": "type",
    "text": "Euler()\n\nThis uses the explicit Euler method for prediction, also known as the tangent predictor.\n\n\n\n"
},

{
    "location": "predictors.html#HomotopyContinuation.Predictors.NullPredictor",
    "page": "Predictors",
    "title": "HomotopyContinuation.Predictors.NullPredictor",
    "category": "type",
    "text": "NullPredictor()\n\nA predictor which does no prediction step, i.e., it just returns the input as its prediction.\n\n\n\n"
},

{
    "location": "predictors.html#Predictors-1",
    "page": "Predictors",
    "title": "Predictors",
    "category": "section",
    "text": "The following predictors are currently implemented.Predictors.RK4\nPredictors.Euler\nPredictors.NullPredictor"
},

]}
