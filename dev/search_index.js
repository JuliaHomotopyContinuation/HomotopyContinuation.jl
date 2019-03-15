var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "Homotopy.Continuation.jl is a package for solving systems of polynomials equations with only finitely many solutions using numerical homotopy continuation. If this is your first time reading this documentation, we recommend you start with the getting started guide."
},

{
    "location": "#Contents-1",
    "page": "Introduction",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\n  \"solving.md\",\n  \"systems.md\",\n  \"homotopies.md\",\n  \"predictors-correctors.md\",\n  \"pathtracking.md\",\n  \"newton.md\",\n  \"sorting.md\",\n  \"norms_distances.md\",\n  \"reference.md\"]"
},

{
    "location": "solving/#",
    "page": "Solving Polynomial Systems",
    "title": "Solving Polynomial Systems",
    "category": "page",
    "text": ""
},

{
    "location": "solving/#Solving-Polynomial-Systems-1",
    "page": "Solving Polynomial Systems",
    "title": "Solving Polynomial Systems",
    "category": "section",
    "text": "At the heart of the package is the solve function. It takes a bunch of different input combinations and returns an AffineResult or ProjectiveResult depending on the input.The solve function works in 3 stages.It takes the input and constructs a homotopy H(xt) such that H(x1)=G(x) and H(x0)=F(x) as well as start solutions mathcalX where for all x_1  mathcalX we have H(x_1 1)  0. This step highly depends on the concrete input you provide.\nThen all start solutions x(1) = x_1  mathcalX will be tracked to solutions x(t_e) such that H(x(t_e) t_e)  0 using a predictor-corrector scheme where t_e is a value between 0 and 1 (by default 01).\nFrom these intermediate solutions x(t_e) we start the so called endgame. This is an algorithm to predict the value x(0) for each intermediate solution.The reason for step 3 is that the final solutions can be singular which provides significant challenges for our gradient based predictor-corrector methods. For more background also check the FAQ."
},

{
    "location": "solving/#HomotopyContinuation.solve",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.solve",
    "category": "function",
    "text": "solve(F; options...)\n\nSolve the system F using a total degree homotopy. F can be\n\nVector{<:MultivariatePolynomials.AbstractPolynomial} (e.g. constructed by @polyvar)\nAbstractSystem (the system has to represent a homogenous polynomial system.)\n\nExample\n\nAssume we want to solve the system F(xy) = (x^2+y^2+1 2x+3y-1).\n\n@polyvar x y\nsolve([x^2+y^2+1, 2x+3y-1])\n\nIf you polynomial system is already homogenous, but you would like to consider it as an affine system you can do\n\n@polyvar x y z\nsolve([x^2+y^2+z^2, 2x+3y-z], homvar=z)\n\nThis would result in the same result as solve([x^2+y^2+1, 2x+3y-1]).\n\nTo solve F by a custom AbstractSystem you can do\n\n@polyvar x y z\n# The system `F` has to be homgoenous system\nF = SPSystem([x^2+y^2+z^2, 2x+3y-z]) # SPSystem <: AbstractSystem\n# To solve the original affine system we have to tell that the homogenization variable has index 3\nsolve(F, homvar=3)\n\nor equivalently (in this case) by\n\nsolve([x^2+y^2+z^2, 2x+3y-z], system=SPSystem)\n\nStart Target Homotopy\n\nsolve(G, F, start_solutions; options...)\n\nSolve the system F by tracking the each provided solution of G (as provided by start_solutions).\n\nExample\n\n@polyvar x y\nG = [x^2-1,y-1]\nF = [x^2+y^2+z^2, 2x+3y-z]\nsolve(G, F, [[1, 1], [-1, 1]])\n\nParameter Homotopy\n\nsolve(F::Vector{<:MultivariatePolynomials.AbstractPolynomial},\n    startsolutions; parameters::Vector{<:MP.AbstractVariable}, p₁, p₀, γ₁=nothing, γ₀=nothing)\n\nSolve the parameter homotopy\n\nH(x t) = F(x (tγ₁p₁+(1-t)γ₀p₀)  (tγ₁+(1-t)γ₀))\n\n, where p₁ and p₀ are a vector of parameter values for F and γ₁ and γ₀ are complex numbers. If γ₁ or γ₀ is nothing, it is assumed that γ₁ and γ₀ are 1. The input parameters specifies the parameter variables of F which should be considered as parameters. Neccessarily, length(parameters) == length(p₁) == length(p₀).\n\nsolve(F::Vector{<:MultivariatePolynomials.AbstractPolynomial},\n        startsolutions; parameters::Vector{<:MP.AbstractVariable},\n        startparameters, targetparameters,\n        startgamma=randn(ComplexF64), targetgamma=randn(ComplexF64))\n\nThis is a non-unicode variant where γ₁=start_parameters, γ₀=target_parameters,     γ₁=start_gamma, γ₀=target_gamma.\n\nExample\n\nWe want to solve a parameter homotopy H(xt) = F(x t1 0+(1-t)2 4) where\n\nF(x a) = (x₁^2-a₁ x₁x₂-a₁+a₂)\n\nand let\'s say we are only intersted in tracking of 11. This can be accomplished as follows\n\n@polyvar x[1:2] a[1:2]\nF = [x[1]^2-a[1], x[1]*x[2]-a[1]+a[2]]\nstartsolutions = [[1, 1]]\nsolve(F, startsolutions, parameters=a, p₁=p₁, p₀=p₀)\n# If you don\'t like unicode this is also possible\nsolve(F, startsolutions, parameters=a, startparameters=p₁, targetparameters=p₀)\n\nAbstract Homotopy\n\nsolve(H::AbstractHomotopy, start_solutions; options...)\n\nSolve the homotopy H by tracking the each solution of H( t) (as provided by start_solutions) from t=1 to t=0. Note that H has to be a homotopy between homogenous polynomial systems. If it should be considered as an affine system indicate which is the index of the homogenization variable, e.g. solve(H, startsolutions, homvar=3) if the third variable is the homogenization variable.\n\nOptions\n\nGeneral options:\n\nsystem::AbstractSystem: A constructor to assemble a AbstractSystem. The default is FPSystem. This constructor is only applied to the input of solve. The constructor is called with system(polynomials, variables) where polynomials is a vector of MultivariatePolynomials.AbstractPolynomials and variables determines the variable ordering.\nhomotopy::AbstractHomotopy: A constructor to construct a AbstractHomotopy. The default is StraightLineHomotopy. The constructor is called with homotopy(start, target) where start and target are homogenous AbstractSystems.\nseed::Int: The random seed used during the computations.\nhomvar::Union{Int,MultivariatePolynomials.AbstractVariable}: This considers the homogenous system F as an affine system which was homogenized by homvar. If F is an AbstractSystem homvar is the index (i.e. Int) of the homogenization variable. If F is an AbstractVariables (e.g. created by @polyvar x) homvar is the actual variable used in the system F.\nendgame_start=0.1: The value of t for which the endgame is started.\nreport_progress=true: Whether a progress bar should be printed to STDOUT.\nthreading=true: Enable or disable multi-threading.\n\nPathtracking specific:\n\ncorrector::AbstractCorrector: The corrector used during in the predictor-corrector scheme. The default is NewtonCorrector.\ncorrector_maxiters=2: The maximal number of correction steps in a single step.\npredictor::AbstractPredictor: The predictor used during in the predictor-corrector scheme. The default is Heun.\nrefinement_maxiters=corrector_maxiters: The maximal number of correction steps used to refine the final value.\nrefinement_tol=1e-8: The precision used to refine the final value.\ntol=1e-7: The precision used to track a value.\ninitial_steplength=0.1: The initial step size for the predictor.\nminimal_steplength=1e-14: The minimal step size. If the size of step is below this the path is considered failed.\nmaxiters=1000: The maximal number of steps per path.\n\nEndgame specific options\n\ncauchy_loop_closed_tolerance=1e-3: The tolerance for which is used to determine whether a loop is closed. The distance between endpoints is normalized by the maximal difference between any point in the loop and the starting point.\ncauchy_samples_per_loop=6: The number of samples used to predict an endpoint. A higher number of samples should result in a better approximation. Note that the error should be roughly t^n where t is the current time of the loop and n is cauchy_samples_per_loop.\negtol=1e-10: This is the tolerance necessary to declare the endgame converged.\nmaxnorm=1e5: If our original problem is affine we declare a path at infinity if the infinity norm with respect to the standard patch is larger than maxnorm.\nmaxwindingnumber=15: The maximal windingnumber we try to find using Cauchys integral formula.\nmax_extrapolation_samples=4: During the endgame a Richardson extrapolation is used to improve the accuracy of certain approximations. This is the maximal number of samples used for this.\nminradius=1e-15: A path is declared false if the endgame didn\'t finished until then.\nsampling_factor=0.5: During the endgame we approach 0 by the geometric series h^kR₀ where h is sampling_factor and R₀ the endgame start provided in runendgame.\nmaxiters_per_step=100: The maximal number of steps bewtween two samples.\n\n\n\n\n\n"
},

{
    "location": "solving/#The-*solve*-function-1",
    "page": "Solving Polynomial Systems",
    "title": "The solve function",
    "category": "section",
    "text": "solve"
},

{
    "location": "solving/#HomotopyContinuation.AffineResult",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.AffineResult",
    "category": "type",
    "text": "AffineResult <: Result\n\nThe result of an (non-homogenous) system of polynomials.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.ProjectiveResult",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.ProjectiveResult",
    "category": "type",
    "text": "ProjectiveResult <: Result\n\nThe result of a homogenous system of polynomials.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.results",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.results",
    "category": "function",
    "text": "results(result; onlyreal=false, realtol=1e-6, onlynonsingular=false, onlysigular=false, singulartol=1e14, onlyfinite=true)\n\nReturn all PathResults for which the given conditions apply.\n\nExample\n\nR = solve(F)\n\n# This gives us all PathResults considered non-singular and real (but still as a complex vector).\nrealsolutions = results(R, onlyreal=true, onlynonsingular=true)\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.mapresults",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.mapresults",
    "category": "function",
    "text": "mapresults(f::Function, result; conditions...)\n\nApply the function f to all PathResults for which the given conditions apply. For the possible conditions see results.\n\nExample\n\n# This gives us all solutions considered real (but still as a complex vector).\nrealsolutions = results(solution, R, onlyreal=true)\n\n\n\n\n\nmapresults(f, result::MonodromyResult; onlyreal=false, realtol=1e-6)\n\nApply the function f to all entries of MonodromyResult for which the given conditions apply.\n\nExample\n\n# This gives us all solutions considered real (but still as a complex vector).\nrealsolutions = mapresults(solution, R, onlyreal=true)\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.solutions",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.solutions",
    "category": "function",
    "text": "solutions(result; conditions...)\n\nReturn all solution (as Vectors) for which the given conditions apply. For the possible conditions see results.\n\nExample\n\njulia> @polyvar x y\njulia> result = solve([(x-2)y, y+x+3]);\njulia> solutions(result)\n[[2.0+0.0im, -5.0+0.0im], [-3.0+0.0im, 0.0+0.0im]]\n\n\n\n\n\nsolutions(loop::Loop)\n\nGet the solutions of the loop.\n\n\n\n\n\nsolutions(result::MonodromyResult; onlyreal=false, realtol=1e-6)\n\nReturn all solutions (as SVectors) for which the given conditions apply.\n\nExample\n\nrealsolutions = solutions(R, onlyreal=true)\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.realsolutions",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.realsolutions",
    "category": "function",
    "text": "realsolutions(result; tol=1e-6, conditions...)\n\nReturn all real solution (as Vectors of reals) for which the given conditions apply. For the possible conditions see results. Note that onlyreal is always true and realtol is now tol.\n\nExample\n\njulia> @polyvar x y\njulia> result = solve([(x-2)y, y+x+3]);\njulia> realsolutions(result)\n[[2.0, -5.0], [-3.0, 0.0]]\n\n\n\n\n\nrealsolutions(res::MonodromyResult; tol=1e-6)\n\nReturns the solutions of res whose imaginary part has norm less than 1e-6.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.uniquesolutions",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.uniquesolutions",
    "category": "function",
    "text": "uniquesolutions(R::Result; tol=1e-6, multiplicities=false, conditions...)\n\nReturn all unique solutions. If multiplicities is true, then all unique solutions with their correspnding multiplicities as pairs (s, m) where s is the solution and m the multiplicity are returned. For the possible conditions see results.\n\nExample\n\njulia> @polyvar x;\njulia> uniquesolutions([(x-3)^3*(x+2)], multiplicities=true)\n[([3.0+0.0im], 3), ([-2.0+0.0im], 1)]\njulia> uniquesolutions([(x-3)^3*(x+2)])\n[[3.0+0.0im], [-2.0+0.0im]]\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.finite",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.finite",
    "category": "function",
    "text": "finite(result::AffineResults; conditions...)\n\nReturn all PathResults for which the solution is finite. This is just a shorthand for results(R; onlyfinite=true, conditions...). For the possible conditions see results.\n\n\n\n\n\n"
},

{
    "location": "solving/#Base.real-Tuple{Union{Result, Array{#s139,1} where #s139<:PathResult}}",
    "page": "Solving Polynomial Systems",
    "title": "Base.real",
    "category": "method",
    "text": "real(result, tol=1e-6)\n\nGet all results where the solutions are real with the given tolerance tol. See isreal for details regarding the determination of \'realness\'.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.atinfinity",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.atinfinity",
    "category": "function",
    "text": "atinfinity(result::AffineResult)\n\nGet all results where the solutions is at infinity.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.singular",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.singular",
    "category": "function",
    "text": "singular(result::Results; conditions...)\n\nReturn all PathResults for which the solution is singular. This is just a shorthand for results(R; onlysingular=true, conditions...). For the possible conditions see results.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.nonsingular",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.nonsingular",
    "category": "function",
    "text": "nonsingular(result::Results; conditions...)\n\nReturn all PathResults for which the solution is non-singular. This is just a shorthand for results(R; onlynonsingular=true, conditions...). For the possible conditions see results.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.failed",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.failed",
    "category": "function",
    "text": "failed(result)\n\nGet all results where the path tracking failed.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.multiplicities-Tuple{Union{Result, Array{#s139,1} where #s139<:PathResult}}",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.multiplicities",
    "category": "method",
    "text": "multiplicities(V::Results; tol=1e-6)\n\nReturns a Vector of Vector{PathResult}s grouping the PathResults whose solutions appear with multiplicities greater 1 in \'V\'. Two solutions are regarded as equal, when their pairwise distance is less than \'tol\'.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.seed",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.seed",
    "category": "function",
    "text": "seed(result)\n\nThe random seed used in the computation.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.nresults",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.nresults",
    "category": "function",
    "text": "nresults(result; onlyreal=false, realtol=1e-6, onlynonsingular=false, singulartol=1e14, onlyfinite=true)\n\nThe number of solutions which satisfy the corresponding predicates.\n\nExample\n\nresult = solve(F)\n# Get all non-singular results where all imaginary parts are smaller than 1e-8\nnresults(result, onlyreal=true, realtol=1e-8, onlynonsingular=true)\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.nfinite",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.nfinite",
    "category": "function",
    "text": "nfinite(affineresult)\n\nThe number of finite solutions.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.nreal",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.nreal",
    "category": "function",
    "text": "nreal(result; tol=1e-6)\n\nThe number of real solutions where all imaginary parts of each solution are smaller than tol.\n\n\n\n\n\nnreal(res::MonodromyResult; tol=1e-6)\n\nCounts how many solutions of res have imaginary part norm less than 1e-6.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.nsingular",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.nsingular",
    "category": "function",
    "text": "nsingular(result; tol=1e10)\n\nThe number of singular solutions. A solution is considered singular if its windingnumber is larger than 1 or the condition number is larger than tol.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.nnonsingular",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.nnonsingular",
    "category": "function",
    "text": "nnonsingular(result; tol=1e-10)\n\nThe number of non-singular solutions.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.natinfinity",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.natinfinity",
    "category": "function",
    "text": "natinfinity(affineresult)\n\nThe number of solutions at infinity.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.nfailed",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.nfailed",
    "category": "function",
    "text": "nafailed(result)\n\nThe number of failed paths.\n\n\n\n\n\n"
},

{
    "location": "solving/#The-result-of-*solve*-1",
    "page": "Solving Polynomial Systems",
    "title": "The result of solve",
    "category": "section",
    "text": "Depending on the input solve returns one of the following typesAffineResult\nProjectiveResultA Result is a wrapper around the results of each single path (PathResult) and it contains some additional informations like the used random seed for the computation.In order to analyze a Result we provide the following helper functionsresults\nmapresults\nsolutions\nrealsolutions\nuniquesolutions\nfinite\nBase.real(::HomotopyContinuation.Results)\natinfinity\nsingular\nnonsingular\nfailed\nmultiplicities(::HomotopyContinuation.Results)\nseedIf you are interested in the number of solutions of a certain kind we also provide the following helper functions.nresults\nnfinite\nnreal\nnsingular\nnnonsingular\nnatinfinity\nnfailed"
},

{
    "location": "solving/#HomotopyContinuation.PathResult",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.PathResult",
    "category": "type",
    "text": "PathResult(startvalue, pathtracker_result, endgamer_result, solver)\n\nA PathResult is the result of the tracking of a path (inclusive endgame). Its fields are\n\nreturncode: One of :success, :at_infinity or any error code from the EndgamerResult\nsolution::Vector{T}: The solution vector. If the algorithm computed in projective space and the solution is at infinity then the projective solution is given. Otherwise an affine solution is given if the startvalue was affine and a projective solution is given if the startvalue was projective.\nresidual::Float64: The value of the infinity norm of H(solution, 0).\nnewton_residual: The value of the 2-norm of J_H(textsolution)^-1H(textsolution 0)\ncondition_number: This is the condition number of the row-equilibrated Jacobian at the solution. A high condition number indicates a singularity.\nwindingnumber: The estimated winding number\nangle_to_infinity: The angle to infinity is the angle of the solution to the hyperplane where the homogenizing coordinate is 0.\nreal_solution: Indicates whether the solution is real given the defined tolerance at_infinity_tol (from the solver options).\nstartvalue: The startvalue of the path\niterations: The number of iterations the pathtracker needed.\nendgame_iterations: The number of steps in the geometric series the endgamer did.\nnpredictions: The number of predictions the endgamer did.\npredictions: The predictions of the endgamer.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.solution",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.solution",
    "category": "function",
    "text": "solution(pathresult)\n\nGet the solution of the path.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.residual",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.residual",
    "category": "function",
    "text": "residual(pathresult)\n\nGet the residual of the solution x of the path, i.e., H(x t)_infty.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.startsolution",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.startsolution",
    "category": "function",
    "text": "startsolution(pathresult)\n\nGet the start solution of the solution x of the path.\n\n\n\n\n\n"
},

{
    "location": "solving/#Base.isreal-Tuple{PathResult}",
    "page": "Solving Polynomial Systems",
    "title": "Base.isreal",
    "category": "method",
    "text": "isreal(pathresult; tol=1e-6)\nisreal(pathresult, tol)\n\nDetermine whether infinity norm of the imaginary part of the so\n\n\n\n\n\n"
},

{
    "location": "solving/#LinearAlgebra.issuccess-Tuple{PathResult}",
    "page": "Solving Polynomial Systems",
    "title": "LinearAlgebra.issuccess",
    "category": "method",
    "text": "issuccess(pathresult)\n\nChecks whether the path is successfull.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.isfailed",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.isfailed",
    "category": "function",
    "text": "isfailed(pathresult)\n\nChecks whether the path failed.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.isaffine",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.isaffine",
    "category": "function",
    "text": "isaffine(pathresult)\n\nChecks whether the path result is affine.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.isprojective",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.isprojective",
    "category": "function",
    "text": "isprojective(pathresult)\n\nChecks whether the path result is projective.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.isatinfinity",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.isatinfinity",
    "category": "function",
    "text": "isatinfinity(pathresult)\n\nChecks whether the path goes to infinity.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.issingular",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.issingular",
    "category": "function",
    "text": "issingular(pathresult;tol=1e14)\nissingular(pathresult, tol)\n\nChecks whether the path result is singular. This is true if the winding number > 1 or if the condition number of the Jacobian is larger than tol.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.isnonsingular",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.isnonsingular",
    "category": "function",
    "text": "isnonsingular(pathresult;tol=1e14)\n\nChecks whether the path result is non-singular. This is true if it is not singular.\n\n\n\n\n\n"
},

{
    "location": "solving/#PathResult-1",
    "page": "Solving Polynomial Systems",
    "title": "PathResult",
    "category": "section",
    "text": "For each path we return a PathResult containing the detailed information about the single path.PathResultThe following helper functions are providedsolution\nresidual\nstartsolution\nBase.isreal(::PathResult)\nLinearAlgebra.issuccess(::PathResult)\nisfailed\nisaffine\nisprojective\nisatinfinity\nissingular\nisnonsingular"
},

{
    "location": "solving/#HomotopyContinuation.monodromy_solve",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.monodromy_solve",
    "category": "function",
    "text": "monodromy_solve(F, sols, p; parameters=..., options..., pathtrackerkwargs...)\n\nSolve a polynomial system F(x;p) with specified parameters and initial solutions sols by monodromy techniques. This makes loops in the parameter space of F to find new solutions.\n\nOptions\n\ntarget_solutions_count=nothing: The computations are stopped if this number of solutions is reached.\ndone_callback=always_false: A callback to end the computation early. This function takes 2 arguments. The first one is the new solution x and the second one are all current solutions (including x). Return true if the compuation is done.\nmaximal_number_of_iterations_without_progress::Int=10: The maximal number of iterations (i.e. loops generated) without any progress.\ngroup_action=nothing: A function taking one solution and returning other solutions if there is a constructive way to obtain them, e.g. by symmetry.\nstrategy: The strategy used to create loops. If F only depends linearly on p this will be Petal. Otherwise this will be Triangle with weights if F is a real system.\nshowprogress=true: Enable a progress meter.\naccuracy::Float64=1e-6: The tolerance with which it is decided whether two solutions are identical.\ngroup_actions=nothing: If there is more than one group action you can use this to chain the application of them. For example if you have two group actions foo and bar you can set group_actions=[foo, bar]. See GroupActions for details regarding the application rules.\ngroup_action_on_all_nodes=false: By default the group_action(s) are only applied on the solutions with the main parameter p. If this is enabled then it is applied for every parameter q.\nparameter_sampler=independent_normal: A function taking the parameter p and returning a new random parameter q. By default each entry of the parameter vector is drawn independently from the univariate normal distribution.\nequivalence_classes=true: This only applies if there is at least one group action supplied. We then consider two solutions in the same equivalence class if we can transform one to the other by the supplied group actions. We only track one solution per equivalence class.\ntimeout=float(typemax(Int)): The maximal number of seconds the computation is allowed to run.\nminimal_number_of_solutions: The minimal number of solutions before a stopping heuristic is applied. By default this is half of target_solutions_count if applicable otherwise 2.\n\n\n\n\n\n"
},

{
    "location": "solving/#Monodromy-Solve-1",
    "page": "Solving Polynomial Systems",
    "title": "Monodromy Solve",
    "category": "section",
    "text": "monodromy_solve"
},

{
    "location": "solving/#HomotopyContinuation.Triangle",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Triangle",
    "category": "type",
    "text": "Triangle(;weights=true)\n\nA triangle is a loop consisting of the main node and two addtional nodes. If weights is true the edges are equipped with additional random weights. Note that this is usually only necessary for real parameters.\n\n\n\n\n\n"
},

{
    "location": "solving/#HomotopyContinuation.Petal",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.Petal",
    "category": "type",
    "text": "Petal()\n\nA petal is a loop consisting of the main node and one other node connected by two edges with different random weights.\n\n\n\n\n\n"
},

{
    "location": "solving/#Strategies-1",
    "page": "Solving Polynomial Systems",
    "title": "Strategies",
    "category": "section",
    "text": "Triangle\nPetal"
},

{
    "location": "solving/#HomotopyContinuation.GroupActions",
    "page": "Solving Polynomial Systems",
    "title": "HomotopyContinuation.GroupActions",
    "category": "type",
    "text": "GroupActions(actions::Function...)\n\nStore a bunch of group actions (f1, f2, f3, ...). Each action has to return a tuple. The actions are applied in the following sense\n\nf1 is applied on the original solution s\nf2 is applied on s and the results of 1\nf3 is applied on s and the results of 1) and 2)\n\nand so on\n\nExample\n\njulia> f1(s) = (s * s,);\n\njulia> f2(s) = (2s, -s, 5s);\n\njulia> f3(s) = (s + 1,);\n\njulia> GroupActions(f1)(3)\n(9,)\n\njulia> GroupActions(f1,f2)(3)\n(9, 18, -9, 45)\n\njulia> GroupActions(f1,f2, f3)(3)\n(9, 18, -9, 45, 10, 19, -8, 46)\n\n\n\n\n\n"
},

{
    "location": "solving/#GroupActions-1",
    "page": "Solving Polynomial Systems",
    "title": "GroupActions",
    "category": "section",
    "text": "If there is a group acting on the solution set of the polynomial system this can provided with the group_action keyword for single group actions or with the group_actions keyword for compositions of group actions. These will be internally transformed into GroupActions.GroupActions"
},

{
    "location": "systems/#",
    "page": "Systems",
    "title": "Systems",
    "category": "page",
    "text": ""
},

{
    "location": "systems/#Polynomial-systems-1",
    "page": "Systems",
    "title": "Polynomial systems",
    "category": "section",
    "text": "Polynomial systems can be represented in numerous ways in a computer and each representation has certain tradeoffs. For our purposes the most important thing is that it is fast to evaluate the system. Therefore we automatically convert an input given by DynamicPolynomials to another representation more suitable for numerically evaluations. The default is currently FPSystem."
},

{
    "location": "systems/#HomotopyContinuation.FPSystem",
    "page": "Systems",
    "title": "HomotopyContinuation.FPSystem",
    "category": "type",
    "text": "FPSystem(polynomials, vars) <: AbstractSystem\n\nCreate a polynomial system using the FixedPolynomials package.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.SPSystem",
    "page": "Systems",
    "title": "HomotopyContinuation.SPSystem",
    "category": "type",
    "text": "SPSystem(polynomials, vars) <: AbstractSystem\n\nCreate a system using the StaticPolynomials package. Note that StaticPolynomials leverages Julias metaprogramming capabilities to automatically generate functions to evaluate the system and its Jacobian. These generated functions are very fast but at the cost of possibly large compile times. The compile time depends on the size of the support of the polynomial system. If you intend to solve a large system or you need to solve a system with the same support but different coefficients even large compile times can be worthwile. As a general rule of thumb this usually is twice as fast as solving the same system using FPSystem.\n\nExample\n\nYou can use SPSystem as follows with solve\n\n@polyvar x y\nF = [x^2+3y^4-2, 2y^2+3x*y+4]\nsolve(F, system=SPSystem)\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.FixedHomotopy",
    "page": "Systems",
    "title": "HomotopyContinuation.FixedHomotopy",
    "category": "type",
    "text": "FixedHomotopy(H, t) <: AbstractSystem\n\nFix a homotopy H(x,t) at t\n\n\n\n\n\n"
},

{
    "location": "systems/#Default-systems-1",
    "page": "Systems",
    "title": "Default systems",
    "category": "section",
    "text": "We provide the following systems by default.FPSystem\nSPSystem\nFixedHomotopy"
},

{
    "location": "systems/#Interface-for-custom-systems-1",
    "page": "Systems",
    "title": "Interface for custom systems",
    "category": "section",
    "text": "The great thing is that you are not limited to the systems provided by default. Maybe your polynomial system has a particular structure which you want to use to efficiently evaluate it. For this you can define your own homotopy by defining a struct with super type AbstractSystem. For this the following interface has to be defined."
},

{
    "location": "systems/#HomotopyContinuation.AbstractSystem",
    "page": "Systems",
    "title": "HomotopyContinuation.AbstractSystem",
    "category": "type",
    "text": "AbstractSystem\n\nRepresenting a system of polynomials.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.AbstractSystemCache",
    "page": "Systems",
    "title": "HomotopyContinuation.AbstractSystemCache",
    "category": "type",
    "text": "AbstractSystemCache\n\nA cache to avoid allocations for the evaluation of an AbstractSystem.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.SystemNullCache",
    "page": "Systems",
    "title": "HomotopyContinuation.SystemNullCache",
    "category": "type",
    "text": "SystemNullCache\n\nAn empty cache if no cache is necessary.\n\n\n\n\n\n"
},

{
    "location": "systems/#Types-1",
    "page": "Systems",
    "title": "Types",
    "category": "section",
    "text": "AbstractSystem\nAbstractSystemCache\nSystemNullCache"
},

{
    "location": "systems/#HomotopyContinuation.cache-Tuple{AbstractSystem,Vararg{Any,N} where N}",
    "page": "Systems",
    "title": "HomotopyContinuation.cache",
    "category": "method",
    "text": "cache(F::AbstractSystem, x)::AbstractSystemCache\n\nCreate a cache for the evaluation (incl. Jacobian) of F with elements of the type of x.\n\ncache(F::AbstractSystem, x, p)::AbstractSystemCache\n\nCreate a cache for the evaluation (incl. Jacobian) of F with elements of the type of x and parameters p.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.evaluate!-Tuple{Any,AbstractSystem,Vararg{Any,N} where N}",
    "page": "Systems",
    "title": "HomotopyContinuation.evaluate!",
    "category": "method",
    "text": "evaluate!(u, F::AbstractSystem, x, cache::AbstractSystemCache)\n\nEvaluate the system F at x and store the result in u.\n\nevaluate!(u, F::AbstractSystem, x, p, cache::AbstractSystemCache)\n\nEvaluate the system F at x and parameters p and store the result in u.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.evaluate",
    "page": "Systems",
    "title": "HomotopyContinuation.evaluate",
    "category": "function",
    "text": "evaluate(F::AbstractSystem, x::AbstractVector, cache=cache(F, x))\n\nEvaluate the system F at x.\n\nevaluate(F::AbstractSystem, x::AbstractVector, p, cache=cache(F, x))\n\nEvaluate the system F at x and parameters p.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.jacobian!-Tuple{Any,AbstractSystem,Vararg{Any,N} where N}",
    "page": "Systems",
    "title": "HomotopyContinuation.jacobian!",
    "category": "method",
    "text": "jacobian!(u, F::AbstractSystem, x , cache::AbstractSystemCache)\n\nEvaluate the Jacobian of the system F at x and store the result in u.\n\njacobian!(u, F::AbstractSystem, x , p, cache::AbstractSystemCache)\n\nEvaluate the Jacobian of the system F at x and parameters p and store the result in u.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.jacobian",
    "page": "Systems",
    "title": "HomotopyContinuation.jacobian",
    "category": "function",
    "text": "jacobian(F::AbstractSystem, x, cache=cache(F, x))\n\nEvaluate the Jacobian of the system F at x.\n\njacobian(F::AbstractSystem, x , p, cache::AbstractSystemCache)\n\nEvaluate the Jacobian of the system F at x and parameters p.\n\n\n\n\n\n"
},

{
    "location": "systems/#Base.size-Tuple{AbstractSystem}",
    "page": "Systems",
    "title": "Base.size",
    "category": "method",
    "text": "Base.size(F::AbstractSystem)\n\nReturns a tuple (m, n) indicating that F is a system of m polynomials m in n variables.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.differentiate_parameters!",
    "page": "Systems",
    "title": "HomotopyContinuation.differentiate_parameters!",
    "category": "function",
    "text": "differentiate_parameters!(u, F::AbstractSystem, x, p, cache::AbstractSystemCache)\n\nEvaluate the Jacobian of the system F at x and parameters p w.r.t. the parameters and store the result in u.\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.differentiate_parameters",
    "page": "Systems",
    "title": "HomotopyContinuation.differentiate_parameters",
    "category": "function",
    "text": "differentiate_parameters(F::AbstractSystem, x, p, cache=cache(F, x))\n\nEvaluate the Jacobian of the system F at x and parameters p w.r.t. the parameters\n\n\n\n\n\n"
},

{
    "location": "systems/#Mandatory-1",
    "page": "Systems",
    "title": "Mandatory",
    "category": "section",
    "text": "The following methods are mandatory to implement.cache(F::AbstractSystem, args...)\nevaluate!(u, F::AbstractSystem, args...)\nevaluate(F::AbstractSystem, x, c::AbstractSystemCache=cache(F, x))\njacobian!(u, F::AbstractSystem, args...)\njacobian(F::AbstractSystem, x, c::AbstractSystemCache=cache(F, x))\nBase.size(::AbstractSystem)Additionally if the system should support parameter homotopies it needs to supportdifferentiate_parameters!\ndifferentiate_parameters"
},

{
    "location": "systems/#HomotopyContinuation.evaluate_and_jacobian!-Tuple{Any,Any,AbstractSystem,Any,AbstractSystemCache}",
    "page": "Systems",
    "title": "HomotopyContinuation.evaluate_and_jacobian!",
    "category": "method",
    "text": "evaluate_and_jacobian!(u, U, F, x , cache::AbstractSystemCache)\n\nEvaluate the system F and its Jacobian at x and store the results in u (evalution) and U (Jacobian).\n\n\n\n\n\n"
},

{
    "location": "systems/#HomotopyContinuation.evaluate_and_jacobian!-Tuple{Any,Any,AbstractSystem,Any,Any,AbstractSystemCache}",
    "page": "Systems",
    "title": "HomotopyContinuation.evaluate_and_jacobian!",
    "category": "method",
    "text": "evaluate_and_jacobian!(u, U, F, x, p, cache::AbstractSystemCache)\n\nEvaluate the system F and its Jacobian at x and parameters p and store the results in u (evalution) and U (Jacobian).\n\n\n\n\n\n"
},

{
    "location": "systems/#Optional-1",
    "page": "Systems",
    "title": "Optional",
    "category": "section",
    "text": "The following methods are mandatory to implement. The following are optional to implement but usually you want to define at least cache.evaluate_and_jacobian!(u, U, F::AbstractSystem, x, cache::AbstractSystemCache)\nevaluate_and_jacobian!(u, U, F::AbstractSystem, x, p, cache::AbstractSystemCache)"
},

{
    "location": "homotopies/#",
    "page": "Homotopies",
    "title": "Homotopies",
    "category": "page",
    "text": ""
},

{
    "location": "homotopies/#Homotopies-1",
    "page": "Homotopies",
    "title": "Homotopies",
    "category": "section",
    "text": "A homotopy is a functionH mathbbC^N  mathbbC  mathbbC^n (xt)  H(xt)where H( t) is a polynomial system for all tmathbbC."
},

{
    "location": "homotopies/#HomotopyContinuation.StraightLineHomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.StraightLineHomotopy",
    "category": "type",
    "text": "StraightLineHomotopy(G, F; gamma=exp(i * 2π*rand()))\n\nConstruct the homotopy H(x t) = γtG(x) + (1-t)F(x).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.FixedPointHomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.FixedPointHomotopy",
    "category": "type",
    "text": "FixedPointHomotopy(F, x₀; gamma=exp(i * 2π*rand()))\n\nConstruct the homotopy H(x t) = (1-t)F(x) + γt(x-x₀).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.ParameterHomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.ParameterHomotopy",
    "category": "type",
    "text": "ParameterHomotopy(F, parameters;\n    variables=setdiff(MP.variables(F), parameters),\n    p₁=randn(ComplexF64, length(parameters)),\n    p₀=randn(ComplexF64, length(parameters)),\n    γ₁=nothing, γ₀=nothing)\n\nConstruct the homotopy\n\nH(x t) = F(x (tγ₁p₁+(1-t)γ₀p₀)  (tγ₁+(1-t)γ₀))\n\nwhere p₁ and p₀ are a vector of parameter values for F and γ₁ and γ₀ are complex numbers. If γ₁ or γ₀ is nothing, it is assumed that γ₁ and γ₀ are 1. The input parameters specifies the parameter variables of F. Neccessarily, length(parameters) == length(p₁) == length(p₀).\n\nNote that p₁ and p₀ are stored as a tuple p of SVectors and γ₁ and γ₀ are stored as a tuple γ or as γ=nothing\n\nParameterHomotopy(F, parameters;\n    variables=setdiff(MP.variables(F), parameters),\n    startparameters=randn(ComplexF64, length(parameters)),\n    targetparameters=randn(ComplexF64, length(parameters)),\n    startgamma=nothing, targetgamma=nothing)\n\nThis is a non-unicode variant where γ₁=startparameters, γ₀=targetparameters, γ₁=startgamma, γ₀=targetgamma.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.PatchedHomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.PatchedHomotopy",
    "category": "type",
    "text": "PatchedHomotopy(H::AbstractHomotopy, patch, v::PVector)\n\nAugment the homotopy H with the given patch v. This results in the system [H(x,t); v ⋅ x - 1]\n\n\n\n\n\n"
},

{
    "location": "homotopies/#Default-homotopies-1",
    "page": "Homotopies",
    "title": "Default homotopies",
    "category": "section",
    "text": "The following homotopies are available by defaultStraightLineHomotopy\nFixedPointHomotopy\nParameterHomotopyWe also provide more specialised homotopies, which are mostly used internally currently but could be useful in conjunction with the PathTracker primitive.PatchedHomotopy"
},

{
    "location": "homotopies/#Interface-for-custom-homotopies-1",
    "page": "Homotopies",
    "title": "Interface for custom homotopies",
    "category": "section",
    "text": "The great thing is that you are not limited to the homotopies provided by default. You can define your own homotopy by defining a struct with super type AbstractHomotopy. For this the following interface has to be defined."
},

{
    "location": "homotopies/#HomotopyContinuation.AbstractHomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.AbstractHomotopy",
    "category": "type",
    "text": "AbstractHomotopy\n\nRepresenting a homotopy.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.AbstractHomotopyCache",
    "page": "Homotopies",
    "title": "HomotopyContinuation.AbstractHomotopyCache",
    "category": "type",
    "text": "AbstractHomotopyCache\n\nA cache to avoid allocations for the evaluation of an AbstractHomotopy.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.HomotopyNullCache",
    "page": "Homotopies",
    "title": "HomotopyContinuation.HomotopyNullCache",
    "category": "type",
    "text": "HomotopyNullCache\n\nThe default AbstractHomotopyCache containing nothing.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#Types-1",
    "page": "Homotopies",
    "title": "Types",
    "category": "section",
    "text": "AbstractHomotopy\nAbstractHomotopyCache\nHomotopyNullCache"
},

{
    "location": "homotopies/#HomotopyContinuation.cache-Tuple{AbstractHomotopy,Any,Any}",
    "page": "Homotopies",
    "title": "HomotopyContinuation.cache",
    "category": "method",
    "text": "cache(H::AbstractHomotopy, x, t)::AbstractHomotopyCache\n\nCreate a cache for the evaluation (incl. Jacobian) of F with elements of the type of x. The default implementation returns HomotopyNullCache.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.evaluate!-Tuple{Any,AbstractHomotopy,Vararg{Any,N} where N}",
    "page": "Homotopies",
    "title": "HomotopyContinuation.evaluate!",
    "category": "method",
    "text": "evaluate!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H at (x, t) and store the result in u.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.jacobian!-Tuple{Any,AbstractHomotopy,Vararg{Any,N} where N}",
    "page": "Homotopies",
    "title": "HomotopyContinuation.jacobian!",
    "category": "method",
    "text": "jacobian!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the Jacobian of the homotopy H at (x, t) and store the result in u.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.dt!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.dt!",
    "category": "function",
    "text": "dt!(u, H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H at (x, t) and store the result in u.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#Base.size-Tuple{AbstractHomotopy}",
    "page": "Homotopies",
    "title": "Base.size",
    "category": "method",
    "text": "Base.size(H::AbstractHomotopy)\n\nReturns a tuple (m, n) indicating that H is a homotopy of m polynomials m in n variables.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#Mandatory-1",
    "page": "Homotopies",
    "title": "Mandatory",
    "category": "section",
    "text": "The following methods are mandatory to implement.cache(H::AbstractHomotopy, x, t)\nevaluate!(u, F::AbstractHomotopy, args...)\njacobian!(u, H::AbstractHomotopy, args...)\ndt!\nBase.size(::AbstractHomotopy)"
},

{
    "location": "homotopies/#HomotopyContinuation.evaluate_and_jacobian!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.evaluate_and_jacobian!",
    "category": "function",
    "text": "evaluate_and_jacobian!(u, U, F, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H and its Jacobian at (x, t) and store the results in u (evalution) and U (Jacobian).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.evaluate_and_jacobian",
    "page": "Homotopies",
    "title": "HomotopyContinuation.evaluate_and_jacobian",
    "category": "function",
    "text": "evaluate_and_jacobian(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H and its Jacobian at (x, t).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.jacobian_and_dt!",
    "page": "Homotopies",
    "title": "HomotopyContinuation.jacobian_and_dt!",
    "category": "function",
    "text": "jacobian_and_dt!(U, u, H, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H and its derivative w.r.t. t at (x, t) and store the results in u (evalution) and v (∂t).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.jacobian_and_dt",
    "page": "Homotopies",
    "title": "HomotopyContinuation.jacobian_and_dt",
    "category": "function",
    "text": "jacobian_and_dt(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H and its derivative w.r.t. t at (x, t).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.evaluate-Tuple{AbstractHomotopy,Any,Any}",
    "page": "Homotopies",
    "title": "HomotopyContinuation.evaluate",
    "category": "method",
    "text": "evaluate(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H at (x, t).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.jacobian-Tuple{AbstractHomotopy,Any,Any}",
    "page": "Homotopies",
    "title": "HomotopyContinuation.jacobian",
    "category": "method",
    "text": "jacobian(H::AbstractHomotopy, x, t, cache::AbstractHomotopyCache)\n\nEvaluate the Jacobian of the homotopy H at (x, t).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.dt-Tuple{AbstractHomotopy,Any,Any}",
    "page": "Homotopies",
    "title": "HomotopyContinuation.dt",
    "category": "method",
    "text": "dt(H::AbstractHomotopy, x::AbstractVector, cache::AbstractHomotopyCache)\n\nEvaluate the homotopy H at (x, t).\n\n\n\n\n\n"
},

{
    "location": "homotopies/#HomotopyContinuation.basehomotopy",
    "page": "Homotopies",
    "title": "HomotopyContinuation.basehomotopy",
    "category": "function",
    "text": "basehomotopy(H::AbstractHomotopy)\n\nReturns the \'proper\' homotopy describing the problem. Any wrapper homotopy recursively calls wrappedhomotopy with the wrapped homotopy as argument.\n\n\n\n\n\n"
},

{
    "location": "homotopies/#Optional-1",
    "page": "Homotopies",
    "title": "Optional",
    "category": "section",
    "text": "evaluate_and_jacobian!(u, U, H::AbstractHomotopy, x, t, c=cache(H, x, t))\nevaluate_and_jacobian(H::AbstractHomotopy, x, t, c=cache(H, x, t))\njacobian_and_dt!(U, u, H::AbstractHomotopy, x, t, c=cache(H, x, t))\njacobian_and_dt(H::AbstractHomotopy, x, t, c=cache(H, x, t))\nevaluate(H::AbstractHomotopy, x, t)\njacobian(H::AbstractHomotopy, x, t)\ndt(H::AbstractHomotopy, x, t)\nbasehomotopy"
},

{
    "location": "predictors-correctors/#",
    "page": "Predictors and Correctors",
    "title": "Predictors and Correctors",
    "category": "page",
    "text": ""
},

{
    "location": "predictors-correctors/#Predictors-and-Correctors-1",
    "page": "Predictors and Correctors",
    "title": "Predictors and Correctors",
    "category": "section",
    "text": "We use a predictor-corrector scheme to track paths. These are the predictors and correctors currently available."
},

{
    "location": "predictors-correctors/#HomotopyContinuation.Euler",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.Euler",
    "category": "type",
    "text": "Euler()\n\nThis uses the explicit Euler method for prediction, also known as the tangent predictor.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.Heun",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.Heun",
    "category": "type",
    "text": "Heun()\n\nThe Heun predictor of order 2.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.Ralston",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.Ralston",
    "category": "type",
    "text": "Ralston()\n\nThe Ralston predictor of order 2.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.RK3",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.RK3",
    "category": "type",
    "text": "RK3()\n\nThe classical Runge-Kutta predictor of order 3.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.RK4",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.RK4",
    "category": "type",
    "text": "RK4()\n\nThe classical Runge-Kutta predictor of order 4.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.Pade21",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.Pade21",
    "category": "type",
    "text": "Pade21()\n\nThis uses a Padé-approximation of type (2,1) for prediction.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.NullPredictor",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.NullPredictor",
    "category": "type",
    "text": "NullPredictor()\n\nA predictor which does no prediction step, i.e., it just returns the input as its prediction.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#Predictors-1",
    "page": "Predictors and Correctors",
    "title": "Predictors",
    "category": "section",
    "text": "The following predictors are currently implemented.Euler\nHeun\nRalston\nRK3\nRK4\nPade21\nNullPredictor"
},

{
    "location": "predictors-correctors/#HomotopyContinuation.NewtonCorrector",
    "page": "Predictors and Correctors",
    "title": "HomotopyContinuation.NewtonCorrector",
    "category": "type",
    "text": "NewtonCorrector(;simplified_last_step=true)\n\nAn ordinary Newton\'s method. If simplified_last_step is true, then for the last iteration the previously Jacobian will be used. This uses an LU-factorization for square systems and a QR-factorization for overdetermined.\n\n\n\n\n\n"
},

{
    "location": "predictors-correctors/#Correctors-1",
    "page": "Predictors and Correctors",
    "title": "Correctors",
    "category": "section",
    "text": "The following correctors are currently implemented.NewtonCorrector"
},

{
    "location": "pathtracking/#",
    "page": "Path tracker",
    "title": "Path tracker",
    "category": "page",
    "text": ""
},

{
    "location": "pathtracking/#HomotopyContinuation.pathtracker_startsolutions",
    "page": "Path tracker",
    "title": "HomotopyContinuation.pathtracker_startsolutions",
    "category": "function",
    "text": "pathtracker_startsolutions(args...; kwargs...)\n\nConstruct a PathTracker and startsolutions in the same way solve does it. This also takes the same input arguments as solve. This is convenient if you want to investigate single paths.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.pathtracker",
    "page": "Path tracker",
    "title": "HomotopyContinuation.pathtracker",
    "category": "function",
    "text": "pathtracker(args...; kwargs...)\n\nConstruct a PathTracker in the same way solve does it. This also takes the same input arguments as solve with the exception that you do not need to specify startsolutions. This is convenient if you want to investigate single paths.\n\nExamples\n\nObtain single solution\n\nWe want to construct a path tracker to track a parameterized system f with parameters p from the parameters a to b.\n\ntracker = pathtracker(f, parameters=p, p₁=a, p₀=b)\n\nYou then can obtain a single solution at b by using\n\nx_b = track(tracker, x_a).x\n\nTrace a path\n\nTo trace a path you can use the iterator method.\n\ntracker = pathtracker(f, parameters=p, p₁=a, p₀=b, maximal_step_size=0.01)\nfor (x, t) in iterator(tracker, x₁)\n    @show (x,t)\nend\n\nIf we want to guarantee smooth traces we can limit the maximal step size.\n\ntracker = pathtracker(f, parameters=p, p₁=a, p₀=b, maximal_step_size=0.01)\nfor (x, t) in iterator(tracker, x₁)\n    @show (x,t)\nend\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#Path-tracking-1",
    "page": "Path tracker",
    "title": "Path tracking",
    "category": "section",
    "text": "We also export a path tracking primitive to make the core path tracking routine available for other applications. At the heart is a PathTracker object which holds all the state. The easiest way to construct a PathTracker is to use the pathtracker_startsolutions routine.pathtracker_startsolutions\npathtracker"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTracker",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTracker",
    "category": "type",
    "text": " PathTracker(H::AbstractHomotopy, x₁, t₁, t₀; options...)::PathTracker\n\nCreate a PathTracker to track x₁ from t₁ to t₀. The homotopy H needs to be homogenous. Note that a PathTracker is also a (mutable) iterator.\n\nPathTrackerOptions\n\ncorrector::AbstractCorrector: The corrector used during in the predictor-corrector scheme. The default is NewtonCorrector.\ncorrector_maxiters=3: The maximal number of correction steps in a single step.\ninitial_step_size=0.1: The step size of the first step.\nmaxiters=10_000: The maximal number of iterations the path tracker has available.\nminimal_step_size=1e-14: The minimal step size.\nmaximal_step_size=Inf: The maximal step size.\npredictor::AbstractPredictor: The predictor used during in the predictor-corrector scheme. The default is Heun()`.\nrefinement_maxiters=corrector_maxiters: The maximal number of correction steps used to refine the final value.\nrefinement_tol=1e-8: The precision used to refine the final value.\ntol=1e-7: The precision used to track a value.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTrackerResult",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTrackerResult",
    "category": "type",
    "text": " PathTrackerResult(tracker)\n\nContaining the result of a tracked path. The fields are\n\nsuccessfull::Bool Indicating whether tracking was successfull.\nreturncode::PathTrackerStatus.states If the tracking was successfull then it is PathTrackerStatus.success.\nx::V The result.\nt::Float64 The t when the path tracker stopped.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.PathTrackerStatus.states",
    "page": "Path tracker",
    "title": "HomotopyContinuation.PathTrackerStatus.states",
    "category": "type",
    "text": "PathTrackerStatus.states\n\nThe possible states the pathtracker can achieve are\n\nPathTrackerStatus.success\nPathTrackerStatus.tracking\nPathTrackerStatus.terminated_maximal_iterations\nPathTrackerStatus.terminated_invalid_startvalue\nPathTrackerStatus.terminated_step_size_too_small\nPathTrackerStatus.terminated_singularity\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#Types-1",
    "page": "Path tracker",
    "title": "Types",
    "category": "section",
    "text": "PathTracker\nPathTrackerResult\nPathTrackerStatus.states"
},

{
    "location": "pathtracking/#HomotopyContinuation.track",
    "page": "Path tracker",
    "title": "HomotopyContinuation.track",
    "category": "function",
    "text": "track(tracker, x₁, t₁=1.0, t₀=0.0; options...)::PathTrackerResult\n\nTrack a value x₁ from t₁ to t₀ using the given PathTracker tracker. This returns a PathTrackerResult. This modifies tracker. See track! for the possible options.\n\n\n\n\n\ntrack(tracker, x::AbstractVector, edge::Edge, loop::Loop, stats::MonodromyStatistics)\n\nTrack x along the edge edge in the loop loop using tracker. Record statistics in stats.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.track!",
    "page": "Path tracker",
    "title": "HomotopyContinuation.track!",
    "category": "function",
    "text": " track!(tracker, x₁, t₁=1.0, t₀=0.0; setup_patch=true, checkstartvalue=true, compute_ẋ=true)\n\nTrack a value x₁ from t₁ to t₀ using the given PathTracker tracker. Returns one of the enum values of PathTrackerStatus.states indicating the status. If the tracking was successfull it is PathTrackerStatus.success. If setup_patch is true then setup! is called at the beginning of the tracking.\n\ntrack!(x₀, tracker, x₁, t₁=1.0, t₀=0.0; options...)\n\nAdditionally also stores the result in x₀ if the tracking was successfull.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.setup!",
    "page": "Path tracker",
    "title": "HomotopyContinuation.setup!",
    "category": "function",
    "text": "setup!(::AbstractAffinePatchState, x::AbstractVector)\n\nSetup the affine patch depending on x and modify x if necessary. This is only called once at the beginning of a tracked path.\n\n\n\n\n\nsetup!(cache::AbstractStatefulPredictorCache, H, x, ẋ, t, fac)\n\nSetup the cache. x is the new path value at t and ẋ is the derivative at t. fac is a factorization of the Jacobian at (x,t). This falls back to calling update.\n\n\n\n\n\nsetup!(pathtracker, x₁, t₁=1.0, t₀=0.0, setup_patch=pathtracker.options.update_patch, checkstartvalue=true, compute_ẋ=true)\n\nSetup pathtracker to track x₁ from t₁ to t₀. Use this if you want to use the pathtracker as an iterator.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.iterator",
    "page": "Path tracker",
    "title": "HomotopyContinuation.iterator",
    "category": "function",
    "text": "iterator(tracker::PathTracker, x₁, t₁=1.0, t₀=0.0; affine=true)\n\nPrepare a tracker to make it usable as a (stateful) iterator. Use this if you want to inspect a specific path. In each iteration the tuple (x,t) is returned. If affine == true then x is the affine solution (internally we compute in projective space).\n\nExample\n\nAssume you have PathTracker tracker and you wan to track x₁ from 1.0 to 0.25:\n\nfor (x,t) in iterator(tracker, x₁, 1.0, 0.25)\n    println(\"x at t=$t:\")\n    println(x)\nend\n\nNote that this is a stateful iterator. You can still introspect the state of the tracker. For example to check whether the tracker was successfull (and did not terminate early due to some problem) you can do\n\nprintln(\"Success: \", currstatus(tracker) == PathTrackerStatus.success)\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#Methods-1",
    "page": "Path tracker",
    "title": "Methods",
    "category": "section",
    "text": "To track from a start to an endpoint with the PathTracker we provide the following routines.track\ntrack!\nsetup!It is also possible to use a PathTracker as an iterator. This can either be done by the high level iterator method or by directly using a PathTracker as an iterator. The recommend approach is simply using iterator.iterator"
},

{
    "location": "pathtracking/#HomotopyContinuation.currx",
    "page": "Path tracker",
    "title": "HomotopyContinuation.currx",
    "category": "function",
    "text": "currx(tracker::PathTracker)\n\nReturn the current value of x.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.currt",
    "page": "Path tracker",
    "title": "HomotopyContinuation.currt",
    "category": "function",
    "text": " currt(tracker::PathTracker)\n\nCurrent t.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.currΔt",
    "page": "Path tracker",
    "title": "HomotopyContinuation.currΔt",
    "category": "function",
    "text": " currΔt(tracker::PathTracker)\n\nCurrent step_size Δt.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.curriters",
    "page": "Path tracker",
    "title": "HomotopyContinuation.curriters",
    "category": "function",
    "text": " curriters(tracker::PathTracker)\n\nCurrent number of iterations.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.currstatus",
    "page": "Path tracker",
    "title": "HomotopyContinuation.currstatus",
    "category": "function",
    "text": " currstatus(tracker::PathTracker)\n\nCurrent status.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#Introspecting-the-current-state-1",
    "page": "Path tracker",
    "title": "Introspecting the current state",
    "category": "section",
    "text": "To introspect the current state we provide the following routines.currx\ncurrt\ncurrΔt\ncurriters\ncurrstatus"
},

{
    "location": "pathtracking/#HomotopyContinuation.corrector_maxiters",
    "page": "Path tracker",
    "title": "HomotopyContinuation.corrector_maxiters",
    "category": "function",
    "text": " corrector_maxiters(tracker::PathTracker)\n\nCurrent correction maxiters.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.set_corrector_maxiters!",
    "page": "Path tracker",
    "title": "HomotopyContinuation.set_corrector_maxiters!",
    "category": "function",
    "text": " set_corrector_maxiters!(tracker::PathTracker, n)\n\nSet the correction maxiters to n.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.maximal_step_size",
    "page": "Path tracker",
    "title": "HomotopyContinuation.maximal_step_size",
    "category": "function",
    "text": " maximal_step_size(tracker::PathTracker)\n\nCurrent maximal step size.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.set_maximal_step_size!",
    "page": "Path tracker",
    "title": "HomotopyContinuation.set_maximal_step_size!",
    "category": "function",
    "text": " set_corrector_maxiters!(tracker::PathTracker, Δs)\n\nSet the maximal step size to Δs.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.refinement_maxiters",
    "page": "Path tracker",
    "title": "HomotopyContinuation.refinement_maxiters",
    "category": "function",
    "text": " refinement_maxiters(tracker::PathTracker)\n\nCurrent refinement maxiters.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.set_refinement_maxiters!",
    "page": "Path tracker",
    "title": "HomotopyContinuation.set_refinement_maxiters!",
    "category": "function",
    "text": " set_refinement_maxiters!(tracker::PathTracker, n)\n\nSet the current refinement maxiters to n.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.refinement_tol",
    "page": "Path tracker",
    "title": "HomotopyContinuation.refinement_tol",
    "category": "function",
    "text": " refinement_tol(tracker::PathTracker)\n\nCurrent refinement tolerance.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.set_refinement_tol!",
    "page": "Path tracker",
    "title": "HomotopyContinuation.set_refinement_tol!",
    "category": "function",
    "text": " set_refinement_maxiters!(tracker::PathTracker, tol)\n\nSet the current refinement tolerance to tol.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.tol",
    "page": "Path tracker",
    "title": "HomotopyContinuation.tol",
    "category": "function",
    "text": " tol(tracker::PathTracker)\n\nCurrent tolerance.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#HomotopyContinuation.set_tol!",
    "page": "Path tracker",
    "title": "HomotopyContinuation.set_tol!",
    "category": "function",
    "text": " set_tol!(tracker::PathTracker, tol)\n\nSet the current tolerance to tol.\n\n\n\n\n\n"
},

{
    "location": "pathtracking/#Changing-options-1",
    "page": "Path tracker",
    "title": "Changing options",
    "category": "section",
    "text": "To change settingscorrector_maxiters\nset_corrector_maxiters!\nmaximal_step_size\nset_maximal_step_size!\nrefinement_maxiters\nset_refinement_maxiters!\nrefinement_tol\nset_refinement_tol!\ntol\nset_tol!"
},

{
    "location": "newton/#",
    "page": "Newton\'s method",
    "title": "Newton\'s method",
    "category": "page",
    "text": ""
},

{
    "location": "newton/#HomotopyContinuation.newton",
    "page": "Newton\'s method",
    "title": "HomotopyContinuation.newton",
    "category": "function",
    "text": "newton(F::AbstractSystem, x₀; tol=1e-6, maxiters=3, simplified_last_step=true)\n\nAn ordinary Newton\'s method. If simplified_last_step is true, then for the last iteration the previously Jacobian will be used. This uses an LU-factorization for square systems and a QR-factorization for overdetermined.\n\n\n\n\n\n"
},

{
    "location": "newton/#HomotopyContinuation.NewtonResult",
    "page": "Newton\'s method",
    "title": "HomotopyContinuation.NewtonResult",
    "category": "type",
    "text": "NewtonResult{T}\n\nStructure holding information about the outcome of the newton function. The fields are.\n\nretcode The return code of the compuation. converged means that accuracy ≤ tol.\naccuracy::T |xᵢ-xᵢ₋₁| for i = iters and x₀,x₁,…,xᵢ₋₁,xᵢ are the Newton iterates.\niters::Int The number of iterations used.\nnewton_update_error::Float64 δx/(|x| ⋅ ϵ), where x is the first Newton update, δx is an estimate for the error |x-x̂| between x and the exact Newton update x̂ and ϵ is the machine precision.\n\n\n\n\n\n"
},

{
    "location": "newton/#HomotopyContinuation.newton!",
    "page": "Newton\'s method",
    "title": "HomotopyContinuation.newton!",
    "category": "function",
    "text": "newton!(out, F::AbstractSystem, x₀, tol, maxiters::Integer, simplified_last_step::Bool,  cache::NewtonCache, newton_update_error = 1.0)\n\nIn-place version of newton. Needs a NewtonCache as input.\n\n\n\n\n\n"
},

{
    "location": "newton/#HomotopyContinuation.NewtonCache",
    "page": "Newton\'s method",
    "title": "HomotopyContinuation.NewtonCache",
    "category": "type",
    "text": "NewtonCache(F::AbstractSystem, x)\n\nCache for the newton function.\n\n\n\n\n\n"
},

{
    "location": "newton/#Newton\'s-method-1",
    "page": "Newton\'s method",
    "title": "Newton\'s method",
    "category": "section",
    "text": "Sometimes it is necessary to refine obtained solutions. For this we provide an interface to Newton\'s method.newton\nNewtonResultFor high performance applications we also provide an in-place version of Newton\'s method which avoids any temporary allocations.newton!\nNewtonCache"
},

{
    "location": "sorting/#",
    "page": "Sorting arrays of solutions",
    "title": "Sorting arrays of solutions",
    "category": "page",
    "text": ""
},

{
    "location": "sorting/#Sorting-arrays-of-solutions-1",
    "page": "Sorting arrays of solutions",
    "title": "Sorting arrays of solutions",
    "category": "section",
    "text": "We provide functions for sorting analyzing arrays of vectors."
},

{
    "location": "sorting/#HomotopyContinuation.UniquePoints",
    "page": "Sorting arrays of solutions",
    "title": "HomotopyContinuation.UniquePoints",
    "category": "type",
    "text": "UniquePoints{V<:AbstractVector, T, F<:Function}\n\nA data structure which holds points of type V where T=real(eltype(V)). This data structure provides an efficient (poly)logarithmic check whether a point already exists where two points u,v are considered equal if F(u,v)<tol, where tol is a tolerance provided through the add! function.\n\nUniquePoints(v::AbstractVector{<:Number}, distance::F)\n\nInitialize the data structure with just one data point v.\n\nUniquePoints(V::Vector{<:AbstractVector{<:Number}}, distance::F; tol=1e-5)\n\nInitialize the data structure with all points in v. These are added in order by add! with the given tolerance tol. In particular, \'UniquePoints\' structure will contain only points for which the pairwise distance given by F is less than tol.\n\nUniquePoints(v) = UniquePoints(v, euclidean_distance)\n\nIf F is not specialized, euclidean_distance is used.\n\nExample\n\njulia> UniquePoints([[1,0.5]; [1,0.5]; [1,1]])\n[[1,0.5], [1,1]]\n\nThis is the same as\n\nUniquePoints([[1,0.5]; [1,0.5]; [1,1]], (x,y) -> LinearAlgebra.norm(x-y))\n\n\n\n\n\n"
},

{
    "location": "sorting/#HomotopyContinuation.points",
    "page": "Sorting arrays of solutions",
    "title": "HomotopyContinuation.points",
    "category": "function",
    "text": "points(data::UniquePoints)\n\nReturn the points stored in data.\n\n\n\n\n\n"
},

{
    "location": "sorting/#HomotopyContinuation.iscontained",
    "page": "Sorting arrays of solutions",
    "title": "HomotopyContinuation.iscontained",
    "category": "function",
    "text": "iscontained(data::UniquePoints{V}, x::V; tol=1e-5)::Bool\n\nCheck whether x is contained in the data by using the tolerance tol to decide for duplicates.\n\niscontained(data::UniquePoints{V}, x::V, Val{true}(); tol=1e-5)::Int\n\nIf x is contained in data by using the tolerance tol return the index of the data point which already exists. If the data point is not existing -1 is returned.\n\n\n\n\n\niscontained(node::Node, x; kwargs...)\n\nCalls iscontained on the points of the Node.\n\n\n\n\n\n"
},

{
    "location": "sorting/#HomotopyContinuation.add!",
    "page": "Sorting arrays of solutions",
    "title": "HomotopyContinuation.add!",
    "category": "function",
    "text": "add!(data::UniquePoints{V}, x::V; tol=1e-5)::Bool\n\nAdd x to data if it doesn\'t already exists by using the tolerance tol to decide for duplicates.\n\nadd!(data::UniquePoints{V}, x::V, Val(true); tol=1e-5)::Int\n\nIf x is contained in data by using the tolerance tol to decide for duplicates return the index of the data point which already exists. If the data point is not existing add it to x and return -1. The element will be the last element of points(data).\n\n\n\n\n\nadd!(node::Node, x; kwargs...)\n\nCalls add! on the points of the Node.\n\n\n\n\n\n"
},

{
    "location": "sorting/#HomotopyContinuation.unsafe_add!",
    "page": "Sorting arrays of solutions",
    "title": "HomotopyContinuation.unsafe_add!",
    "category": "function",
    "text": "unsafe_add!(data::UniquePoints{V}, x::V)::Bool\n\nSimilarly to add! but assumes that it was already checked that there is no duplicate with iscontained. This has to be called directly after iscontained with the same value of x.\n\n\n\n\n\n"
},

{
    "location": "sorting/#Base.empty!",
    "page": "Sorting arrays of solutions",
    "title": "Base.empty!",
    "category": "function",
    "text": "empty!(collection) -> collection\n\nRemove all elements from a collection.\n\nExamples\n\njulia> A = Dict(\"a\" => 1, \"b\" => 2)\nDict{String,Int64} with 2 entries:\n  \"b\" => 2\n  \"a\" => 1\n\njulia> empty!(A);\n\njulia> A\nDict{String,Int64} with 0 entries\n\n\n\n\n\nempty!(data::UniquePoints)\n\nRemove all points from data.\n\n\n\n\n\n"
},

{
    "location": "sorting/#Computing-unique-points-in-an-array-of-vectors-1",
    "page": "Sorting arrays of solutions",
    "title": "Computing unique points in an array of vectors",
    "category": "section",
    "text": "UniquePointsWe provide several helper functions for UniquePoints.points\niscontained\nadd!\nunsafe_add!\nempty!"
},

{
    "location": "sorting/#HomotopyContinuation.multiplicities",
    "page": "Sorting arrays of solutions",
    "title": "HomotopyContinuation.multiplicities",
    "category": "function",
    "text": "multiplicities(vectors, distance=euclidean_distance; tol::Real = 1e-5)\n\nReturns an array of arrays of integers. Each vector w in \'v\' contains all indices i,j such that w[i] and w[j] have distance at most tol.\n\nmultiplicities(v; tol::Real = 1e-5) = multiplicities(v, euclidean_distance, tol = tol)\n\nIf distance is not specified, euclidean_distance is used.\n\njulia> multiplicities([[1,0.5]; [1,0.5]; [1,1]])\n[[1,2]]\n\nThis is the same as\n\nmultiplicities([[1,0.5]; [1,0.5]; [1,1]], (x,y) -> LinearAlgebra.norm(x-y))\n\n\n\n\n\nmultiplicities(V::Results; tol=1e-6)\n\nReturns a Vector of Vector{PathResult}s grouping the PathResults whose solutions appear with multiplicities greater 1 in \'V\'. Two solutions are regarded as equal, when their pairwise distance is less than \'tol\'.\n\n\n\n\n\n"
},

{
    "location": "sorting/#Computing-points-in-an-array-of-vectors-which-appear-multiple-times-1",
    "page": "Sorting arrays of solutions",
    "title": "Computing points in an array of vectors which appear multiple times",
    "category": "section",
    "text": "If instead of unique points, the user wants to have the information which points in an array of points appear with multiplicity, they should use the next function.multiplicitiesThe multiplicities functions may also be applied to AffineResult and ProjectiveResult structures; see here: multiplicities(::HomotopyContinuation.Results)."
},

{
    "location": "norms_distances/#",
    "page": "Norms and Distances arrays of solutions",
    "title": "Norms and Distances arrays of solutions",
    "category": "page",
    "text": ""
},

{
    "location": "norms_distances/#HomotopyContinuation.euclidean_distance",
    "page": "Norms and Distances arrays of solutions",
    "title": "HomotopyContinuation.euclidean_distance",
    "category": "function",
    "text": "euclidean_distance(u, v)\n\nCompute ||u-v||₂.\n\n\n\n\n\n"
},

{
    "location": "norms_distances/#HomotopyContinuation.euclidean_norm",
    "page": "Norms and Distances arrays of solutions",
    "title": "HomotopyContinuation.euclidean_norm",
    "category": "function",
    "text": "euclidean_norm(u)\n\nCompute ||u||₂.\n\n\n\n\n\n"
},

{
    "location": "norms_distances/#HomotopyContinuation.infinity_distance",
    "page": "Norms and Distances arrays of solutions",
    "title": "HomotopyContinuation.infinity_distance",
    "category": "function",
    "text": "infinity_distance(u, v)\n\nCompute the ∞-norm of u-v.\n\n\n\n\n\n"
},

{
    "location": "norms_distances/#HomotopyContinuation.infinity_norm",
    "page": "Norms and Distances arrays of solutions",
    "title": "HomotopyContinuation.infinity_norm",
    "category": "function",
    "text": "infinity_norm(z)\n\nCompute the ∞-norm of z. If z is a complex vector this is more efficient than norm(z, Inf).\n\ninfinity_norm(z₁, z₂)\n\nCompute the ∞-norm of z₁-z₂.\n\n\n\n\n\n"
},

{
    "location": "norms_distances/#HomotopyContinuation.fubini_study",
    "page": "Norms and Distances arrays of solutions",
    "title": "HomotopyContinuation.fubini_study",
    "category": "function",
    "text": "fubini_study(x, y)\n\nComputes the Fubini-Study distance between x and y.\n\n\n\n\n\n"
},

{
    "location": "norms_distances/#Distances-and-norms-1",
    "page": "Norms and Distances arrays of solutions",
    "title": "Distances and norms",
    "category": "section",
    "text": "We provide functions for computing norms and distances.euclidean_distance\neuclidean_norminfinity_distance\ninfinity_normfubini_study"
},

{
    "location": "reference/#",
    "page": "Reference",
    "title": "Reference",
    "category": "page",
    "text": ""
},

{
    "location": "reference/#Reference-1",
    "page": "Reference",
    "title": "Reference",
    "category": "section",
    "text": ""
},

{
    "location": "reference/#Input-1",
    "page": "Reference",
    "title": "Input",
    "category": "section",
    "text": "We support any polynomials which follow the MultivariatePolynomials interface. By default we export the routines @polyvar, PolyVar, differentiate and variables from the DynamicPolynomials implementation. With these you can simply create variables# Create variables x, y, z\n@polyvar x y z\nf = x^2+y^2+z^2\n\n# You can also create an array of variables\n@polyvar x[1:3] # This creates x1, x2, x3 accessed by x[1], x[2], x[3]\nf = dot(x, x) # = x[1]^2+x[2]^2+x[3]^2\n\n# Also you can create matrices of variables\n# This creates x1_1, x1_2, x2_1, x2_2 accessed by\n# x[1,1], x[1,2], x[2,1], x[2,2]\n@polyvar x[1:2, 1:2]"
},

{
    "location": "reference/#HomotopyContinuation.bezout_number",
    "page": "Reference",
    "title": "HomotopyContinuation.bezout_number",
    "category": "function",
    "text": "bezout_number(F::MPPolys; variable_groups=[variables(F)], homvars=nothing, parameters=nothing)\nbezout_number(multidegrees, groups::VariableGroups)\n\nCompute the multi-homogenous bezout number associated to the given system and variable groups.\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.ishomogenous",
    "page": "Reference",
    "title": "HomotopyContinuation.ishomogenous",
    "category": "function",
    "text": "ishomogenous(f::MP.AbstractPolynomialLike)\n\nChecks whether f is homogenous.\n\nishomogenous(f::MP.AbstractPolynomialLike, vars)\n\nChecks whether f is homogenous in the variables vars with possible weights.\n\n\n\n\n\nishomogenous(F::Vector{MP.AbstractPolynomialLike}, variables)\n\nChecks whether each polynomial in F is homogenous in the variables variables.\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.uniquevar",
    "page": "Reference",
    "title": "HomotopyContinuation.uniquevar",
    "category": "function",
    "text": "uniquevar(f::MP.AbstractPolynomialLike, tag=:x0)\nuniquevar(F::MPPolys, tag=:x0)\n\nCreates a unique variable.\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.homogenize",
    "page": "Reference",
    "title": "HomotopyContinuation.homogenize",
    "category": "function",
    "text": "homogenize(f::MP.AbstractPolynomial, variable=uniquevar(f))\n\nHomogenize the polynomial f by using the given variable variable.\n\nhomogenize(F::Vector{<:MP.AbstractPolynomial}, variable=uniquevar(F))\n\nHomogenize each polynomial in F by using the given variable variable.\n\nhomogenize(f::MP.AbstractPolynomial, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(f))\n\nHomogenize the variables v in the polynomial f by using the given variable variable.\n\nhomogenize(F::Vector{<:MP.AbstractPolynomial}, v::Vector{<:MP.AbstractVariable}, variable=uniquevar(F))\n\nHomogenize the variables v in each polynomial in F by using the given variable variable.\n\n\n\n\n\n"
},

{
    "location": "reference/#Utilities-1",
    "page": "Reference",
    "title": "Utilities",
    "category": "section",
    "text": "bezout_number\nishomogenous\nuniquevar\nhomogenize"
},

{
    "location": "reference/#HomotopyContinuation.OrthogonalPatch",
    "page": "Reference",
    "title": "HomotopyContinuation.OrthogonalPatch",
    "category": "type",
    "text": "OrthogonalPatch()\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.EmbeddingPatch",
    "page": "Reference",
    "title": "HomotopyContinuation.EmbeddingPatch",
    "category": "type",
    "text": "EmbeddingPatch()\n\nHolds an PVector onto its affine patch. With this the effect is basically the same as tracking in affine space.\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.RandomPatch",
    "page": "Reference",
    "title": "HomotopyContinuation.RandomPatch",
    "category": "type",
    "text": "RandomPatch()\n\nA random patch. The vector has norm 1.\n\n\n\n\n\n"
},

{
    "location": "reference/#HomotopyContinuation.FixedPatch",
    "page": "Reference",
    "title": "HomotopyContinuation.FixedPatch",
    "category": "type",
    "text": "FixedPatch()\n\n\n\n\n\n"
},

{
    "location": "reference/#AffinePatches-1",
    "page": "Reference",
    "title": "AffinePatches",
    "category": "section",
    "text": "Affine patches are there to augment projective system such that they can be considered as (locally) affine system. By default the following patches are definedOrthogonalPatch\nEmbeddingPatch\nRandomPatch\nFixedPatch"
},

]}
